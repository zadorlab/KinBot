from kinbot.species_routing import routing_name
import os
import logging
import numpy as np
from ase.data import atomic_numbers
import re
from copy import deepcopy
from dataclasses import dataclass
import hashlib
import shutil
import shlex
from pathlib import Path

from kinbot import kb_path
from kinbot.theory import fortran_threshold
from kinbot.execution import file_sha256
from kinbot import constants
import math

logger = logging.getLogger('KinBot')


def validate_molpro_convergence(output, *, allow_scf_guess_failure=False):
    """Reject explicit native solver failures even after normal termination.

    For multireference energies, a failed HF guess is usable only when its
    subsequent MULTI optimization explicitly converged. Check every program
    invocation, including displaced geometries in a numerical force job.
    """
    failure = r'NO CONVERGENCE\s*!|No convergence in (?:[ru]?hf(?:pr)?|multi|mcscf|[ru]?ccsd)\b|(?:[RU]?HF|MCSCF|[RU]?CCSD)[^\n]*not converged'
    programs = list(re.finditer(r'PROGRAM\s*\*\s*(\S+)', output, re.I))
    sections = [output[:programs[0].start()]] if programs else [output]
    for index, match in enumerate(programs):
        end = programs[index + 1].start() if index + 1 < len(programs) else len(output)
        block = output[match.start():end]
        if allow_scf_guess_failure and re.match(r'[RU]?HF\b', match[1], re.I):
            for following in range(index + 1, len(programs)):
                name = programs[following][1]
                if re.match(r'[RU]?HF\b', name, re.I):
                    break
                if name.upper() == 'MULTI':
                    stop = programs[following + 1].start() if following + 1 < len(programs) else len(output)
                    multi = output[programs[following].end():stop]
                    if (re.search(r'(?im)^\s*CONVERGENCE REACHED!', multi)
                            and not re.search(failure, multi, re.I)):
                        block = ''
                    break
        sections.append(block)
    if re.search(failure, '\n'.join(sections), re.I):
        raise ValueError('Molpro output contains an unconverged calculation.')


def molpro_wfu_name(stem):
    """Return a lowercase checkpoint basename within Molpro's 31-char limit.

    Keep full identities in metadata. FILE names are lowercased by Molpro and
    long names have been observed silently truncated on the cluster.
    """
    name = f'{stem}.wfu'
    if len(name) <= 31 and re.fullmatch(r'[a-z0-9][a-z0-9_.-]*', name):
        return name
    return 'kb_' + hashlib.sha256(name.encode()).hexdigest()[:20] + '.wfu'


def compact_molpro_titles(contents):
    """Keep execution titles short enough for Molpro's native XML exporter.

    Long generated job titles can overflow Molpro's internal XML strings even
    after a successful electronic calculation. The 32-character bound leaves
    conservative headroom; it is not a documented Molpro input-title limit.
    Only title lines change. Full job identities remain in paths and metadata.
    """
    def compact(match):
        title = match[2]
        if len(title) <= 32 and title.isascii():
            return match[0]
        digest = hashlib.sha256(title.encode('utf-8')).hexdigest()[:16]
        return match[1] + 'KinBot_' + digest

    return re.sub(r'(?m)^([ \t]*\*\*\*[ \t]*,)([^\r\n]*)', compact, contents)


def validate_input_blocks(contents, label='Molpro input'):
    """Reject unmatched generated program/geometry delimiters before launch.

    This is a structural check, not a replacement for Molpro's input parser.
    Quoted paths and comments are ignored while counting block delimiters.
    """
    depth = 0
    for number, line in enumerate(contents.splitlines(), 1):
        quote = None
        for character in line:
            if quote:
                if character == quote:
                    quote = None
                continue
            if character in ('"', "'"):
                quote = character
            elif character == '!':
                break
            elif character == '{':
                depth += 1
            elif character == '}':
                depth -= 1
                if depth < 0:
                    raise ValueError(f'{label}: unmatched closing brace at line {number}.')
    if depth:
        raise ValueError(f'{label}: {depth} unclosed Molpro input block(s).')


_VRC_MULTIREFERENCE = re.compile(
    r'(caspt2|mrci)(?:-f12([ab])?)?(\+q)?\((\d+),(\d+)\)')

_AVAS_ORBITAL = re.compile(
    r'(?:(?P<principal>[1-9]\d*))?'
    r'(?P<kind>s|p|d|f|g|h|px|py|pz|d0|d1[+-]|d2[+-])\Z',
    re.IGNORECASE)


# Main-group valence shells used by the graph-derived VRC active-space
# builder.  The table describes elements, rather than named molecules or
# reaction channels.  Hydrogen orbitals are not AVAS projection centres when
# a heavy atom is present: the resulting molecular valence orbitals still
# contain their H 1s contributions, while the projection selects the compact
# heavy-atom valence manifold used in dissociation work.
_MAIN_GROUP_VALENCE = {
    'H': (1, ('1s',)), 'He': (2, ('1s',)),
    'Li': (1, ('2s',)), 'Be': (2, ('2s',)),
    'B': (3, ('2s', '2p')), 'C': (4, ('2s', '2p')),
    'N': (5, ('2s', '2p')), 'O': (6, ('2s', '2p')),
    'F': (7, ('2s', '2p')), 'Ne': (8, ('2s', '2p')),
    'Na': (1, ('3s',)), 'Mg': (2, ('3s',)),
    'Al': (3, ('3s', '3p')), 'Si': (4, ('3s', '3p')),
    'P': (5, ('3s', '3p')), 'S': (6, ('3s', '3p')),
    'Cl': (7, ('3s', '3p')), 'Ar': (8, ('3s', '3p')),
    'K': (1, ('4s',)), 'Ca': (2, ('4s',)),
    'Ga': (3, ('4s', '4p')), 'Ge': (4, ('4s', '4p')),
    'As': (5, ('4s', '4p')), 'Se': (6, ('4s', '4p')),
    'Br': (7, ('4s', '4p')), 'Kr': (8, ('4s', '4p')),
}


def full_valence_active_space(symbols, charge=0, *, additional_centers=None,
                              include_light_atoms=False):
    """Derive one whole-molecule main-group valence CAS and AVAS target.

    The electron count contains every main-group valence electron and the
    molecular charge.  Projection centres are all heavy atoms, or all atoms
    for an all-hydrogen system.  This is intentionally independent of which
    bond is later selected as the VRC coordinate.
    """
    symbols = [str(symbol) for symbol in symbols]
    if not symbols:
        raise ValueError('A full-valence active space requires atoms.')
    unsupported = sorted(set(symbols) - set(_MAIN_GROUP_VALENCE))
    if unsupported:
        raise ValueError(
            'Automatic full-valence VRC spaces support main-group elements; '
            'provide an explicit active-space model for: '
            + ', '.join(unsupported))
    if isinstance(charge, bool) or not isinstance(charge, (int, np.integer)):
        raise ValueError('Molecular charge must be an integer.')
    electrons = sum(_MAIN_GROUP_VALENCE[symbol][0]
                    for symbol in symbols) - int(charge)
    selected = [index for index, symbol in enumerate(symbols, 1)
                if include_light_atoms or symbol not in ('H', 'He')]
    if not selected:
        selected = list(range(1, len(symbols) + 1))
    centers = []
    orbitals = 0
    for index in selected:
        labels = list(_MAIN_GROUP_VALENCE[symbols[index - 1]][1])
        centers.append({'atom': index, 'orbitals': labels})
        orbitals += sum(_avas_orbital_count(label) for label in labels)
    for addition in additional_centers or []:
        atom = int(addition['atom'])
        if atom < 1 or atom > len(symbols):
            raise ValueError('An additional AVAS atom is outside the molecule.')
        target = next((center for center in centers
                       if center['atom'] == atom), None)
        if target is None:
            target = {'atom': atom, 'orbitals': []}
            centers.append(target)
        for label in addition['orbitals']:
            token = str(label).strip().casefold()
            _avas_orbital_count(token)
            if token in target['orbitals']:
                raise ValueError('An additional AVAS orbital duplicates the '
                                 'whole-molecule valence space.')
            target['orbitals'].append(token)
            orbitals += _avas_orbital_count(token)
    if electrons < 1 or electrons > 2 * orbitals:
        raise ValueError(
            f'The derived full-valence CAS({electrons},{orbitals}) is not '
            'compatible with the selected heavy-atom valence manifold.')
    return {
        'active_electrons': int(electrons),
        'active_orbitals': int(orbitals),
        'orbital_selection': {
            'mode': 'avas', 'centers': centers,
            'threshold': 0.3, 'minao_basis': 'minao',
        },
        'selection': 'whole_molecule_main_group_valence',
        'additional_orbitals': deepcopy(additional_centers or []),
    }


def method_with_active_space(method, active_space):
    """Replace a VRC multireference method's CAS by a derived space."""
    source = str(method).strip()
    match = _VRC_MULTIREFERENCE.fullmatch(source.casefold())
    if match is None:
        return source
    prefix = source[:source.rfind('(')]
    return (f'{prefix}({int(active_space["active_electrons"])},'
            f'{int(active_space["active_orbitals"])})')


def automatic_vrc_state_model(symbols, charge, multiplicities, *,
                              roots=None, selected_roots=None,
                              active_space=None):
    """Build one C1 SA-CASSCF model from the molecular valence manifold."""
    space = (full_valence_active_space(symbols, charge)
             if active_space is None else deepcopy(active_space))
    multiplicities = list(dict.fromkeys(int(value)
                                        for value in multiplicities))
    if not multiplicities or any(value < 1 for value in multiplicities):
        raise ValueError('Automatic state multiplicities must be positive.')
    roots = {} if roots is None else {int(key): int(value)
                                      for key, value in roots.items()}
    averages = []
    for multiplicity in multiplicities:
        count = roots.get(multiplicity, 1)
        averages.append({'multiplicity': multiplicity, 'roots': count,
                         'weights': [1. / count] * count})
    model = {
        'state_average': averages,
        'selected_roots': ({} if selected_roots is None
                           else dict(selected_roots)),
        'canonical_orbitals': True,
    }
    if space.get('orbital_selection') is not None:
        model['orbital_selection'] = space['orbital_selection']
    return validate_vrc_state_model(
        model, active_electrons=space['active_electrons'],
        active_orbitals=space['active_orbitals'])


def _avas_orbital_count(label):
    """Return the number of atomic functions requested by one AVAS token."""
    match = _AVAS_ORBITAL.fullmatch(str(label).strip())
    if match is None:
        raise ValueError(f'Invalid Molpro AVAS orbital specification: {label!r}.')
    kind = match.group('kind').casefold()
    return {'s': 1, 'p': 3, 'd': 5, 'f': 7, 'g': 9, 'h': 11}.get(
        kind, 1)


def validate_vrc_state_model(model, *, active_electrons=None,
                             active_orbitals=None, nelectron=None,
                             atom_count=None):
    """Validate the electronic-state and active-orbital contract for VRC.

    The model is deliberately independent of a particular molecule.  It can
    describe singlet, doublet, triplet, quartet, or higher-spin manifolds in
    C1 symmetry, and an optional AVAS projection makes the orbital identity
    chemical rather than an unchecked energy-window choice.
    """
    if model is None:
        return None
    if not isinstance(model, dict):
        raise ValueError('vrc_tst_state_model must be an object.')
    allowed = {'orbital_selection', 'state_average', 'selected_roots',
               'canonical_orbitals'}
    unknown = set(model) - allowed
    if unknown:
        raise ValueError('Unknown vrc_tst_state_model option(s): '
                         + ', '.join(sorted(unknown)))

    result = dict(model)
    canonical = result.get('canonical_orbitals', True)
    if not isinstance(canonical, bool):
        raise ValueError('canonical_orbitals must be true or false.')
    result['canonical_orbitals'] = canonical

    averages = result.get('state_average')
    if not isinstance(averages, list) or not averages:
        raise ValueError('state_average must contain at least one spin manifold.')
    seen = set()
    normalized = []
    for entry in averages:
        if (not isinstance(entry, dict)
                or set(entry) - {'multiplicity', 'roots', 'weights'}):
            raise ValueError('Each state_average entry accepts multiplicity, '
                             'roots, and optional weights.')
        multiplicity = entry.get('multiplicity')
        roots = entry.get('roots')
        if (isinstance(multiplicity, bool)
                or not isinstance(multiplicity, int) or multiplicity < 1
                or isinstance(roots, bool) or not isinstance(roots, int)
                or roots < 1):
            raise ValueError('State-average multiplicity and roots must be '
                             'positive integers.')
        if multiplicity in seen:
            raise ValueError('Each state-average multiplicity must be unique.')
        if nelectron is not None and (
                multiplicity - 1 > nelectron
                or (nelectron - multiplicity + 1) % 2):
            raise ValueError('A state-average multiplicity is incompatible '
                             'with the molecular electron count.')
        weights = entry.get('weights')
        if weights is None:
            weights = [1. / roots] * roots
        try:
            weights = [float(value) for value in weights]
        except (TypeError, ValueError) as exc:
            raise ValueError('State-average weights must be numeric.') from exc
        if (len(weights) != roots or any(
                not np.isfinite(value) or value <= 0. for value in weights)
                or not np.isclose(sum(weights), 1., atol=1.e-12, rtol=0.)):
            raise ValueError('State-average weights must be positive, match '
                             'the root count, and sum to one.')
        normalized.append({'multiplicity': multiplicity, 'roots': roots,
                           'weights': weights})
        seen.add(multiplicity)
    result['state_average'] = normalized

    selected = result.get('selected_roots', {})
    if not isinstance(selected, dict):
        raise ValueError('selected_roots must map multiplicities to roots.')
    selected_normalized = {}
    roots_by_mult = {entry['multiplicity']: entry['roots']
                     for entry in normalized}
    for key, value in selected.items():
        try:
            multiplicity = int(key)
        except (TypeError, ValueError) as exc:
            raise ValueError('selected_roots keys must be multiplicities.') from exc
        if str(multiplicity) != str(key) and key != multiplicity:
            raise ValueError('selected_roots keys must be integer multiplicities.')
        if (multiplicity not in roots_by_mult or isinstance(value, bool)
                or not isinstance(value, int) or value < 1
                or value > roots_by_mult[multiplicity]):
            raise ValueError('A selected root is outside its declared spin manifold.')
        selected_normalized[multiplicity] = value
    for multiplicity in roots_by_mult:
        selected_normalized.setdefault(multiplicity, 1)
    result['selected_roots'] = selected_normalized

    orbitals = result.get('orbital_selection')
    if orbitals is not None:
        if not isinstance(orbitals, dict):
            raise ValueError('orbital_selection must be an object.')
        mode = str(orbitals.get('mode', '')).casefold()
        if mode == 'molpro_record':
            allowed_record = {'mode', 'path', 'record', 'sha256'}
            unknown_record = set(orbitals) - allowed_record
            if unknown_record:
                raise ValueError('Unknown Molpro-record option(s): '
                                 + ', '.join(sorted(unknown_record)))
            path = Path(str(orbitals.get('path', ''))).expanduser().resolve()
            record = str(orbitals.get('record', '')).strip()
            digest = str(orbitals.get('sha256', '')).strip().casefold()
            if not path.is_file():
                raise ValueError(f'Molpro orbital record is unavailable: {path}.')
            if (not re.fullmatch(r'\d+\.2', record)
                    or int(record.split('.')[0]) <= 2000):
                raise ValueError('Molpro orbital record must be an orbital '
                                 'dump on file 2, such as 2140.2.')
            if not re.fullmatch(r'[0-9a-f]{64}', digest):
                raise ValueError('Molpro orbital record requires a SHA-256.')
            actual = file_sha256(path)
            if actual != digest:
                raise ValueError('Molpro orbital-record SHA-256 mismatch.')
            result['orbital_selection'] = {
                'mode': mode, 'path': str(path), 'record': record,
                'sha256': digest,
            }
            orbitals = None
        elif mode != 'avas':
            raise ValueError('orbital_selection mode must be avas or '
                             'molpro_record.')
        if orbitals is None:
            pass
        else:
            allowed_orbitals = {'mode', 'centers', 'threshold', 'minao_basis'}
            if set(orbitals) - allowed_orbitals:
                raise ValueError('Unknown orbital_selection option(s): '
                                 + ', '.join(sorted(set(orbitals)
                                                   - allowed_orbitals)))
            centers = orbitals.get('centers')
            if not isinstance(centers, list) or not centers:
                raise ValueError('AVAS orbital selection requires centers.')
            target_count = 0
            normalized_centers = []
            occupied = set()
            for center in centers:
                if (not isinstance(center, dict)
                        or set(center) != {'atom', 'orbitals'}):
                    raise ValueError('Each AVAS center requires atom and orbitals.')
                atom = center['atom']
                labels = center['orbitals']
                if (isinstance(atom, bool) or not isinstance(atom, int)
                        or atom < 1 or atom_count is not None and atom > atom_count
                        or not isinstance(labels, list) or not labels):
                    raise ValueError('Invalid AVAS atom index or orbital list.')
                normalized_labels = []
                for label in labels:
                    token = str(label).strip().casefold()
                    if (atom, token) in occupied:
                        raise ValueError('Duplicate AVAS center/orbital request.')
                    occupied.add((atom, token))
                    target_count += _avas_orbital_count(token)
                    normalized_labels.append(token)
                normalized_centers.append({'atom': atom,
                                           'orbitals': normalized_labels})
            if active_orbitals is not None and target_count != active_orbitals:
                raise ValueError(
                    f'AVAS targets {target_count} orbitals but the method requests '
                    f'{active_orbitals}.')
            threshold = orbitals.get('threshold', 0.3)
            if (isinstance(threshold, bool)
                    or not isinstance(threshold, (int, float))
                    or not np.isfinite(threshold) or not 0. < threshold < 1.):
                raise ValueError('AVAS threshold must be between zero and one.')
            minao = str(orbitals.get('minao_basis', 'minao')).strip()
            if not minao or not re.fullmatch(r'[A-Za-z0-9_+./*()-]+', minao):
                raise ValueError('AVAS minao_basis must be one safe basis identifier.')
            result['orbital_selection'] = {
                'mode': 'avas', 'centers': normalized_centers,
                'threshold': float(threshold), 'minao_basis': minao,
            }
    if active_electrons is not None and nelectron is not None:
        inactive = nelectron - active_electrons
        if inactive < 0 or inactive % 2:
            raise ValueError('The active electron count is incompatible with '
                             'the molecular electron count.')
    if active_electrons is not None and active_orbitals is not None:
        for multiplicity in seen:
            spin = multiplicity - 1
            if (active_electrons < 1 or active_orbitals < 1
                    or (active_electrons - spin) % 2
                    or spin > min(active_electrons,
                                  2 * active_orbitals - active_electrons)):
                raise ValueError('The active space cannot represent a '
                                 'state-average multiplicity.')
    return result


def stage_vrc_state_model(model, directory, stem):
    """Copy a hash-verified Molpro orbital record into one calculation."""
    model = validate_vrc_state_model(model)
    if model is None or model.get('orbital_selection', {}).get('mode') \
            != 'molpro_record':
        return model
    selection = model['orbital_selection']
    source = Path(selection['path'])
    actual = file_sha256(source)
    if actual != selection['sha256']:
        raise ValueError('Molpro orbital-record SHA-256 changed after input '
                         'validation.')
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / molpro_wfu_name(f'{stem}.initial')
    if source.resolve() != target.resolve():
        temporary = target.with_name(target.name + f'.tmp.{os.getpid()}')
        try:
            shutil.copyfile(source, temporary)
            if file_sha256(temporary) != actual:
                raise ValueError('Staged Molpro orbital record failed its '
                                 'SHA-256 check.')
            os.replace(temporary, target)
        finally:
            temporary.unlink(missing_ok=True)
    staged = deepcopy(model)
    staged['orbital_selection']['path'] = str(target.resolve())
    return staged


def _render_vrc_state_model(model, *, nelectron, active_electrons,
                            active_orbitals, spin):
    """Render AVAS and SA-CASSCF directives for one selected spin surface."""
    model = validate_vrc_state_model(
        model, active_electrons=active_electrons,
        active_orbitals=active_orbitals, nelectron=nelectron)
    if model is None:
        inactive = nelectron - active_electrons
        closed = inactive // 2
        occupied = closed + active_orbitals
        return ('', '{multi,\n'
                f'occ,{occupied}\nclosed,{closed}\n'
                f'wf,{nelectron},1,{spin}\nmaxit,50;}}\n\n', None)
    multiplicity = spin + 1
    roots = {entry['multiplicity']: entry['roots']
             for entry in model['state_average']}
    if multiplicity not in roots:
        raise ValueError(f'Multiplicity {multiplicity} is absent from the '
                         'state-average model.')
    avas = ''
    orbital_selection = model.get('orbital_selection')
    inactive = nelectron - active_electrons
    closed = inactive // 2
    occupied = closed + active_orbitals
    if (orbital_selection is not None
            and orbital_selection['mode'] == 'avas'):
        avas = ('{avas,nela=' + str(active_electrons)
                + f',thr={orbital_selection["threshold"]:g}'
                + ',minao_basis=' + orbital_selection['minao_basis'] + '\n')
        for center in orbital_selection['centers']:
            avas += ('center,' + str(center['atom']) + ','
                     + ','.join(center['orbitals']) + '\n')
        avas += '}\n\n'
        space = ''
    elif (orbital_selection is not None
          and orbital_selection['mode'] == 'molpro_record'):
        # Execution callers stage this verified record in the work directory.
        # FILE rejects path separators and can truncate long names; retain
        # the absolute source path only in the state-model provenance.
        filename = Path(orbital_selection['path']).name
        if (len(filename) > 31 or
                not re.fullmatch(r'[a-z0-9][a-z0-9_.-]*', filename)):
            raise ValueError('Stage the Molpro orbital record under a short lowercase basename before rendering.')
        avas = (f'! KinBot orbital source SHA-256: {orbital_selection["sha256"]}\n'
                f'file,2,{filename}\n')
        space = (f'start,{orbital_selection["record"]}\n'
                 f'occ,{occupied}\nclosed,{closed}\n'
                 f'diab,{orbital_selection["record"]}\n')
    else:
        space = f'occ,{occupied}\nclosed,{closed}\n'
    multi = '{multi,\n' + space
    for entry in model['state_average']:
        multi += (f'wf,{nelectron},1,{entry["multiplicity"] - 1}\n'
                  f'state,{entry["roots"]}\n'
                  'weight,' + ','.join(f'{value:.16g}'
                                        for value in entry['weights']) + '\n')
    if model['canonical_orbitals']:
        multi += 'canonical,2140.2\n'
    multi += 'maxit,50;}\n\n'
    return avas, multi, model['selected_roots'][multiplicity]


def vrc_multireference_energy_contract(method):
    """Return the exact Molpro energy selected for a VRC MR method.

    Molpro stores the F12A and F12B results of an unqualified F12 command in
    indices 1 and 2, respectively.  The literature shorthand ``CASPT2-F12``
    and ``MRCI-F12+Q`` is interpreted as F12B here.  Users can still request
    F12A explicitly.  MRCI ``+Q`` always means the relaxed-reference Davidson
    correction (``OPTION,CLUSTER=0``), including the default CABS-singles
    correction in an F12 calculation.
    """
    source = str(method).strip().casefold()
    match = _VRC_MULTIREFERENCE.fullmatch(source)
    if match is None:
        return None
    family, variant, correction, electrons_text, orbitals_text = \
        match.groups()
    if family == 'caspt2' and correction is not None:
        return None
    if family == 'mrci' and correction != '+q':
        return None
    f12 = '-f12' in source
    variant = variant.upper() if f12 and variant else ('B' if f12 else None)
    index = 2 if variant == 'B' else 1
    root = 'CASPT2' if family == 'caspt2' else 'MRCI'
    label = root + (f'-F12{variant}' if f12 else '')
    if family == 'mrci':
        label += '+Q'
    return {
        'method': label,
        'active_electrons': int(electrons_text),
        'active_orbitals': int(orbitals_text),
        'f12_variant': variant,
        'energy_variable': (f'ENERGD({index})' if family == 'mrci'
                            else f'ENERGY({index})'),
        'cabs_singles_included': (
            True if f12 and family == 'mrci' else (False if not f12 else None)),
        'davidson_correction': (
            'relaxed_reference' if family == 'mrci' else None),
        'molpro_cluster_option': 0 if family == 'mrci' else None,
    }


def _vrc_multireference_method(method, nelectron, spin=0, *, shift=None,
                               ipea=None, state_model=None):
    """Render a supported VRC active-space method and its energy variable.

    KinBot uses ``CASPT2(ne,no)`` for the inexpensive ROTD sampling surface
    and can use ``MRCI+Q(ne,no)`` for the trusted one-dimensional correction.
    Molpro obtains the Davidson-corrected MRCI energy from ``ENERGD``; ``MRCI``
    followed by ``energy`` would silently select the variational MRCI value.
    """
    contract = vrc_multireference_energy_contract(method)
    if contract is None:
        return None
    family = 'caspt2' if contract['method'].startswith('CASPT2') else 'mrci'
    active_electrons = contract['active_electrons']
    active_orbitals = contract['active_orbitals']
    if active_electrons < 1 or active_orbitals < 1:
        raise ValueError('VRC active electrons and orbitals must be positive.')
    if ((active_electrons - spin) % 2
            or spin < 0
            or spin > min(active_electrons, 2 * active_orbitals - active_electrons)):
        raise ValueError('The VRC active space cannot represent the requested spin.')
    inactive_electrons = int(nelectron) - active_electrons
    if inactive_electrons < 0 or inactive_electrons % 2:
        raise ValueError(
            f'{method} is incompatible with {nelectron} electrons: the '
            'inactive electron count must be nonnegative and even.')
    avas, multi, selected_root = _render_vrc_state_model(
        state_model, nelectron=int(nelectron),
        active_electrons=active_electrons,
        active_orbitals=active_orbitals, spin=int(spin))
    f12 = contract['f12_variant'] is not None
    selected = contract['energy_variable'].casefold()
    if family == 'caspt2':
        program = ('rs2-f12' if f12
                   else ('rs2c' if int(spin) == 0 else 'rs2'))
        options = []
        if shift is None:
            shift = 0.3
        if shift is not None:
            options.append(f'shift={float(shift):g}')
        if ipea is not None:
            options.append(f'ipea={float(ipea):g}')
        command = ','.join([program, *options])
        state = ('' if selected_root is None else
                 f'wf,{int(nelectron)},1,{int(spin)};'
                 f'state,1,{selected_root};')
        audit = (f'! KinBot energy contract: {contract["method"]} = '
                 f'{contract["energy_variable"]}\n'
                 f'KBVRCSEL={selected}\n')
        body = command if not state else f'{command};{state}'
        return avas + multi + f'{{{body}}}\n' + audit, 'KBVRCSEL'
    program = ('mrci-f12,scalf12=1,cabs_singles=1' if f12 else 'mrci')
    variant_index = 2 if contract['f12_variant'] == 'B' else 1
    audit = (f'! KinBot energy contract: {contract["method"]} = '
             f'{contract["energy_variable"]}; Davidson relaxed reference\n'
             f'KBVRCVAR=energy({variant_index})\n'
             f'KBVRCSEL={selected}\n')
    state = ('' if selected_root is None else
             f'wf,{int(nelectron)},1,{int(spin)};'
             f'state,1,{selected_root};')
    return (avas + multi
            + f'{{{program};{state}option,cluster=0;}}\n' + audit,
            'KBVRCSEL')


def validate_vrc_single_reference_root(method, state_model, multiplicity):
    """Do not label a ground-state SCF/CC result as an excited CAS root."""
    if state_model is None or vrc_multireference_energy_contract(method) is not None:
        return
    roots = state_model.get('selected_roots', {})
    selected = roots.get(multiplicity, roots.get(str(multiplicity), 1))
    if selected != 1:
        raise ValueError(f'Single-reference VRC method {method!r} cannot select '
                         f'excited CAS root {selected} for multiplicity {multiplicity}. '
                         'Use a root-resolved MR anchor, or a reviewed ground-root CC anchor.')



class Molpro:
    """
    Class to write and read molpro file and to run molpro
    """
    def __init__(self, species, par):
        self.species = species
        self.par = par

    def create_molpro_input(self, bls=0, name='', shift_vec=None, natom1=None,
                            do_vdW=False, VTS=False, sample=False,
                            vrc_level=None):
        """
        Create the input for molpro based on the template,
        which is either the one in the system, or provided
        by the user.
        : bls is whether it is a bls rection
        : name is an alternative name for the file
        : shift_vec is for bls to define the direction of shift in prod scan
        : natom1 is the number of atoms in the second fragment
        """
        if bls == 1:
            if shift_vec is None:
                tpl_file = self.par['barrierless_saddle_single_point_template']
            elif shift_vec is not None:
                tpl_file = self.par['barrierless_saddle_prod_single_point_template']
        elif VTS:
            if self.par['vrc_tst_scan_molpro_tpl'] == '':
                tpl_file = f'{kb_path}/tpl/molpro_vts.tpl'
            else:
                tpl_file = self.par['vrc_tst_scan_molpro_tpl']  # currently not really supported
        elif self.par['single_point_template'] == '':
            tpl_file = f'{kb_path}/tpl/molpro.tpl'
        else:
            tpl_file = self.par['single_point_template']

        with open(tpl_file) as f:
            tpl = f.read()

        fname = self.get_name(name, from_name=do_vdW)

        geom = ''
        nelectron = 0
        for i, at in enumerate(self.species.atom):
            x, y, z = self.species.geom[i]
            geom += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
            nelectron += atomic_numbers[str(at)]

        geometry_block = 'geometry={ \n'
        geometry_block += f"{self.species.natom};\n{self.species.name};\n{geom}\n"
        geometry_block += "}\n" 

        nelectron -= self.species.charge
        symm = self.molpro_symm()
        spin = self.species.mult - 1
        if vrc_level is not None:
            multiplicity = vrc_level.get('multiplicity', self.species.mult)
            if multiplicity is None:
                multiplicity = self.species.mult
            if (isinstance(multiplicity, bool)
                    or not isinstance(multiplicity, int)
                    or multiplicity < 1):
                raise ValueError('VRC multiplicity must be a positive integer.')
            spin = multiplicity - 1
            if spin > nelectron or (nelectron - spin) % 2:
                raise ValueError('VRC multiplicity is incompatible with the '
                                 'molecular electron count.')

        if not bls and not VTS:
            with open('molpro/' + fname + '.inp', 'w') as outf:
                outf.write(tpl.format(name=fname,
                                      natom=self.species.natom,
                                      geom=geom,
                                      nelectron=nelectron,
                                      symm=symm,
                                      spin=spin,
                                      charge=self.species.charge
                                      ))

        elif bls:
            closed = (nelectron - self.par['barrierless_saddle_nelectron']) / 2
            if closed.is_integer() is not True:
                logger.warning("The number of closed orbitals is not an integer,\n\
                             the CASPT2-like calculation will crash, but\n\
                             KinBot carries on for now. Revise your input,\n\
                             barrierless_saddle_nelectron is incorrect.")
            else:
                closed = int(closed)

            occ = closed + self.par['barrierless_saddle_norbital']

            if shift_vec is None:
                with open('molpro/' + fname + '.inp', 'w') as outf:
                    outf.write(tpl.format(name=fname,
                                          natom=self.species.natom,
                                          geom=geom,
                                          nelectron=nelectron,
                                          symm=symm,
                                          spin=spin,
                                          charge=self.species.charge,
                                          state=self.par['barrierless_saddle_nstate'],
                                          closed=closed,
                                          occ=occ
                                          ))
            else:
                shift_vec = shift_vec / np.linalg.norm(shift_vec) * 0.5  # step of 0.5 A
                geom0 = ''
                scancoo = ''
                scanstart = ''
                shift = ''
                for i, at in enumerate(self.species.atom):
                    x, y, z = self.species.geom[i]
                    if i < self.species.natom - natom1:
                        geom0 += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
                    else:
                        scancoo += '{} s{} s{} s{}\n'.format(at, 3 * i, 3 * i + 1, 3 * i + 2)
                        scanstart += 's{} = {:.8f}\ns{}= {:.8f}\ns{}= {:.8f}\n'.\
                                     format(3 * i, x, 3 * i + 1, y, 3 * i + 2, z)
                        shift += 's{0} = s{0} + {1:.8f}\n'.format(3 * i, shift_vec[0])
                        shift += 's{0} = s{0} + {1:.8f}\n'.format(3 * i + 1, shift_vec[1])
                        shift += 's{0} = s{0} + {1:.8f}\n'.format(3 * i + 2, shift_vec[2])
                with open('molpro/' + fname + '.inp', 'w') as f:
                    f.write(tpl.format(name=fname,
                                       natom=self.species.natom,
                                       geom=geom0,
                                       scanstart=scanstart,
                                       scancoo=scancoo,
                                       shift=shift,
                                       nelectron=nelectron,
                                       symm=symm,
                                       spin=spin,
                                       charge=self.species.charge,
                                       state=self.par['barrierless_saddle_nstate'],
                                       closed=closed,
                                       occ=occ
                                       ))
        elif VTS:
            scf_threshold = fortran_threshold(
                self.par.get('vrc_tst_scf_convergence', 1.e-10),
                'vrc_tst_scf_convergence')
            correlation_threshold = fortran_threshold(
                self.par.get('vrc_tst_correlation_convergence', 1.e-8),
                'vrc_tst_correlation_convergence')
            options = 'memory,3200,M\n' \
                      'GPRINT,ORBITALS,ORBEN,CIVECTOR\n' \
                      f'GTHRESH,energy={scf_threshold},' \
                      f'coeff={correlation_threshold},' \
                      f'civec={correlation_threshold}\n' \
                      'angstrom\n' \
                      'orient,noorient\n' \
                      'nosym\n' \
                      f'set,charge={self.species.charge}\nset,spin={spin}'

            if vrc_level is None:
                vrc_level = {
                    'method': (self.par['vrc_tst_sample_method'] if sample
                               else self.par['vrc_tst_high_method']),
                    'basis': (self.par['vrc_tst_sample_basis'] if sample
                              else self.par['vrc_tst_high_basis']),
                }
            basis = f"basis = {vrc_level['basis']}"
            l3_method = vrc_level['method']
            state_model = vrc_level.get(
                'state_model', self.par.get('vrc_tst_state_model'))
            state_model = stage_vrc_state_model(
                state_model, 'vrctst/molpro', fname)
            validate_vrc_single_reference_root(l3_method, state_model, spin + 1)

            method = ''
            energy_expression = 'energy(1)'
            multireference = _vrc_multireference_method(
                l3_method, nelectron, spin=spin,
                shift=vrc_level.get('shift'), ipea=vrc_level.get('ipea'),
                state_model=state_model)
            if multireference is not None:
                method, energy_expression = multireference
                # Keep the orbitals used by the correction calculation for
                # direct inspection of an external validation run.
                multi_start = method.index('{multi,')
                multi_end = method.index('}', multi_start) + 1
                method = (method[:multi_end]
                          + f'\nput,xml,{fname}.orbitals.xml\n'
                          + method[multi_end:])
            else:
                match regex_in(str(l3_method).casefold()):
                    case "uwb97xd":
                        method = " {rhf;wf," + f"{nelectron},{symm},{spin},{self.species.charge}" + "}\n\n"
                        method += " omega=0.2    !range-separation parameter\n"
                        method += " srx=0.222036 !short-range exchange\n"
                        # Same as Gaussian fine grid, important for long range
                        method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                        method += " {int; ERFLERFC,mu=$omega,srfac=$srx}\n"
                        method += " uks,HYB_GGA_XC_WB97X_D\n"
                    case "wb97xd":
                        method = " {rhf;wf," + f"{nelectron},{symm},{spin},{self.species.charge}" + "}\n\n"
                        method += " omega=0.2    !range-separation parameter\n"
                        method += " srx=0.222036 !short-range exchange\n"
                        method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                        method += " {int; ERFLERFC,mu=$omega,srfac=$srx}\n"
                        method += " ks,HYB_GGA_XC_WB97X_D\n"
                    case "ub3lyp":
                        method = " {rhf;wf," + f"{nelectron},{symm},{spin},{self.species.charge}" + "}\n\n"
                        method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                        method += " uks,HYB_GGA_XC_B3LYP\n"
                    case "b3lyp":
                        method = " {rhf;wf," + f"{nelectron},{symm},{spin},{self.species.charge}" + "}\n\n"
                        method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                        method += " ks,HYB_GGA_XC_B3LYP\n"
                    case _:
                        method += (" {rhf;wf," +
                                   f"{nelectron},{symm},{spin}" + "}\n")
                        method += "{}\n".format(l3_method)
                        energy_expression = 'energy'

            with open('vrctst/molpro/' + fname + '.inp', 'w') as f:
                contents = tpl.format(options=options,
                                      fname=fname,
                                      basis=basis,
                                      geometry_block=geometry_block,
                                      methods=method,
                                      key=self.par['vrc_tst_scan_molpro_key'].upper(),
                                      energy_expression=energy_expression)
                f.write(compact_molpro_titles(contents))
                                       
        return 0

    def get_molpro_energy(self, key, name='', do_vdW=False, VTS=False):
        """
        Verify if there is a molpro output file and if yes, read the energy
        key is the keyword for the energy we want to read
        returns 1, energy if successful
        returns 0, -1 if the energy or the file was not there
        A non-object-oriented version is used in pes.py
        """
        fname = self.get_name(name, do_vdW)
        if (self.species.natom == 1 and self.species.charge ==
                atomic_numbers[str(self.species.atom[0])]):
            return 1, 0.0
        if not VTS:
            status = os.path.exists('molpro/' + fname + '.out')
            if status:
                molpro_dir = "molpro/"
            else:
                status = os.path.exists('../molpro/' + fname + '.out')
                if status:
                    molpro_dir = "../molpro/"
                else:
                    molpro_dir = f"{fname.split('_')[0]}/molpro/"
                    status = os.path.exists(molpro_dir + fname + '.out')
        else:
            status = os.path.exists('vrctst/molpro/' + fname + '.out')
            if status:
                molpro_dir = 'vrctst/molpro/'
        if status:
            with open(f"{molpro_dir}{fname}.out") as f:
                lines = f.readlines()
            output = ''.join(lines)
            if 'Molpro calculation terminated' not in output or 'ERROR EXIT' in output:
                return 0, -1
            try:
                validate_molpro_convergence(output, allow_scf_guess_failure=VTS)
            except ValueError:
                return 0, -1
            def setting(name):
                pattern = re.compile(
                    rf'^\s*SETTING\s+{re.escape(name)}\s*=\s*'
                    r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][-+]?\d+)?)')
                for line in reversed(lines):
                    found = pattern.search(line)
                    if found:
                        return float(found.group(1).replace('D', 'E')
                                     .replace('d', 'e'))
                return None

            energy = setting(key)
            if energy is not None:
                # Generated multireference inputs bind the selected
                # ENERGY/ENERGD array element to this audit variable first.
                # Reject a result if the public key does not reproduce it.
                input_path = f"{molpro_dir}{fname}.inp"
                if os.path.isfile(input_path):
                    with open(input_path, encoding='utf-8') as source:
                        audited = 'KBVRCSEL' in source.read()
                    if audited:
                        selected = setting('KBVRCSEL')
                        if (selected is None or not np.isclose(
                                energy, selected, rtol=0., atol=1.e-10)):
                            return 0, -1
                return 1, energy
        return 0, -1

    def create_molpro_submit(self, name='', do_vdW=False, VTS=False):
        """
        write a pbs or slurm file for the molpro input file
        """
        fname = self.get_name(name, from_name=do_vdW)

        # open the template head and template
        if self.par['q_temp_l3'] == '':
            molpro_head = f'{kb_path}/tpl/{self.par["queuing"]}.tpl'
        else:
            molpro_head = self.par['q_temp_l3']
        with open(molpro_head) as f:
            tpl_head = f.read()
        from kinbot.qc import _ordinary_slurm_walltime, _slurm_exclude_directive
        walltime_directive = (_ordinary_slurm_walltime(
            self.par, self.par['queue_name'])
            if '{walltime_directive}' in tpl_head else '')
        slurm_exclude_directive = _slurm_exclude_directive(self.par)
        molpro_tpl = f'{kb_path}/tpl/{self.par["queuing"]}_molpro.tpl'
        with open(molpro_tpl) as f:
            tpl = f.read()
        # substitution
        if not VTS:
            file_string = f'molpro/{fname}.{self.par["queuing"]}'
        else:
            file_string = f'vrctst/molpro/{fname}.{self.par["queuing"]}'
        command = self.par['single_point_command']
        if VTS:
            command += ' -W ' + shlex.quote(str(Path(file_string).parent.resolve()))
        with open(file_string, 'w') as f:
            if self.par['queuing'] == 'pbs':
                f.write((tpl_head + tpl).format(
                        name=fname,
                        ppn=self.par['single_point_ppn'],
                        queue_name=self.par['queue_name'],
                        errdir='molpro',
                        command=command))
            elif self.par['queuing'] == 'slurm':
                f.write((tpl_head + tpl).format(
                        name=fname,
                        ppn=self.par['single_point_ppn'],
                        queue_name=self.par['queue_name'],
                        errdir='.',
                        command=command,
                        slurm_feature=self.par['slurm_feature'],
                        slurm_exclude_directive=slurm_exclude_directive,
                        walltime_directive=walltime_directive))
        if self.par['queuing'] == 'slurm':
            from kinbot.slurm_retry import wrap_script
            path = Path(file_string)
            path.write_text(wrap_script(path.read_text(), identity=fname, prefixes=[fname]))
        return 0

    def molpro_symm(self):
        if np.array_equal(self.species.atom, ['O']) and self.species.mult == 3:
            return 4
        if np.array_equal(self.species.atom, ['S']) and self.species.mult == 3:
            return 4
        if np.array_equal(self.species.atom, ['O', 'O']) \
                and self.species.mult == 3:
            return 4
        if np.array_equal(sorted(self.species.atom), sorted(['O', 'H'])) \
                and self.species.mult == 2:
            return 2
        return 1

    def get_name(self, name, from_name=False):
        if name != '':
            fname = name
        elif self.species.wellorts or from_name:
            fname = self.species.name
        else:
            fname = routing_name(self.species)
        return fname

@dataclass
class regex_in:
    string: str

    def __eq__(self, other: str | re.Pattern):
        if isinstance(other, str):
            other = re.compile(other)
        assert isinstance(other, re.Pattern)
        # TODO extend for search and match variants
        return other.fullmatch(self.string) is not None
