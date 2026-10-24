"""Verified L3 thermal models without overwriting accepted L2 calculations.

Molpro prints Cartesian (mass-unweighted) normal modes. Their complete
mass-weighted spectral decomposition supplies the vibrational Hessian, to
native print precision, for the same internal-rotor projection used by KinBot.
Gaussian modes are used only for correspondence of L2 VPT2 shifts.
"""
from copy import copy, deepcopy
import re

import numpy as np
from ase.data import atomic_numbers
from ase.io import read
from ase.units import Bohr
from scipy.optimize import linear_sum_assignment

from kinbot import constants, frequencies
from kinbot.calculation import array_fingerprint


def _number(value):
    value = value.replace('D', 'E').replace('d', 'e')
    return -abs(float(value[:-1])) if value.lower().endswith('i') else float(value)


def _rotation(source, target):
    x, y = source-source.mean(0), target-target.mean(0)
    u, _, vt = np.linalg.svd(x.T @ y)
    sign = np.eye(3); sign[-1, -1] = np.linalg.det(u @ vt)
    return u @ sign @ vt


def _mass_modes(vectors, symbols, *, tolerance):
    if not len(vectors):
        raise ValueError('Native normal-mode table is empty.')
    q = np.asarray(vectors, float) * np.sqrt(np.repeat(
        [constants.exact_mass[str(s)] for s in symbols], 3))
    if not np.isfinite(q).all() or np.any(np.linalg.norm(q, axis=1) == 0):
        raise ValueError('Native normal modes are nonfinite or zero.')
    q /= np.linalg.norm(q, axis=1)[:, None]
    if np.max(abs(q @ q.T-np.eye(len(q)))) > tolerance:
        raise ValueError('Native normal modes are not mass-orthogonal within print precision.')
    return q


def molpro_modes(text, symbols, parsed):
    """Read all stable/imaginary vectors, excluding rigid low/zero modes."""
    section = text.split('FREQUENCIES * CALCULATION OF NORMAL MODES', 1)
    if len(section) != 2:
        raise ValueError('Molpro output lacks the native normal-mode section.')
    section = section[1]
    coords = re.search(r'Atomic Coordinates\s+Nr\s+Atom\s+Charge\s+X\s+Y\s+Z\s*\n((?:\s*\d+\s+\S+\s+[^\n]+\n)+)', section)
    if coords is None:
        raise ValueError('Molpro normal-mode geometry is absent.')
    rows = [row.split() for row in coords[1].splitlines()]
    if [r[1] for r in rows] != list(symbols):
        raise ValueError('Molpro normal-mode atom order differs.')
    positions = np.array([list(map(_number, r[3:6])) for r in rows])*Bohr
    values, vectors = [], []
    sections = re.split(r'^\s*Normal Modes([^\n]*)\n', section, flags=re.M)
    for title, body in zip(sections[1::2], sections[2::2]):
        if 'low' in title.lower() or 'zero' in title.lower():
            continue
        for block in re.split(r'^\s*Wavenumbers \[cm-1\]\s*', body, flags=re.M)[1:]:
            lines = block.splitlines()
            freq = [_number(v) for v in lines[0].split()]
            if 'imaginary' in title.lower():
                freq = [-abs(v) for v in freq]
            columns = []
            for line in lines[1:]:
                row = re.match(r'^\s*([A-Za-z]+)([XYZ])(\d+)\s+(.+)$', line)
                if row:
                    index = len(columns)
                    if (row[1] != symbols[index//3] or row[2] != 'XYZ'[index%3]
                            or int(row[3]) != index//3+1):
                        raise ValueError('Molpro normal-mode Cartesian labels differ.')
                    column = [_number(v) for v in row[4].split()]
                    if len(column) != len(freq):
                        raise ValueError('Molpro normal-mode block is incomplete.')
                    columns.append(column)
                    if len(columns) == 3*len(symbols):
                        break
            if len(columns) != 3*len(symbols):
                raise ValueError('Molpro normal-mode Cartesian rows are incomplete.')
            values.extend(freq); vectors.extend(np.array(columns).T.tolist())
    expected = sorted(parsed['wavenumbers_cm_inverse']+parsed['imaginary_wavenumbers_cm_inverse'])
    if len(values) != len(expected) or not np.allclose(sorted(values), expected, atol=.03, rtol=0):
        raise ValueError('Molpro native vectors and accepted frequency list disagree.')
    if len(values) not in (3*len(symbols)-6, 3*len(symbols)-5):
        raise ValueError('Molpro normal modes do not span the vibrational space.')
    order = np.argsort(values)
    q = _mass_modes(vectors, symbols, tolerance=.003)[order]
    # Remove only printed-vector roundoff; do not repair a wrong subspace.
    u, _, vt = np.linalg.svd(q, full_matrices=False)
    return positions, np.array(values)[order], u @ vt


def gaussian_modes(text, symbols, target):
    """Read the final harmonic vector table and verify its coordinate frame."""
    headers = list(re.finditer(r'Harmonic frequencies \(cm\*\*-1\)', text))
    if not headers:
        raise ValueError('L3/VPT2 matching needs Gaussian harmonic normal vectors.')
    start = headers[-1].start()
    orientations = list(re.finditer(r'(?:Input|Standard) orientation:\s*\n', text[:start]))
    if not orientations:
        raise ValueError('Gaussian harmonic mode orientation is absent.')
    body = text[orientations[-1].end():start]
    rows = re.findall(r'^\s*(\d+)\s+(\d+)\s+-?\d+\s+([-+.\dDEde]+)\s+([-+.\dDEde]+)\s+([-+.\dDEde]+)\s*$', body, re.M)
    # The first orientation table; later properties can also contain numbers.
    rows = rows[:len(symbols)]
    if [int(r[1]) for r in rows] != [atomic_numbers[s] for s in symbols]:
        raise ValueError('Gaussian harmonic mode atom order differs.')
    coords = np.array([[float(v.replace('D','E')) for v in r[2:]] for r in rows])
    rotation = _rotation(coords, target)
    if np.max(np.linalg.norm((coords-coords.mean(0)) @ rotation-(target-target.mean(0)), axis=1)) > 2.e-4:
        raise ValueError('Gaussian harmonic vectors do not belong to the accepted L2 geometry.')
    values, vectors = [], []
    for block in re.split(r'^\s*Frequencies --\s*', text[start:], flags=re.M)[1:]:
        lines = block.splitlines(); freq = list(map(_number, lines[0].split()))
        table = re.split(r'\bAtom\s+AN\s+', '\n'.join(lines[1:]), maxsplit=1)
        if len(table) != 2:
            raise ValueError('Gaussian Cartesian normal-mode table is absent.')
        rows = []
        for line in table[1].splitlines()[1:]:
            fields = line.split()
            if len(fields) != 2+3*len(freq) or not fields[0].isdigit():
                if rows: break
                continue
            if int(fields[0]) != len(rows)+1 or int(fields[1]) != atomic_numbers[symbols[len(rows)]]:
                raise ValueError('Gaussian normal-mode Cartesian labels differ.')
            rows.append(list(map(_number, fields[2:])))
            if len(rows) == len(symbols): break
        if len(rows) != len(symbols):
            raise ValueError('Gaussian Cartesian normal-mode rows are incomplete.')
        v = np.array(rows).reshape(len(symbols),len(freq),3).transpose(1,0,2) @ rotation
        values.extend(freq); vectors.extend(v.reshape(len(freq),-1))
    q = _mass_modes(vectors, symbols, tolerance=.08)
    return np.array(values), q


def handoff_options(run):
    """Pick the highest completed L3 geometry and Hessian independently."""
    from kinbot.anl.dispatch import _load
    run, spec, state = _load(run)
    complete = {k for k,v in state['tasks'].items() if v['status']=='complete'}
    geometry = next((k for k in ('l3_geometry','l3_geometry_tz') if k in complete), None)
    harmonic = next((k for k in ('harmonic_qz','harmonic','harmonic_tz') if k in complete), None)
    if geometry is None and harmonic is None:
        return None
    return dict(run=str(run), geometry_task=geometry, harmonic_task=harmonic)


def attach_model(point, options):
    from kinbot.anl.dispatch import _load
    from kinbot.anl.ts_reuse import aligned_permutation
    from kinbot.anl.validation import molecule_from_completed_run
    from kinbot.anl.workflow import _verified_task_result
    run, spec, state = _load(options['run'])
    initial = spec['molecule']
    if state['tasks'].get('l2_geometry', {}).get('status') == 'complete':
        initial = molecule_from_completed_run(run, 'l2_geometry')
    symbols = initial['symbols']
    if (initial.get('charge',0), initial.get('multiplicity',1)) != (point.charge,point.mult):
        raise ValueError('ANL thermal model electronic state differs.')
    x, y = np.array(initial['positions']), np.array(point.geom)
    mapping = aligned_permutation(list(point.atom), y-y.mean(0), symbols, x-x.mean(0), 2.e-4)
    if mapping is None:
        raise ValueError('ANL thermal model cannot map its accepted L2 atoms to this species.')
    permutation = mapping[1]
    model = dict(schema=1, run=str(run), atom_map=permutation, l2_geometry_sha256=array_fingerprint(point.geom),
                 geometry=None, harmonic=None, vpt2=None)
    geometry_task = options.get('geometry_task')
    if geometry_task:
        molecule = molecule_from_completed_run(run, geometry_task)
        coords = np.array(molecule['positions'])[permutation]
        rot = _rotation(coords,y)
        model['geometry'] = dict(positions=((coords-coords.mean(0)) @ rot+y.mean(0)).tolist(),
                                 source=molecule['source'])
        from kinbot.symmetry import conformer_symmetry_numbers
        model['geometry']['sigma_ext'] = conformer_symmetry_numbers(
            point, model['geometry']['positions'])['sigma_ext']
    task_id = options.get('harmonic_task')
    if task_id:
        _, _, task, execution, native, parsed = _verified_task_result(run, task_id)
        if parsed['kind'] != 'molpro_harmonic' or parsed.get('review_required'):
            raise ValueError('L3 thermal Hessian is not an accepted Molpro harmonic result.')
        coords, freq, q = molpro_modes(native.read_text(errors='replace'), symbols, parsed)
        staged = read(run/'tasks'/task_id/'geometry.xyz')
        rot = _rotation(coords,staged.positions)
        if np.max(np.linalg.norm((coords-coords.mean(0)) @ rot-(staged.positions-staged.positions.mean(0)),axis=1)) > 2.e-5:
            raise ValueError('Molpro normal-mode geometry differs from its verified task.')
        coords = coords[permutation]; rot = _rotation(coords,y)
        q = (q.reshape(len(q),len(symbols),3)[:,permutation] @ rot).reshape(len(q),-1)
        model['harmonic'] = dict(positions=((coords-coords.mean(0)) @ rot+y.mean(0)).tolist(),
            frequencies=freq.tolist(), mass_weighted_modes=q.tolist(), method=parsed['method'], basis=parsed['basis'],
            task_id=task_id, native_output=str(native), source_sha256=execution['artifacts'][task['result_parser']['file']],
            representation='complete native normal-mode spectral Hessian; printed precision')
    # Attachment is reconstructed from the verified artifacts, not accumulated.
    point.__dict__.pop('anl_thermochemistry_frequencies', None)
    point.__dict__.pop('anl_thermochemistry_frequency_source', None)
    point.anl_thermal_model = model
    update_model(point)
    return provenance(point)


def attach_vpt2(point, parsed, native, staged_geometry, *, mode_map=None, match_tolerance_cm_inverse=100.):
    """Retain the L2 shift source; match by mass-weighted vector overlap."""
    model = point.anl_thermal_model
    permutation = model['atom_map']
    symbols = staged_geometry.get_chemical_symbols()
    f, q = gaussian_modes(native.read_text(errors='replace'), symbols, staged_geometry.positions)
    harmonic = np.asarray(parsed['harmonic_fundamentals_cm_inverse'])
    rows, cols = linear_sum_assignment(abs(harmonic[:,None]-f[None,:]))
    if len(rows) != len(f) or len(harmonic) != len(f) or max(abs(harmonic[rows]-f[cols])) > .1:
        raise ValueError('Gaussian VPT2 table cannot be matched to its harmonic vectors.')
    q = q[cols]
    coords = staged_geometry.positions[permutation]
    if [symbols[i] for i in permutation] != list(point.atom):
        raise ValueError('VPT2 atom order differs from the accepted thermal model.')
    rot = _rotation(coords,np.asarray(point.geom))
    if np.max(np.linalg.norm((coords-coords.mean(0)) @ rot
                            -(point.geom-point.geom.mean(0)), axis=1)) > 2.e-4:
        raise ValueError('VPT2 geometry differs from the accepted L2 thermal model.')
    q = (q.reshape(len(q),len(symbols),3)[:,permutation] @ rot).reshape(len(q),-1)
    model['vpt2'] = dict(shifts=(np.asarray(parsed['anharmonic_fundamentals_cm_inverse'])-harmonic).tolist(),
        harmonic=harmonic.tolist(), mass_weighted_modes=q.tolist(), mode_map=mode_map,
        match_tolerance_cm_inverse=float(match_tolerance_cm_inverse),
        method=parsed['method'], basis=parsed['basis'], dispersion=parsed['dispersion'])
    update_model(point)
    return point.anl_thermochemistry_frequencies


class _SelectedRotors:
    def __init__(self, indices): self.indices = indices
    def is_valid_rotor(self, index): return index in self.indices
    def invalid_rotor_reason(self, index): return 'not in the accepted L2 rotor model'


def update_model(point):
    """Reproject if the owning MESS writer repairs/demotes an L2 rotor."""
    model = point.anl_thermal_model
    if array_fingerprint(point.geom) != model['l2_geometry_sha256']:
        raise ValueError('L2 geometry changed after ANL thermal handoff; reconstruct the handoff.')
    h = model['harmonic']
    if h is None:
        return
    projection = getattr(point,'rotor_projection',None) or {}
    selected = [r['rotor_index'] for r in projection.get('rotors',[]) if r['projected']]
    if not projection and len(point.reduced_freqs) != len(point.freq):
        raise ValueError('L3 thermal projection needs the accepted L2 rotor identities.')
    signature = (tuple(selected), tuple(tuple(r) for r in getattr(point,'dihed',[])), repr(model['vpt2']))
    if model.get('_cache') == signature:
        return
    q = np.asarray(h['mass_weighted_modes']); values = np.asarray(h['frequencies'])
    # Invert the same unit conversion used by KinBot's projected Hessians.
    scale = frequencies.convert_to_wavenumbers(1.)
    eig = np.sign(values)*(values/scale)**2
    hmw = (q.T*eig) @ q
    view = copy(point); view.geom = np.asarray(h['positions']); view.hir = _SelectedRotors(selected)
    view.dihed = deepcopy(getattr(point,'dihed',[]))
    _, reduced, _, modes = frequencies.get_frequencies(view,hmw,view.geom,massweighted=True,return_modes=True)
    saddle = bool(getattr(point,'wellorts',False))
    if (sum(v<0 for v in reduced) != int(saddle) or any(v==0 for v in reduced)
            or len(reduced) != len(values)-len(selected)):
        raise ValueError('L3 rotor-projected spectrum has the wrong dimension or stationary-point character.')
    source = dict(model='l3_harmonic', harmonic={k:v for k,v in h.items()
                  if k not in ('positions','mass_weighted_modes')},
                  projected_l2_rotors=selected, l2_geometry_sha256=model['l2_geometry_sha256'])
    if model['vpt2'] is not None:
        v = model['vpt2']; shifts = np.array(v['shifts']); target = np.array(reduced)
        overlaps = (modes @ np.array(v['mass_weighted_modes']).T)**2
        if v['mode_map'] is None:
            rows, cols = linear_sum_assignment(-overlaps)
            if len(cols) != len(target): raise ValueError('VPT2 has fewer modes than the retained L3 spectrum.')
            for i,j in zip(rows,cols):
                rivals = [k for k in range(len(shifts)) if k != j and overlaps[i,k] >= overlaps[i,j]-.1]
                if overlaps[i,j] < .5 or any(abs(shifts[k]-shifts[j])>1. for k in rivals):
                    raise ValueError('Ambiguous L3/L2 VPT2 mode overlap; provide a reviewed explicit mode map.')
        else:
            mapping = v['mode_map']
            if (len(mapping)!=len(target) or any(type(i) is not int for i in mapping)
                    or len(set(mapping))!=len(mapping) or min(mapping)<1 or max(mapping)>len(shifts)):
                raise ValueError('VPT2 mode map must select distinct valid modes for the L3 spectrum.')
            cols = np.array(mapping)-1
        mismatch = abs(target-np.array(v['harmonic'])[cols])
        if (not np.isfinite(v['match_tolerance_cm_inverse'])
                or v['match_tolerance_cm_inverse'] <= 0
                or np.any(mismatch > v['match_tolerance_cm_inverse'])):
            raise ValueError('L3/L2 harmonic mismatch exceeds the reviewed VPT2 matching tolerance.')
        corrected = target+shifts[cols]
        if np.any(corrected <= 0): raise ValueError('L3 plus L2 VPT2 shifts gives nonpositive stable frequencies.')
        reduced = corrected.tolist()
        source.update(model='l3_harmonic_plus_l2_vpt2_mode_shifts', mode_map=(cols+1).tolist(),
            squared_overlaps=overlaps[np.arange(len(cols)),cols].tolist(),
            l2_mode_shifts_cm_inverse=shifts[cols].tolist(),
            maximum_harmonic_mismatch_cm_inverse=float(max(mismatch)),
            vpt2_method=v['method'], vpt2_basis=v['basis'],
            **v.get('evidence', {}))
    point.anl_thermochemistry_frequencies = tuple(reduced)
    point.anl_thermochemistry_frequency_source = source
    model['_cache'] = signature


def provenance(point):
    model = point.anl_thermal_model
    return dict(run=model['run'], atom_map=model['atom_map'],
                geometry=deepcopy(model['geometry']),
                frequencies=deepcopy(getattr(point,'anl_thermochemistry_frequency_source',None)),
                hindered_rotor_geometry='accepted L2 scan geometry')


def member_spectrum(point):
    """The selected conformer's RRHO spectrum, without HR mode removal."""
    if not point.anl_thermal_model.get('harmonic'):
        return point.freq
    view = copy(point)
    view.anl_thermal_model = deepcopy(point.anl_thermal_model)
    view.rotor_projection = {'rotors': []}
    update_model(view)
    return getattr(view, 'anl_thermochemistry_frequencies', point.freq)


def enrich_reference_model(model, run, *, vpt2_review=None):
    """Use the same verified thermal model for an external MESSPF reference.

    Legacy exports identify rotors by their atom axes. Their L2 potential and
    inertia geometry are retained verbatim; only the main RRHO is replaced.
    """
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.anl.dispatch import _load
    from kinbot.anl.workflow import attach_task_vpt2_frequencies, _verified_task_result
    options = handoff_options(run)
    _, _, state = _load(run)
    if options is None and state['tasks'].get('gaussian_vpt2', {}).get('status') != 'complete':
        return deepcopy(model)
    out = deepcopy(model)
    point = StationaryPoint(name=model.get('name', 'reference'),
        charge=model['charge'], mult=model['multiplicity'], smiles=model['smiles'],
        natom=len(model['symbols']), atom=model['symbols'], geom=np.array(model['positions']))
    point.characterize()
    blocks = model.get('hindered_rotor_blocks', '')
    axes = [tuple(int(n)-1 for n in pair) for pair in
            re.findall(r'^\s*Axis\s+(\d+)\s+(\d+)\s*$', blocks, re.M)]
    selected = []
    for axis in axes:
        found = [i for i,rot in enumerate(point.dihed) if set(rot[1:3]) == set(axis)]
        if len(found) != 1 or found[0] in selected:
            raise ValueError('External L2 rotor axes cannot be mapped uniquely to the L3 model.')
        selected.extend(found)
    point.rotor_projection = {'rotors': [dict(rotor_index=i, projected=True) for i in selected]}
    point.reduced_freqs = list(model['frequencies_cm_inverse'])
    point.freq = point.reduced_freqs + [0.]*len(selected)
    if options:
        attach_model(point, options)
    if state['tasks'].get('gaussian_vpt2', {}).get('status') == 'complete':
        parsed = _verified_task_result(run, 'gaussian_vpt2')[-1]
        if vpt2_review or not parsed.get('review_required'):
            attach_task_vpt2_frequencies(point, run, 'gaussian_vpt2', review_file=vpt2_review)
    geom = getattr(point, 'anl_thermal_model', {}).get('geometry')
    if geom:
        from kinbot.symmetry import conformer_symmetry_numbers
        l2_sigma = conformer_symmetry_numbers(point, point.geom)['sigma_ext']
        out['positions'] = geom['positions']
        out['symmetry_factor'] *= geom['sigma_ext']/l2_sigma
    out['frequencies_cm_inverse'] = list(getattr(point, 'anl_thermochemistry_frequencies', point.reduced_freqs))
    if blocks:
        # An explicit Geometry belongs to each rotor, not the enclosing RRHO.
        geometry = ('\n        Geometry[angstrom] ' + str(point.natom) + '\n' +
            '\n'.join('          '+s+' '+' '.join(f'{x:.10f}' for x in p)
                      for s,p in zip(point.atom, point.geom)))
        parts = re.split(r'(^\s*Rotor\s+(?:Hindered|Free)\b[^\n]*)', blocks, flags=re.M)
        for i in range(1, len(parts), 2):
            if not re.search(r'\bGeometry\[', parts[i+1]):
                parts[i] += geometry
        out['hindered_rotor_blocks'] = ''.join(parts)
    if options:
        out['thermal_source'] = provenance(point)
    elif hasattr(point, 'anl_thermochemistry_frequency_source'):
        out['frequency_source'] = point.anl_thermochemistry_frequency_source
    return out
