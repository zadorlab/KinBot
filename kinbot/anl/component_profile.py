"""Versioned executable ANL component choices, shared by planning and assembly.

Version 1 graphs keep their original scientific meaning. Version 2 uses F12b
for eligible CCSD(T) components and the tier-specific DBOC prescription.
"""

from copy import deepcopy
import re


COMPONENT_VERSION = 2


def cv_auxiliary(overrides=None):
    """Explicit weighted-CV fitting sets; no nonexistent pwCVnZ/OptRI alias."""
    result = {n: {'df_basis': f'aug-cc-pwCV{n}Z/mp2fit',
                  'df_basis_exch': 'aug-cc-pV5Z/jkfit',
                  'ri_basis': 'aug-cc-pwCV5Z/mp2fit'} for n in ('T', 'Q', '5')}
    overrides = {} if overrides is None else overrides
    if not isinstance(overrides, dict) or set(overrides)-set(result):
        raise ValueError('ANL CV auxiliary overrides accept T, Q, 5 cardinal keys.')
    for cardinal, values in overrides.items():
        if (not isinstance(values, dict) or set(values)-set(result[cardinal])
                or any(not isinstance(v, str) or not re.fullmatch(r'[A-Za-z0-9_+./*-]+', v)
                       for v in values.values())):
            raise ValueError('ANL CV auxiliaries require safe df_basis/df_basis_exch/ri_basis names.')
        result[cardinal].update(values)
    return result


def component_profile(name, augmented=False, *, cv_f12_auxiliary=None):
    if name not in ('ANL0', 'ANL0-F12', 'ANL1', 'ANL1-F12'):
        raise ValueError('Unknown ANL component profile.')
    if not isinstance(augmented, bool):
        raise ValueError('ANL augmentation must be boolean.')
    high = name.startswith('ANL1')
    f12 = name.endswith('-F12')
    prefix = 'aug-' if augmented else ''
    method = 'CCSD(T)-F12b' if f12 else 'CCSD(T)'
    suffix = '-F12' if f12 else ''
    geometry = tuple(f'{prefix}cc-pV{x}Z{suffix}' for x in ('T', 'Q'))
    # Weighted core-valence ORBITAL bases are shared by the conventional and
    # F12b differences. The method suffix does not imply a "pwCVnZ-F12"
    # orbital family. AE and FC must use identical orbital/auxiliary sets.
    cv_cardinals = ('Q', '5') if high else ('T', 'Q')
    cv = tuple(f'{prefix}cc-pwCV{x}Z' for x in cv_cardinals)
    rel = (tuple(f'{prefix}cc-pwCV{x}Z-DK' for x in ('Q', '5')) if high
           else (f'{prefix}cc-pCVTZ-DK',))
    return {
        'version': COMPONENT_VERSION, 'name': name, 'augmented': augmented,
        'method': method, 'geometry_bases': geometry,
        'geometry_basis': geometry[1 if high else 0],
        'reference_bases': tuple(f'{prefix}cc-pV{x}Z' for x in
                                 (('5', '6') if high else ('Q', '5'))),
        'f12_reference_bases': tuple(f'{prefix}cc-pV{x}Z-F12' for x in
                                     (('Q', '5') if high else ('T', 'Q'))),
        'core_valence_bases': cv,
        'core_valence_cardinals': (4, 5) if high else (3, 4),
        'core_valence_auxiliary': ({basis: cv_auxiliary(cv_f12_auxiliary)[n]
                                   for n, basis in zip(cv_cardinals, cv)} if f12 else {}),
        # DKH/X2C + F12 is explicitly unsupported by Molpro. BOTH sides of
        # this difference must remain conventional and use the same basis.
        'relativistic_method': 'CCSD(T)',
        'relativistic_bases': rel,
        'relativistic_basis': f'CBS({rel[0]},{rel[1]})' if high else rel[0],
        'dboc_method': 'CCSD' if f12 else 'HF',
        'dboc_basis': f'{prefix}cc-pV{"T" if high else "D"}Z',
        'higher_bases': (f'{prefix}cc-pVDZ', f'{prefix}cc-pVTZ'),
    }


def _change_basis(task, basis):
    parser = task.get('result_parser', {})
    if 'basis' in parser:
        old = parser['basis']
        task['input_template'] = task['input_template'].replace(old, basis)
        parser['basis'] = basis
    elif task.get('backend') == 'mrcc':
        raise ValueError('MRCC task lacks a basis-specific parser.')


def apply_component_profile(tasks, name, augmented=False, *, multiplicity=1,
                            cv_f12_auxiliary=None):
    """Change new task definitions before hashing/staging; never live files."""
    p = component_profile(name, augmented, cv_f12_auxiliary=cv_f12_auxiliary)
    f12 = name.endswith('-F12')
    high = name.startswith('ANL1')
    tasks = deepcopy(tasks)
    renamed = ({f'cv_{core}_{old}': f'cv_{core}_{new}'
                for core in ('ae', 'fc') for old, new in (('tz', 'qz'), ('qz', '5z'))}
               if high else {})
    result = []
    for task in tasks:
        ident = task['id']
        if ident in ('l3_geometry', 'l3_geometry_tz'):
            basis = p['geometry_bases'][0] if ident.endswith('_tz') else p['geometry_basis']
            profile = task['profile']
            profile.update(method=p['method'], basis=basis,
                           label=p['method'] + '/' + basis)
            kw = profile.setdefault('calculator_kwargs', {})
            if f12 and kw.get('effective_method') == 'CCSD':
                kw['effective_method'] = 'CCSD-F12b'
        elif ident in ('harmonic', 'harmonic_tz', 'harmonic_qz'):
            _change_basis(task, p['geometry_bases'][1 if ident.endswith('_qz') else 0])
            if f12:
                _make_f12(task)
        elif ident.startswith('cv_'):
            _change_basis(task, p['core_valence_bases'][ident.endswith('_qz')])
            if f12:
                _make_f12(task)
                if task['result_parser'].get('effective_method') != 'HF':
                    auxiliary = p['core_valence_auxiliary'][task['result_parser']['basis']]
                    task['result_parser']['f12_auxiliary'] = dict(auxiliary)
                    options = ','.join(f'{key}={value}' for key, value in auxiliary.items())
                    for command in ('uccsd(t)-f12b,', 'uccsd-f12b,'):
                        task['input_template'] = task['input_template'].replace(
                            command, command + options + ',')
        elif ident in ('rel_dkh', 'rel_nonrel'):
            task['result_parser']['explicit_hamiltonian'] = True
            # Disable auto-DKH explicitly on the NR arm even on Molpro
            # versions that infer a Hamiltonian from the -DK basis suffix.
            if ident == 'rel_nonrel':
                task['input_template'] = task['input_template'].replace(
                    '\nrhf\n', '\nset,dkroll=0\nrhf\n')
            if high:
                for cardinal, basis in zip(('qz', '5z'), p['relativistic_bases']):
                    sibling = deepcopy(task)
                    _change_basis(sibling, basis)
                    _rename_task(sibling, ident+'_'+cardinal)
                    sibling['failure_skips'] = [f'rel_{arm}_{z}'
                        for arm in ('dkh', 'nonrel') for z in ('qz', '5z')
                        if f'rel_{arm}_{z}' != sibling['id']]
                    result.append(sibling)
                continue
            _change_basis(task, p['relativistic_bases'][0])
        elif ident == 'l3_sp_l2':
            _change_basis(task, p['geometry_bases'][0])
            if f12:
                _make_f12(task)
        elif ident == 'cfour_dboc':
            _change_basis(task, p['dboc_basis'])
            task['result_parser']['level'] = p['dboc_method']
            # CFOUR supports ROHF DBOC only at HF. UHF is an explicit,
            # component-local exception for open-shell CCSD DBOC alone.
            reference = ('RHF' if multiplicity == 1 else 'UHF' if f12 else 'ROHF')
            task['result_parser']['reference'] = reference
            task['input_template'] = task['input_template'].replace(
                'DBOC=ON\n', f'DBOC=ON\nREFERENCE={reference}\nFROZEN_CORE=OFF\n')
            if f12:
                task['input_template'] = task['input_template'].replace(
                    'CALC=SCF', 'CALC=CCSD\nCC_CONV=10\nLINEQ_CONV=10')
        elif ident in ('anl0_qz', 'anl0_5z', 'anl1_5z', 'anl1_6z'):
            _change_basis(task, p['reference_bases'][ident in ('anl0_5z', 'anl1_6z')])
        elif ident.startswith('f12_'):
            cardinal = ident.split('_')[1].upper().replace('Z', '')
            _change_basis(task, ('aug-' if augmented else '') + f'cc-pV{cardinal}Z-F12')
        elif ident in ('ccsdt_dz', 'ccsdt_tz', 'ccsdtq_dz', 'ccsdtq_tz', 'ccsdtqp_dz'):
            _change_basis(task, p['higher_bases'][ident.endswith('_tz')])
        # Gaussian VPT2 remains the specified L2 vibrational correction.
        # Gaussian/CFOUR do not provide a CCSD(T)-F12b VPT2 implementation.
        for key in ('depends_on', 'failure_skips', 'fallback_on'):
            if key in task:
                task[key] = [renamed.get(value, value) for value in task[key]]
        if ident in renamed:
            _rename_task(task, renamed[ident])
        result.append(task)
    return result


def _rename_task(task, name):
    """Rename native artifacts with their cardinal, never just the UI label."""
    old = task['id']
    task['id'] = name
    for key in ('input_name', 'input_template'):
        task[key] = task[key].replace(old, name)
    task['required_outputs'] = [v.replace(old, name) for v in task['required_outputs']]
    for key in ('result_parser', 'success_marker'):
        if key in task:
            task[key]['file'] = task[key]['file'].replace(old, name)


def _make_f12(task):
    parser = task['result_parser']
    parser['method'] = 'CCSD(T)-F12b'
    text = task['input_template']
    if parser.get('effective_method') == 'CCSD':
        parser['effective_method'] = 'CCSD-F12b'
        text = text.replace('uccsd,', 'uccsd-f12b,').replace('uccsd;', 'uccsd-f12b;')
    elif parser.get('effective_method') != 'HF':
        text = text.replace('uccsd(t)', 'uccsd(t)-f12b,scale_trip=1')
    task['input_template'] = text
