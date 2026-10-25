"""Explicit, evidence-bound SCF guesses for compact reconnaissance CASSCF.

These guesses do not select production models. Keep source canonical orbitals
on FILE,3 so a new RHF calculation cannot replace the reviewed source indices.
"""
from pathlib import Path
import json
import re

from kinbot.execution import file_sha256
from kinbot.molpro import active_orbital_swaps
from kinbot.molpro_diagnostics import _orbital_tables


def orbitals(text):
    text = re.sub(r'(?im)^\s*ELECTRON ORBITALS\s*$', 'CANONICAL ORBITALS', text)
    tables = _orbital_tables(text.splitlines(), 1, None, None)
    if not tables:
        raise ValueError('SCF seed has no printed canonical orbital table.')
    return {int(o['index'].split('.')[0]): o for o in tables[-1]['orbitals']
            if o['index'].endswith('.1')}


def validate(seed, electrons, active_electrons, active_orbitals):
    indices = seed.get('active_indices')
    if (not isinstance(indices, list) or len(indices) != active_orbitals
            or any(type(i) is not int or i < 1 for i in indices)
            or len(set(indices)) != len(indices)):
        raise ValueError('SCF seed needs distinct explicit active_indices matching the compact CAS.')
    for path_key, hash_key in (('path', 'sha256'), ('source_output', 'source_output_sha256'),
                               ('source_input', 'source_input_sha256')):
        if not isinstance(seed.get(path_key), str) or not isinstance(seed.get(hash_key), str):
            raise ValueError('SCF seed requires hashed checkpoint, input and output evidence.')
        path = Path(seed[path_key])
        if not path.is_file() or file_sha256(path) != seed[hash_key]:
            raise ValueError('SCF seed evidence changed: '+str(path))
    text = Path(seed['source_input']).read_text()
    if ('orbital,2100.2' not in text.lower() or re.search(r'(?i)\{multi\b', text)):
        raise ValueError('SCF source must be a reconnaissance RHF-only record 2100.2.')
    table = orbitals(Path(seed['source_output']).read_text(errors='replace'))
    if any(i not in table for i in indices):
        raise ValueError('Selected SCF orbital is absent from the native character evidence.')
    occupations = {i: o.get('occupation') for i, o in table.items()}
    if (any(v not in (0., 1., 2.) for v in occupations.values())
            or abs(sum(occupations.values())-electrons) > 1.e-6
            or abs(sum(occupations[i] for i in indices)-active_electrons) > 1.e-6
            or any(v == 1. and i not in indices for i, v in occupations.items())):
        raise ValueError('SCF selection must preserve the electron count and every singly occupied orbital.')
    return table


def start(seed, electrons, active_electrons):
    closed = (electrons-active_electrons)//2
    swaps = active_orbital_swaps(seed['active_indices'], closed)
    return ('start,2100.2,type=canonical,ovl\n! Source active indices: '
            + ','.join(map(str, seed['active_indices']))
            + ''.join(f'\nrotate,{a}.1,{b}.1,0' for a, b in swaps))


def separate(text, name):
    if not re.fullmatch(r'[a-z0-9][a-z0-9_.-]{0,30}', name):
        raise ValueError('Unsafe SCF seed filename.')
    if re.search(r'(?im)^\s*file,3,', text):
        raise ValueError('SCF seed needs unused FILE,3.')
    text, count = re.subn(r'(?im)^(file,2,[A-Za-z0-9_.-]+)\s*$',
                          lambda m: m[1]+',new\nfile,3,'+name, text)
    if count != 1:
        raise ValueError('SCF seed requires one destination FILE,2.')
    # Same basis and atom order: retain the coefficient ordering on the
    # moving atom-centred basis, using Molpro's normal restart. Forced OVL
    # projects spatially fixed far-separated orbitals onto the new centres;
    # large translations can lose rank before RHF or MULTI starts. This does
    # not reselect indices from the newly optimized RHF record on FILE,2.
    text, count = re.subn(r'(?im)^start,2100\.2,type=canonical(?:,ovl)?$',
                          'start,2100.3,type=canonical', text)
    if count != 1:
        raise ValueError('SCF seed requires one explicit canonical MULTI start.')
    text, count = re.subn(r'(?i)(\{rhf(?:,[^;{}]*)?;)(wf,)',
                          r'\1start,2100.3,type=canonical;\2', text, count=1)
    if count != 1:
        raise ValueError('SCF seed requires the generated RHF block.')
    return text


def projection_failure(job, *, recheck_completed=False):
    """Identify only the obsolete forced-OVL RHF initialization failure.

    A same-basis restart that itself fails is not eligible. Other native
    failures and unbound/stale attempts are never reset. Completed diagnostic
    probes require an explicit, job-ID-scoped recheck by the review owner.
    The caller owns the manifest lock and checks scheduler liveness.
    """
    space = job.get('probe_active_space') or {}
    seed = job.get('initial_record') or {}
    if (space.get('orbital_source') != 'scf'
            or (job.get('status') == 'complete' and not recheck_completed)
            or not job.get('job_id') or job.get('submission_status') == 'submitting'
            or seed.get('basis', '').lower() != job.get('basis', '').lower()
            or any(h.get('input_transform') == 'native_xml_scf_seed_file_v7'
                   for h in job.get('attempt_history', []))):
        return None
    directory = Path(job['directory'])
    record = directory/'execution.json'
    if not record.exists():
        return None
    execution = json.loads(record.read_text())
    native = directory/(job['name']+'.native.inp')
    output = directory/(job['name']+'.out')
    expected_status = 'complete' if recheck_completed else 'failed'
    if (execution.get('status') != expected_status
            or str(execution.get('slurm_job_id')) != str(job['job_id'])
            or execution.get('input_transform') != 'native_xml_scf_seed_file_v7'
            or execution.get('input_sha256') != job['input_sha256']
            or execution.get('native_input') != str(native)
            or not native.is_file() or not output.is_file()
            or file_sha256(native) != execution.get('native_input_sha256')
            or any(execution.get('scf_seed', {}).get(k) != v for k, v in seed.items())):
        return None
    text = output.read_text(errors='replace')
    programs = re.findall(r'(?im)^\s*PROGRAM\s*\*\s*([^\n]+)', text)
    if (not any(re.match(r'(?:RHF(?:-SCF)?\b|Restricted Hartree-Fock\b)', p, re.I) for p in programs)
            or 'start,2100.3,type=canonical,ovl' not in native.read_text().lower()):
        return None
    if recheck_completed:
        if (job.get('status') != 'complete' or execution.get('returncode') != 0
                or file_sha256(output) != job.get('output_sha256')
                or not any(p.upper().startswith('MULTI') for p in programs)
                or not re.search(r'(?im)^\s*start,2100\.3,type=canonical,ovl\s*$', native.read_text())):
            return None
    elif (any(p.upper().startswith('MULTI') for p in programs)
          or not re.search(r'(?im)^\s*\?\s*ERROR IN VIRTORB:\s*INCORRECT NUMBER OF ORBITALS:', text)):
        return None
    validate(seed, job['electron_count'], space['active_electrons'], space['active_orbitals'])
    return dict(kind=('reviewed_seed_recheck' if recheck_completed else 'same_basis_seed_restart'),
                previous_transform=execution['input_transform'],
                input_sha256=job['input_sha256'], native_input_sha256=file_sha256(native),
                output_sha256=file_sha256(output), source_sha256=seed['sha256'],
                reason=('Explicit review recheck of forced-OVL SCF initialization; retain original result.'
                        if recheck_completed else
                        'VIRTORB failed during forced-OVL RHF initialization; use normal same-basis restart.'))
