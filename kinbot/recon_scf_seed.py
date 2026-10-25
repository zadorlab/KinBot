"""Explicit, evidence-bound SCF guesses for compact reconnaissance CASSCF.

These guesses do not select production models. Keep source canonical orbitals
on FILE,3 so a new RHF calculation cannot replace the reviewed source indices.
"""
from pathlib import Path
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
    text, count = re.subn(r'(?im)^start,2100\.2,type=canonical,ovl$',
                          'start,2100.3,type=canonical,ovl', text)
    if count != 1:
        raise ValueError('SCF seed requires one explicit canonical MULTI start.')
    text, count = re.subn(r'(?i)(\{rhf(?:,[^;{}]*)?;)(wf,)',
                          r'\1start,2100.3,type=canonical,ovl;\2', text, count=1)
    if count != 1:
        raise ValueError('SCF seed requires the generated RHF block.')
    return text
