"""Reuse an accepted relaxed-point energy only at the same physical level."""
from copy import deepcopy
import json
from pathlib import Path
import re

import numpy as np
from ase.io import read
from ase.units import Hartree

from kinbot.execution import file_sha256
from kinbot.molpro import vrc_multireference_energy_contract


def _model(value):
    value = deepcopy(value or {})
    # The accepted scan propagates its reviewed subspace. Check its electronic
    # definition separately from the evolving orbital checkpoint location.
    value.pop('orbital_selection', None)
    return json.loads(json.dumps(value))


def _same_subspace(model, wanted):
    """Follow hash-bound propagated records back to the requested selection."""
    from kinbot.anl import dispatch
    wanted_selection = (wanted or {}).get('orbital_selection')
    seen = set()
    for _ in range(256):
        if _model(model) != _model(wanted):
            return False
        selection = (model or {}).get('orbital_selection')
        if selection == wanted_selection:
            return True
        if not selection or selection.get('mode') != 'molpro_record':
            return False
        checkpoint = Path(selection['path']).resolve()
        if checkpoint in seen or checkpoint.parent.name != 'relaxed_scan' or checkpoint.parent.parent.name != 'tasks':
            return False
        seen.add(checkpoint)
        graph = checkpoint.parents[2]
        _, spec, state = dispatch._load(graph)
        entry = state['tasks'].get('relaxed_scan', {})
        execution = json.loads((checkpoint.parent/'execution.json').read_text())
        if (entry.get('status') != 'complete' or execution.get('status') != 'executed'
                or execution.get('artifacts', {}).get(checkpoint.name) != selection.get('sha256')):
            return False
        prior = next(t for t in spec['tasks'] if t['id'] == 'relaxed_scan')
        counts = re.search(r'\((\d+),(\d+)\)', prior['profile']['method'])
        if not counts:
            return False
        from ase.data import atomic_numbers
        electrons = sum(atomic_numbers[s] for s in spec['molecule']['symbols'])-spec['molecule']['charge']
        closed = (electrons-int(counts[1]))//2
        if (selection.get('record') != '2140.2' or selection.get('orbital_type') != 'natural'
                or selection.get('active_indices') != list(range(closed+1, closed+int(counts[2])+1))):
            return False
        model = prior['profile']['calculator_kwargs'].get('state_model')
    return False


def accepted_energy(row, component, species):
    """Return hash-bound source evidence, or None for an ordinary single point.

    No synthetic Molpro output is written. The row's legacy energy field is
    never used: older publishers attached relaxed energies to frozen rows too.
    """
    data = row.data
    if data.get('vrc_scan_role') not in ('relaxed', 'frozen') or not data.get('vrc_scan_dispatch'):
        return None
    graph = Path(data['vrc_scan_dispatch'])
    from kinbot.anl import dispatch
    _, spec, state = dispatch._load(graph)
    task = next((t for t in spec['tasks'] if t['id'] == 'relaxed_scan'), None)
    entry = state['tasks'].get('relaxed_scan', {})
    if task is None or entry.get('status') != 'complete':
        return None
    profile = task['profile']
    kwargs = profile.get('calculator_kwargs', {})
    if (profile.get('calculator') != 'molpro'
            or profile['method'].casefold() != component['method'].casefold()
            or profile['basis'].casefold() != component['basis'].casefold()
            or spec['molecule']['charge'] != int(species.charge)
            or spec['molecule']['multiplicity'] != (component.get('multiplicity') or int(species.mult))
            or kwargs.get('caspt2_shift', .1) != component.get('shift', .1)
            or kwargs.get('caspt2_ipea', .25) != component.get('ipea', .25)
            or not _same_subspace(kwargs.get('state_model'), component.get('state_model'))):
        return None
    directory = graph/'tasks'/'relaxed_scan'
    atoms = read(directory/'final.xyz')
    if (atoms.get_chemical_symbols() != list(species.atom)
            or not np.allclose(atoms.positions, species.geom, atol=5.e-9, rtol=0)
            or dispatch._geometry_hash(atoms) != entry.get('final_geometry_sha256')):
        return None
    dispatch._verify_stage_files(graph, task, entry)
    execution = json.loads((directory/'execution.json').read_text())
    if execution.get('status') != 'executed' or execution['geometry_sha256'] != entry['geometry_sha256']:
        return None
    files = execution.get('details', {}).get('generated_files', [])
    inputs = [n for n in files if re.fullmatch(r'relaxed_scan_step_\d+\.inp', n)]
    if not inputs:
        return None
    inp = directory/max(inputs, key=lambda n: int(re.search(r'_(\d+)\.inp$', n)[1]))
    out = inp.with_suffix('.out')
    # Bounded reads: only the final central evaluation, not every displacement
    # and multi-megabyte orbital checkpoint from the entire optimization.
    hashes = {}
    for p in (inp, out, directory/'final.xyz'):
        hashes[p.name] = file_sha256(p)
        if execution.get('artifacts', {}).get(p.name) != hashes[p.name]:
            raise ValueError(f'Accepted scan energy source changed: {p}')
    block = re.search(r'(?ims)^geometry=\{\s*\n(.*?)\n\}', inp.read_text())
    if block is None:
        return None
    from io import StringIO
    central = read(StringIO(block[1]+'\n'), format='xyz')
    if (central.get_chemical_symbols() != atoms.get_chemical_symbols()
            or not np.allclose(central.positions, atoms.positions, atol=5.e-9, rtol=0)):
        return None
    text = out.read_text(errors='replace')
    matches = re.findall(r'^\s*KB_GEOM_ENERGY\s*=\s*([-+\d.EeDd]+)\s+AU', text, re.M)
    energy = execution['details']['energy_ev']/Hartree
    if (not matches or 'Molpro calculation terminated' not in text
            or not np.isfinite(energy)
            or abs(float(matches[-1].replace('D', 'E'))-energy) > 1.e-9):
        raise ValueError('Final relaxed energy and native central evaluation disagree.')
    return dict(source_kind='accepted_relaxed_scan', source_graph=str(graph.resolve()),
        source_input=str(inp.resolve()), source_output=str(out.resolve()),
        input_sha256=hashes[inp.name], output_sha256=hashes[out.name],
        geometry_sha256=entry['final_geometry_sha256'],
        execution_sha256=file_sha256(directory/'execution.json'), energy_hartree=energy,
        energy_contract=vrc_multireference_energy_contract(component['method']))
