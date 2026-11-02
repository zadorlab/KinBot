"""Refresh cached population labels/declared levels from their own geometry.

No energies, geometries, frequencies or reactive-surface spin weights are
changed here. The PES energy reader applies the same declared scalar shift.
"""
import re
import numpy as np

from kinbot.mess_networks import split_model, _HEADER


def _refresh_segment(text, header, parameters, selected_states):
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.species_routing import connectivity_name
    from kinbot.reference_identity import same_atomic_connectivity
    from kinbot.electronic_states import model_for_species

    fields = header[2].split('!')
    if len(fields) < 2:
        return text, None
    identity = fields[1].strip()
    key = connectivity_name(identity)
    if not key.isdigit():
        return text, None
    geom = re.search(r'Geometry\[angstrom\]\s+(\d+)\s*\n', text)
    atom = re.search(r'Atom\s+Name\s+([A-Z][a-z]?)\b', text)
    if geom:
        rows = [line.split() for line in text[geom.end():].splitlines()[:int(geom[1])]]
        if any(len(row) != 4 for row in rows):
            raise ValueError('Incomplete cached MESS geometry.')
        symbols, positions = [r[0] for r in rows], np.array([r[1:] for r in rows], float)
    elif atom:
        symbols, positions = [atom[1]], np.zeros((1, 3))
    else:
        return text, None
    previous = fields[-1].strip() if len(fields) > 2 else None
    state = (selected_states or {}).get(identity)
    charge = None if state is None else state['charge']
    if previous:
        from rdkit import Chem
        molecule = Chem.MolFromSmiles(previous)
        if molecule is not None and charge is None:
            charge = Chem.GetFormalCharge(molecule)
    if charge is None:
        # A bare atomic symbol/chemid carries no charge. Do not infer neutral
        # fragments from a neutral parent; ion-pair exits are possible.
        return text, None
    point = StationaryPoint(identity, charge, int(key[-1]),
                            atom=symbols, geom=positions)
    point.characterize()
    # Preserve the configured Lewis representation if it has the same graph.
    # Reperceived labels cannot change the isotope/spin/graph routing identity.
    if previous and same_atomic_connectivity(previous, point.smiles):
        point.smiles = previous
    if previous != point.smiles:
        replacement = '!'.join(fields[:2]).rstrip() + ' ! ' + point.smiles
        text = text[:header.start(2)] + replacement + text[header.end(2):]
    model = model_for_species(point, parameters)
    if model is None or not model['levels']:
        return text, point.smiles
    if point.mult != model['ground_multiplicity']:
        raise ValueError(f'{identity}: cached QC spin differs from its declared electronic model.')
    pattern = re.compile(r'(ElectronicLevels\[1/cm\]\s+)\d+[^\S\n]*\n'
                         r'(?:[ \t]+[-+\d.eE]+[ \t]+\d+[^\S\n]*\n)+')
    levels = model['levels']
    text, count = pattern.subn(lambda m: m[1] + str(len(levels)) + '\n' + ''.join(
        f'          {v["energy_cm_inverse"]:.12g} {v["degeneracy"]}\n' for v in levels), text)
    if not count:
        raise ValueError(f'{identity}: cached population has no electronic levels.')
    return text, point.smiles


def refresh_cached_population_models(contents, parameters, *, selected_states=None):
    """Rebuild display labels and configured levels, retaining all raw caches."""
    prefix, blocks, end = split_model(contents)
    updated = []
    fragment = re.compile(r'^([ \t]*Fragment)[ \t]+([^\n]+)\n', re.M)
    for block in blocks:
        header = _HEADER.search(block)
        if header[1] == 'Well':
            block, _ = _refresh_segment(block, header, parameters, selected_states)
        elif header[1] == 'Bimolecular':
            matches = list(fragment.finditer(block))
            names = []
            for index in reversed(range(len(matches))):
                start = matches[index].start()
                finish = matches[index+1].start() if index+1 < len(matches) else len(block)
                segment = block[start:finish]
                segment, name = _refresh_segment(segment, fragment.search(segment), parameters, selected_states)
                names.insert(0, name)
                block = block[:start]+segment+block[finish:]
            if names and all(names):
                header = _HEADER.search(block)
                fields = header[2].split('!')
                if len(fields) >= 2:
                    value = '!'.join(fields[:2]).rstrip()+' ! '+' '.join(names)
                    block = block[:header.start(2)]+value+block[header.end(2):]
        updated.append(block)
    return prefix+''.join(updated)+end
