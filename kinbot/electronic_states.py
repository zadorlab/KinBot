"""Validated molecular spin states shared by products, VRC, and MESS."""

from copy import deepcopy
from numbers import Real

import numpy as np


def _canonical_smiles(value):
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(str(value))
    if molecule is None:
        raise ValueError(f'Invalid electronic-state SMILES: {value!r}.')
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def electronic_state_models(parameters):
    """Return validated species electronic states from a user parameter map."""
    records = parameters.get('electronic_state_models', [])
    if records is None:
        records = []
    if not isinstance(records, list):
        raise ValueError('electronic_state_models must be a list.')
    result = {}
    for record in records:
        allowed = {'smiles', 'charge', 'ground_multiplicity', 'levels',
                   'source', 'smiles_aliases', 'term_orbital_degeneracy',
                   'spin_orbit_correction_hartree', 'spin_free_reference'}
        if (not isinstance(record, dict) or set(record) - allowed
                or not {'smiles', 'charge', 'ground_multiplicity'}
                <= set(record)):
            raise ValueError('Each electronic-state model requires smiles, '
                             'charge, ground_multiplicity, and level or term evidence.')
        charge = record['charge']
        ground = record['ground_multiplicity']
        if (isinstance(charge, bool) or not isinstance(charge, int)
                or isinstance(ground, bool) or not isinstance(ground, int)
                or ground < 1):
            raise ValueError('Electronic-state charge and multiplicity must '
                             'be integers, with positive multiplicity.')
        levels = record.get('levels', [])
        if not isinstance(levels, list) or not levels and 'term_orbital_degeneracy' not in record:
            raise ValueError('Electronic-state models need levels or a declared spatial term degeneracy.')
        normalized_levels = []
        for level in levels:
            allowed_level = {'energy_cm_inverse', 'multiplicity',
                             'orbital_degeneracy', 'degeneracy', 'label'}
            if (not isinstance(level, dict) or set(level) - allowed_level
                    or not {'energy_cm_inverse', 'multiplicity'} <= set(level)):
                raise ValueError('Each electronic level requires energy and '
                                 'multiplicity.')
            energy = level['energy_cm_inverse']
            multiplicity = level['multiplicity']
            if (isinstance(energy, bool) or not isinstance(energy, Real)
                    or not np.isfinite(energy) or energy < 0.
                    or isinstance(multiplicity, bool)
                    or not isinstance(multiplicity, int)
                    or multiplicity < 1):
                raise ValueError('Electronic-level energies must be finite '
                                 'and nonnegative; multiplicities must be '
                                 'positive integers.')
            orbital = level.get('orbital_degeneracy', 1)
            supplied = level.get('degeneracy')
            if (isinstance(orbital, bool) or not isinstance(orbital, int)
                    or orbital < 1):
                raise ValueError('Orbital degeneracy must be a positive integer.')
            degeneracy = multiplicity * orbital
            if supplied is not None:
                if (isinstance(supplied, bool) or not isinstance(supplied, int)
                        or supplied < 1 or ('orbital_degeneracy' in level
                                          and supplied != degeneracy)):
                    raise ValueError('Explicit electronic degeneracy must equal '
                                     'multiplicity times orbital_degeneracy.')
                # A resolved fine-structure level has degeneracy 2J+1,
                # which need not be divisible by the spin multiplicity.
                degeneracy = supplied
            normalized_levels.append({
                'energy_cm_inverse': float(energy),
                'multiplicity': multiplicity,
                **({'orbital_degeneracy': orbital}
                   if supplied is None or 'orbital_degeneracy' in level else {}),
                'degeneracy': degeneracy,
                **({'label': str(level['label'])} if 'label' in level else {}),
            })
        if normalized_levels and normalized_levels[0]['energy_cm_inverse'] != 0.:
            raise ValueError('The first electronic level must be at zero energy.')
        if any(right['energy_cm_inverse'] < left['energy_cm_inverse']
               for left, right in zip(normalized_levels,
                                      normalized_levels[1:])):
            raise ValueError('Electronic levels must be ordered by energy.')
        if normalized_levels and normalized_levels[0]['multiplicity'] != ground:
            raise ValueError('The zero-energy level must have the declared '
                             'ground-state multiplicity.')
        canonical = _canonical_smiles(record['smiles'])
        from rdkit import Chem
        molecule = Chem.AddHs(Chem.MolFromSmiles(canonical))
        electrons = sum(atom.GetAtomicNum() for atom in molecule.GetAtoms()) - charge
        if (Chem.GetFormalCharge(molecule) != charge
                or ground - 1 > electrons
                or (electrons - ground + 1) % 2
                or any((electrons - level['multiplicity'] + 1) % 2
                       or level['multiplicity'] - 1 > electrons
                       for level in normalized_levels)):
            raise ValueError('Electronic-state charge or spin disagrees with '
                             'the molecular electron count.')
        aliases = record.get('smiles_aliases', [])
        if not isinstance(aliases, list):
            raise ValueError('smiles_aliases must be a list.')
        canonical_aliases = [_canonical_smiles(alias) for alias in aliases]
        # An alias can supply an open-shell Lewis representation for the
        # same connectivity (e.g. [S][S] for KinBot's S=S). Never match merely
        # by formula or silently alias a different molecular skeleton.
        def skeleton(smiles):
            import networkx as nx
            mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
            graph = nx.Graph()
            for atom in mol.GetAtoms():
                graph.add_node(atom.GetIdx(), element=atom.GetAtomicNum(),
                               charge=atom.GetFormalCharge(),
                               isotope=atom.GetIsotope())
            graph.add_edges_from((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                                 for bond in mol.GetBonds())
            return graph
        import networkx as nx
        graph = skeleton(canonical)
        for alias in canonical_aliases:
            if not nx.is_isomorphic(graph, skeleton(alias),
                                    node_match=lambda left, right: left == right):
                raise ValueError('Electronic-state SMILES aliases must have '
                                 'the same atoms, charge, and connectivity.')
        model = {
            'smiles': canonical, 'charge': charge,
            'ground_multiplicity': ground,
            'levels': normalized_levels,
            'source': str(record.get('source', 'user-specified')),
        }
        term = record.get('term_orbital_degeneracy')
        if term is not None:
            if isinstance(term, bool) or not isinstance(term, int) or term < 1:
                raise ValueError('term_orbital_degeneracy must be a positive integer.')
            model['term_orbital_degeneracy'] = term
        correction = record.get('spin_orbit_correction_hartree')
        convention = record.get('spin_free_reference')
        if convention is not None:
            if convention != 'term_centroid' or correction is not None:
                raise ValueError('Use term_centroid or an explicit SO correction, not both.')
            if (term is None or any(level['multiplicity'] != ground for level in normalized_levels)
                    or sum(level['degeneracy'] for level in normalized_levels) != ground * term):
                raise ValueError('A spin-free term centroid requires a complete ground-term level set.')
            from ase.units import invcm, Hartree
            correction = -sum(level['energy_cm_inverse'] * level['degeneracy']
                              for level in normalized_levels) / (ground * term) * invcm / Hartree
        if correction is not None:
            if isinstance(correction, bool) or not isinstance(correction, Real) or not np.isfinite(correction):
                raise ValueError('SO correction must be a finite energy in hartree.')
            model['spin_orbit_correction'] = {
                'value_hartree': float(correction), 'source': model['source'],
                'reference': convention or 'user-specified spin-free reference',
            }
        for smiles in [canonical, *canonical_aliases]:
            key = (smiles, charge)
            if key in result:
                raise ValueError('Duplicate electronic-state species declaration.')
            result[key] = model
    return result


def model_for_species(species, parameters):
    """Return the configured electronic model matching one species."""
    models = electronic_state_models(parameters)
    if not models:
        return None
    smiles = getattr(species, 'smiles', None)
    if not smiles:
        return None
    try:
        canonical = _canonical_smiles(smiles)
    except ValueError:
        # Some legacy KinBot structures carry a connectivity string which is
        # useful internally but is not a sanitizable RDKit SMILES.  An
        # electronic-state model is optional for those species, and its
        # absence must retain the historical multiplicity-only MESS block.
        return None
    return models.get((canonical, int(getattr(species, 'charge', 0))))


def audit_electronic_states(states, parameters):
    """Inventory supplied states without inferring SO data from composition.

    This is a postprocessing coverage report, not an electronic-state solver.
    Missing open-shell models remain explicitly unreviewed. Molecular/atomic
    distinctions follow SMILES, so C (methane) cannot stand in for [C].
    """
    models = electronic_state_models(parameters)
    rows = {}
    for state in states:
        smiles = _canonical_smiles(state['smiles'])
        charge, mult = int(state['charge']), int(state['multiplicity'])
        key = (smiles, charge, mult)
        row = rows.setdefault(key, dict(smiles=smiles, charge=charge,
                                       multiplicity=mult, roles=[]))
        role = state.get('role', 'unspecified')
        if role not in row['roles']:
            row['roles'].append(role)
        model = models.get((smiles, charge))
        if model is None:
            row['status'] = 'singlet_no_first_order_so' if mult == 1 else 'unreviewed_open_shell'
            continue
        if model['ground_multiplicity'] != mult:
            row['status'] = 'configured_multiplicity_mismatch'
        elif 'spin_orbit_correction' in model:
            row['status'] = 'configured_correction'
            row['spin_orbit_correction'] = deepcopy(model['spin_orbit_correction'])
        elif (model.get('term_orbital_degeneracy') == 1
              and len(model['levels']) == 1):
            row['status'] = 'declared_unsplit_ground_term'
        else:
            row['status'] = 'levels_only_no_energy_correction'
        row['levels'] = deepcopy(model['levels'])
        row['source'] = model['source']
    entries = [rows[key] for key in sorted(rows)]
    return {'schema': 1, 'scope': 'first-order ground-term SO and configured electronic levels',
            'states': entries,
            'requires_review': [dict(smiles=r['smiles'], charge=r['charge'],
                                     multiplicity=r['multiplicity'], status=r['status'])
                                for r in entries if r['status'] in (
                                    'unreviewed_open_shell', 'configured_multiplicity_mismatch')],
            'note': 'Coverage is not a guarantee of thermochemical accuracy; '
                    'no new correction or excited-state ladder is inferred.'}


def apply_electronic_state_model(species, parameters):
    """Apply an explicit ground state before QC job names are constructed."""
    model = model_for_species(species, parameters)
    if model is None:
        return None
    species.mult = model['ground_multiplicity']
    species.smiles = model['smiles']
    species.electronic_levels = deepcopy(model['levels'])
    species.electronic_state_source = model['source']
    if 'spin_orbit_correction' in model:
        species.spin_orbit_correction = deepcopy(model['spin_orbit_correction'])
    elif hasattr(species, 'spin_orbit_correction'):
        del species.spin_orbit_correction
    if hasattr(species, 'calc_chemid'):
        species.calc_chemid()
    return model


def state_average_roots(products, multiplicities, parameters=None):
    """Derive spatial root counts from declared separated electronic terms."""
    spatial_roots = 1
    for product in products:
        model = model_for_species(product, parameters or {})
        if model is not None:
            ground = model['levels'][0] if model['levels'] else {}
            spatial_roots *= model.get(
                'term_orbital_degeneracy', ground.get('orbital_degeneracy', 1))
    return {int(multiplicity): spatial_roots for multiplicity in multiplicities}


def mess_electronic_levels(species, parameters, *, ground_only=False):
    """Return MESS level count and lines for a species or reactive surface."""
    levels = electronic_levels_for_species(
        species, parameters, ground_only=ground_only)
    lines = '\n'.join(
        f'{level["energy_cm_inverse"]:.12g} {level["degeneracy"]}'
        for level in levels)
    return len(levels), lines


def electronic_levels_for_species(species, parameters, *, ground_only=False):
    """Return normalized MESS electronic levels for one molecular state.

    Stable species can carry several measured or calculated levels.  A
    transition-state or VRC dividing surface requests ``ground_only`` so its
    reactive-surface multiplicity is counted once while the separated
    fragment levels remain in the reactant partition functions.
    """
    model = None if ground_only else model_for_species(species, parameters)
    if model is not None and int(species.mult) != model['ground_multiplicity']:
        raise ValueError('Species multiplicity disagrees with its electronic-state model.')
    if model is None or not model['levels']:
        levels = [{'energy_cm_inverse': 0.,
                   'degeneracy': int(species.mult)}]
    else:
        levels = model['levels']
    return deepcopy(levels)
