"""Validate configuration evidence for reuse and preserve unsupported observations.

Configured names and verified read aliases live in species_routing.
"""
import copy
import json
import os
from pathlib import Path
import tempfile
import numpy as np
from ase import Atoms
from kinbot.stereo_identity import canonical_identity
from kinbot import constants


class StereoRoutingError(ValueError):
    pass


def _json_value(value):
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    return value


def _chemical_context(species):
    """Keep supplied labels/graphs even when canonical assignment is unsupported."""
    return _json_value({key: getattr(species, key, None) for key in (
        'isotopes', 'formal_charges', 'bond', 'bond01', 'bonds', 'rads',
        'optical_reference', 'optical_population', 'optical_counting_scope')})


def _raw_input(species):
    """Exact requested QC input, not a chemical/configurational equivalence key.

    Optical bookkeeping may be attached later; it is saved in the evidence but
    does not alter the coordinates or atom labels requested for this QC job.
    """
    context = _chemical_context(species)
    return {'atoms': list(map(str, species.atom)),
            'geometry_angstrom': np.asarray(species.geom).tolist(),
            'charge': int(species.charge), 'multiplicity': int(species.mult),
            'chemical_labels': {key: value for key, value in context.items()
                                if not key.startswith('optical_')}}


def _write_reference(qc, requested, name, identity):
    qc.db.write(Atoms(requested.atom, positions=requested.geom), name=name,
                data={'status': 'input_reference', 'identity': identity,
                      'input_charge': requested.charge, 'input_multiplicity': requested.mult,
                      'chemical_context': _chemical_context(requested),
                      'raw_input': _raw_input(requested)})


def preserve_observations(reason, species_list):
    """Write a unique evidence report before refusing unsupported counting/routing."""
    records = []
    for species in species_list:
        records.append({'name': str(getattr(species, 'name', 'observation')),
            'chemid': str(getattr(species, 'chemid', 'unknown')),
            'charge': int(getattr(species, 'charge', 0)),
            'multiplicity': int(getattr(species, 'mult', 1)),
            'atoms': list(map(str, species.atom)),
            'geometry_angstrom': np.asarray(species.geom).tolist(),
            'source_job': getattr(species, 'source_job', None),
            'source_row_id': getattr(species, 'source_row_id', None),
            'electronic_energy_hartree': getattr(species, 'energy', None),
            'zpe_hartree': getattr(species, 'zpe', None),
            'raw_frequencies_cm-1': getattr(species, 'freq', []),
            'identity': canonical_identity(species),
            'chemical_context': _chemical_context(species),
            'calculation_state_evidence': getattr(species, 'calculation_state_evidence', None),
            'conformer_inventory': [r.as_dict() for r in getattr(species, 'conformer_inventory', ())]})
    directory = Path('unsupported_observations')
    directory.mkdir(exist_ok=True)
    fd, filename = tempfile.mkstemp(prefix='observation_', suffix='.json', dir=directory)
    with os.fdopen(fd, 'w') as handle:
        json.dump(_json_value({'schema_version': 1, 'reason': reason,
            'rate_model_ready': False, 'observations': records}), handle, indent=2, allow_nan=False)
        handle.write('\n')
    return str(Path(filename).resolve())


def refuse_routing(reason, species_list):
    path = preserve_observations(reason, species_list)
    raise StereoRoutingError(f'{reason} Observations preserved in {path}. '
                             'This configured population/routing model is unsupported.')


def require_same_configuration(first, second, context):
    """Guard a proposed legacy alias, preserving no-RDKit non-MC compatibility."""
    if first is second:
        return
    left, right = canonical_identity(first), canonical_identity(second)
    if 'unavailable' in (left['status'], right['status']):
        first.stereo_routing_status = second.stereo_routing_status = 'legacy; RDKit unavailable; unverified'
        return
    if left['status'] == right['status'] == 'assigned' and left['id'] == right['id']:
        return
    if _in_declared_mirror_population(first, right):
        return
    if left['status'] != 'assigned' and right['status'] != 'assigned':
        from kinbot.species_routing import same_chemical_graph
        from kinbot.stereo_identity import legacy_stereo_warning
        if same_chemical_graph(first, second):
            legacy_stereo_warning(first, left.get('reason'))
            legacy_stereo_warning(second, right.get('reason'))
            return
    refuse_routing(f'{context}: legacy identity cannot establish the same configured species', [first, second])


def _in_declared_mirror_population(species, identity):
    reference = getattr(species, 'optical_reference', {}) or {}
    return (getattr(species, 'optical_population', 'specified') == 'racemic'
            and reference.get('status') == identity.get('status') == 'assigned'
            and identity.get('id') in {reference['id'], reference['mirror_id']}
            and canonical_identity(species).get('id') in {reference['id'], reference['mirror_id']})


def _observed_species(species, atoms, geom):
    from kinbot.stationary_pt import StationaryPoint
    observed = StationaryPoint('cached observation', species.charge, species.mult, atom=atoms, geom=geom)
    observed.bond_mx()
    observed.calc_chemid()
    return observed


def _row_species(species, row, job, input_reference=None):
    reference = copy.copy(species)
    state = {}
    containers = [('data', row.data), ('identity', row.data.get('identity', {})),
                  ('calculator_parameters', getattr(row, 'calculator_parameters', {}) or {})]
    for attribute, fields in (('charge', ('charge', 'input_charge')),
                              ('mult', ('mult', 'multiplicity', 'input_multiplicity'))):
        observations = [{'source': source + '.' + field, 'value': container[field]}
                        for source, container in containers for field in fields
                        if container.get(field) is not None]
        fallback = getattr(species, attribute)
        source = 'caller assumption; result state unreported'
        if input_reference is not None:
            field = 'input_charge' if attribute == 'charge' else 'input_multiplicity'
            fallback = input_reference.data.get(field, fallback)
            source = 'trusted requested input; result state unreported'
        value = observations[0]['value'] if observations else fallback
        setattr(reference, attribute, int(value))
        state[attribute] = {'value_used': value,
                            'source': observations[0]['source'] if observations else source,
                            'observations': observations}
    point = _observed_species(reference, row.symbols, row.positions)
    point.calculation_state_evidence = state
    if input_reference is not None and list(row.symbols) == list(input_reference.symbols):
        # A trusted requested-input record supplies fixed atom labels absent
        # from ordinary ASE result rows. Rebuild the optimized graph from the
        # result geometry; do not overwrite it with the input bond matrix.
        supplied = input_reference.data.get('chemical_context', {})
        for key in ('isotopes', 'formal_charges', 'optical_reference',
                    'optical_population', 'optical_counting_scope'):
            if supplied.get(key) is not None:
                setattr(point, key, copy.deepcopy(supplied[key]))
    for key, value in row.data.get('chemical_context', {}).items():
        if value is not None:
            if key in ('bond', 'bond01'):
                value = np.asarray(value)
            elif key in ('bonds', 'rads'):
                value = [np.asarray(item) for item in value]
            setattr(point, key, copy.deepcopy(value))
    point.source_job, point.source_row_id = job, row.id
    energy = row.data.get('energy')
    point.energy = energy * constants.EVtoHARTREE if energy is not None else None
    point.zpe = row.data.get('zpe')
    point.freq = row.data.get('frequencies', [])
    return point


def guard_well_job(qc, species, geom, job):
    """Check both input-reference and completed cache geometry before job reuse.

    Ordinary optimized fragments may dissociate; changed connectivity remains
    handled by the existing product workflow. Same-chemid stereo changes must
    not be mistaken for the requested configured species.
    """
    requested = copy.copy(species)
    requested.geom = np.asarray(geom).copy()
    requested.source_job = job
    identity = canonical_identity(requested)
    if identity['status'] == 'unavailable':
        species.stereo_routing_status = 'legacy; RDKit unavailable; unverified'
        return
    declared = getattr(requested, 'optical_reference', {}) or {}
    if (declared.get('status') == identity.get('status') == 'assigned'
            and declared['id'] != identity['id']
            and not _in_declared_mirror_population(requested, identity)):
        refuse_routing(f'{job}: requested geometry is outside its declared configuration', [requested])
    name = f'stereochemistry/{job}'
    references = list(qc.db.select(name=name))
    rows = list(qc.db.select(name=job))
    cached = None
    if rows and rows[-1].data.get('status') == 'normal':
        cached = _row_species(species, rows[-1], job, references[-1] if references else None)
        if any(observation['value'] != getattr(requested, attribute)
               for attribute, evidence in cached.calculation_state_evidence.items()
               for observation in evidence['observations']):
            refuse_routing(f'{job}: completed cache reports a contradictory charge or multiplicity',
                           [requested, cached])
    if (identity['status'] != 'assigned'
            or references and references[-1].data['identity']['status'] != 'assigned'):
        from kinbot.stereo_identity import legacy_stereo_warning
        from kinbot.species_routing import same_chemical_graph
        legacy_stereo_warning(species, identity.get('reason'))
        if references:
            previous = _row_species(species, references[-1], job)
            previous_identity = references[-1].data['identity']
            if ((identity['status'] == 'assigned') !=
                    (previous_identity['status'] == 'assigned')):
                refuse_routing(f'{job}: unsupported identity cannot replace a known configured input',
                               [requested, previous])
            if not same_chemical_graph(requested, previous):
                refuse_routing(f'{job}: saved input has different chemical labels or connectivity',
                               [requested, previous])
        elif cached is not None:
            if not same_chemical_graph(requested, cached):
                refuse_routing(f'{job}: saved result has different chemical labels or connectivity',
                               [requested, cached])
            _write_reference(qc, requested, name, identity)
        else:
            _write_reference(qc, requested, name, identity)
        return
    if references:
        reference = references[-1]
        from kinbot.species_routing import require_cache_atom_order
        require_cache_atom_order(requested, _row_species(species, reference, job), job)
        if (reference.data['identity']['id'] != identity['id']
                and not _in_declared_mirror_population(requested, reference.data['identity'])):
            previous = _row_species(species, reference, job)
            refuse_routing(f'{job}: cached input belongs to a different configured species', [requested, previous])
    if cached is not None:
        from kinbot.species_routing import require_cache_atom_order
        require_cache_atom_order(requested, cached, job)
        if cached.chemid == species.chemid:
            require_same_configuration(requested, cached, f'{job}: completed cache')
    if not references:
        _write_reference(qc, requested, name, identity)
    species.stereo_routing_status = 'verified configured request'


def guard_pes_input(species, filename):
    """Prevent a same-chemid PES input from being overwritten by another isomer."""
    path = Path(filename)
    if not path.exists():
        return
    from kinbot.stationary_pt import StationaryPoint
    data = json.loads(path.read_text())
    previous = StationaryPoint('existing PES input', data.get('charge', 0), data.get('mult', 1),
                               structure=data.get('structure'), smiles=data.get('smiles') or None)
    previous.bond_mx()
    previous.calc_chemid()
    require_same_configuration(species, previous, f'{filename}: PES input collision')
