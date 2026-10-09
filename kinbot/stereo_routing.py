"""Validate current calculation inputs and results before reuse."""
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import numpy as np
from ase import Atoms
from kinbot.stereo_identity import canonical_identity, require_supported_identity
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
    return qc.db.write(Atoms(requested.atom, positions=requested.geom), name=name,
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
    """Check current input/result identity, including a declared mirror pair."""
    left, right = require_supported_identity(first), require_supported_identity(second)
    if left['id'] == right['id'] or _in_declared_mirror_population(first, right):
        return
    refuse_routing(f'{context}: result belongs to a different stereoisomer', [first, second])


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


def guard_well_job(qc, species, geom, job, initial_product=False):
    """Check both input-reference and completed cache geometry before job reuse.

    A newly formed product is named before its first optimization. Its final
    connectivity and stereochemistry are identified by the product workflow;
    they need not agree with that preliminary name. Other well calculations
    retain the requested stereoisomer checks.
    """
    requested = copy.copy(species)
    requested.geom = np.asarray(geom).copy()
    requested.source_job = job
    name = f'stereochemistry/{job}'
    reference = next(qc.db.select(name=name, sort='-id', limit=1), None)
    result = next(qc.db.select(name=job, sort='-id', limit=1), None)
    request = hashlib.sha256(json.dumps(_json_value({
        'input': _raw_input(requested), 'chemid': requested.chemid,
        'wellorts': getattr(requested, 'wellorts', 0),
        'initial_product': bool(initial_product),
        **{key: getattr(requested, key, None) for key in (
            'optical_reference', 'optical_population',
            'stereo_ignored_atoms', 'stereo_ignored_bonds')}
    }), sort_keys=True).encode()).hexdigest()
    cache = getattr(qc, '_verified_well_jobs', {})
    state = (request, _row_revision(reference), _row_revision(result))
    previous_check = cache.get(job)
    if previous_check is not None and previous_check[0] is qc.db and previous_check[1] == state:
        species.stereo_routing_status = 'verified configured request'
        return
    identity = require_supported_identity(requested)
    declared = getattr(requested, 'optical_reference', {}) or {}
    if (declared.get('status') == identity.get('status') == 'assigned'
            and declared['id'] != identity['id']
            and not _in_declared_mirror_population(requested, identity)):
        refuse_routing(f'{job}: requested geometry is outside its declared configuration', [requested])
    cached = None
    if result is not None and result.data.get('status') == 'normal':
        cached = _row_species(species, result, job, reference)
        if any(observation['value'] != getattr(requested, attribute)
               for attribute, evidence in cached.calculation_state_evidence.items()
               for observation in evidence['observations']):
            refuse_routing(f'{job}: completed cache reports a contradictory charge or multiplicity',
                           [requested, cached])
    if reference is not None:
        from kinbot.species_routing import require_cache_atom_order
        previous = _row_species(species, reference, job)
        require_cache_atom_order(requested, previous, job)
        saved_identity = (reference.data.get('chemical_context', {}).get('optical_reference')
                          or reference.data['identity'])
        if saved_identity['id'] != (declared or identity)['id']:
            refuse_routing(f'{job}: cached input has a different full stereoisomer identity',
                           [requested, previous])
        if (reference.data['identity']['id'] != identity['id']
                and not _in_declared_mirror_population(requested, reference.data['identity'])):
            refuse_routing(f'{job}: cached input belongs to a different configured species', [requested, previous])
    if cached is not None:
        from kinbot.species_routing import require_cache_atom_order
        if initial_product:
            if list(requested.atom) != list(cached.atom):
                refuse_routing(f'{job}: product result changes the supplied atom order',
                               [requested, cached])
        else:
            require_cache_atom_order(requested, cached, job)
        if not initial_product and cached.chemid == species.chemid:
            require_same_configuration(requested, cached, f'{job}: completed cache')
    if reference is None:
        reference = qc.db.get(_write_reference(qc, requested, name, identity))
    # ASE writes and updates change the row revision. A newly completed result,
    # changed request, or changed reference must pass the full checks again.
    cache[job] = (qc.db, (request, _row_revision(reference), _row_revision(result)))
    qc._verified_well_jobs = cache
    species.stereo_routing_status = 'verified configured request'


def _row_revision(row):
    """Track both appended results and in-place ASE database updates."""
    return None if row is None else (row.id, row.unique_id, row.mtime, row.data.get('status'))


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
    from kinbot.species_routing import apply_input_reference
    apply_input_reference(previous, data)
    requested_identity = getattr(species, 'optical_reference', None) or require_supported_identity(species)
    saved_identity = getattr(previous, 'optical_reference', None) or require_supported_identity(previous)
    if requested_identity['id'] != saved_identity['id']:
        refuse_routing(f'{filename}: PES input has a different full stereoisomer identity',
                       [species, previous])
    require_same_configuration(species, previous, f'{filename}: PES input collision')
