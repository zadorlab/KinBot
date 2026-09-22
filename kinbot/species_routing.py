"""Configured species names layered on top of the existing connectivity ID.

Ordinary names are unchanged. A stereo suffix contains the existing canonical
key, with no underscore so the PES fragment/reaction delimiters remain valid.
Verified legacy QC names are read aliases, not new chemical identities.
"""
import json
import hashlib
import copy
import logging
from pathlib import Path
import re
from ase.db import connect as ase_connect
from kinbot.stereo_identity import canonical_identity


_KEY = re.compile(r'^(\d+)(?:-s([0-9a-f]{64}))?$')
_JOB = re.compile(r'(^|/)(\d+-s[0-9a-f]{64})(?=_|$)')


def routing_key(species):
    """Return the legacy key or a configured key; never change ``chemid``."""
    if getattr(species, 'wellorts', 0):
        return species.name
    identity = getattr(species, 'optical_reference', None)
    if identity is None:
        identity = canonical_identity(species)
    if identity.get('status') == 'assigned' and any(
            any(tag in graph for tag in ('@', '/', '\\'))
            for graph in identity.get('canonical_graphs', ())):
        return f"{species.chemid}-s{identity['id']}"
    return species.chemid


def routing_name(species):
    return str(routing_key(species))


def mess_filename(key, iteration):
    """Keep complete network keys in content without overlong file components."""
    stem = str(key)
    if len(stem.encode()) > 220:
        stem = 'network-' + hashlib.sha256(stem.encode()).hexdigest()
    return f'{stem}_{int(iteration):04d}.mess'


def require_canonical_summary(job):
    """An old PES serialization must be regenerated, never silently dropped."""
    root = Path(job)
    old = connectivity_name(root.name)
    if (old != root.name and not (root / f'summary_{root.name}.out').exists()
            and (root / f'summary_{old}.out').exists()):
        raise ValueError(f'{job}: legacy PES summary needs configured output regeneration. '
                         'Run a normal cached KinBot restart before no-kinbot postprocessing; '
                         'the original summary and calculation files are retained.')


def connectivity_name(name):
    match = _KEY.fullmatch(str(name))
    return match.group(1) if match else str(name)


def is_species_name(name):
    return _KEY.fullmatch(str(name)) is not None


def matches_name(name, choices):
    """Honor legacy connectivity selectors as well as exact configured names."""
    legacy = _JOB.sub(lambda match: match.group(1) + match.group(2).split('-s')[0], str(name))
    return str(name) in choices or legacy in choices


def configured_selection(mapping, species, default=None):
    """Prefer the configured input selector, retaining connectivity defaults."""
    key = routing_name(species)
    return mapping[key] if key in mapping else mapping.get(str(species.chemid), default)


def same_species(first, second, context='species reuse'):
    """Distinct supported configurations are independent objects, not errors."""
    if first.chemid != second.chemid:
        return False
    left, right = canonical_identity(first), canonical_identity(second)
    if left['status'] == right['status'] == 'assigned':
        if routing_key(first) != routing_key(second):
            return False
        return (left['id'] == right['id'] or
                getattr(first, 'optical_population', 'specified') ==
                getattr(second, 'optical_population', 'specified') == 'racemic'
                and left['mirror_id'] == right['id'])
    # Unknown stereo does not erase a known configured distinction. Ordinary
    # unsupported graphs can still use the legacy treatment, but chemid alone
    # is insufficient because different molecular graphs can share that ID.
    if (left['status'] == 'assigned') != (right['status'] == 'assigned'):
        return False
    from kinbot.stereo_identity import legacy_stereo_warning
    legacy_stereo_warning(first, left.get('reason'))
    legacy_stereo_warning(second, right.get('reason'))
    return same_chemical_graph(first, second)


def same_chemical_graph(first, second):
    """Compare labelled connectivity without assuming unsupported stereo."""
    if (first.charge, first.mult) != (second.charge, second.mult):
        return False
    import networkx as nx
    from kinbot.molecular_symmetry import chemical_graph
    return nx.is_isomorphic(chemical_graph(first), chemical_graph(second),
                            node_match=lambda a, b: a['label'] == b['label'],
                            edge_match=lambda a, b: a['label'] == b['label'])


def reusable_cached_product(qc, species):
    """Adopt a verified product's saved atom ordering before optimization.

    A product reached by another H-transfer path may have the same configured
    identity but different H indices. Rebuild an unoptimized product using the
    trusted input order; subsequent normal reads then load coordinates and all
    indexed tensors in that order. Discovery endpoint snapshots stay untouched.
    """
    if not hasattr(qc, 'db') or getattr(species, 'wellorts', 0):
        return species
    from kinbot.stereo_routing import _row_species
    raw = _raw(qc.db)
    job = resolve_job(qc.db, routing_name(species) + '_well')
    references = list(raw.select(name='stereochemistry/' + job))
    rows = list(raw.select(name=job))
    reference = references[-1] if references else None
    row = reference or (rows[-1] if rows and rows[-1].data.get('status') == 'normal' else None)
    if row is None:
        return species
    saved = _row_species(species, row, job, reference)
    if not same_species(species, saved, 'cached product'):
        return species
    import numpy as np
    if (list(species.atom) == list(saved.atom)
            and np.array_equal(species.bond01, saved.bond01)):
        return species
    # Rebuild rather than reorder an existing Hessian, rotor, or conformer
    # object. This helper is used only for initial isolated product objects.
    saved.characterize(bond_mx=saved.bond)
    saved.name = species.name
    for attribute in ('optical_reference', 'optical_population', 'optical_counting_scope'):
        if hasattr(species, attribute):
            setattr(saved, attribute, copy.deepcopy(getattr(species, attribute)))
    logging.getLogger('KinBot').warning(
        '%s: reusing the verified product in its saved atom order; the original '
        'IRC endpoint is retained separately.', job)
    return saved


def require_cache_atom_order(requested, observed, context):
    """Identity is permutation invariant; indexed cache tensors are not remapped."""
    if routing_name(requested) == str(requested.chemid):
        return
    import numpy as np
    elements_match = list(requested.atom) == list(observed.atom)
    graph_match = (requested.chemid != observed.chemid or
                   np.array_equal(requested.bond01, observed.bond01))
    if not elements_match or not graph_match:
        from kinbot.stereo_routing import refuse_routing
        refuse_routing(f'{context}: cached atom indexing differs; restart with the '
                       'original atom ordering or use a fresh calculation directory',
                       [requested, observed])


def _raw(db):
    return db.raw if isinstance(db, RoutingDatabase) else db


def _has_job(db, name):
    db = _raw(db)
    if next(db.select(name=name), None) is not None:
        return True
    root = Path(db.filename).parent
    return any((root / (name + suffix)).exists() for suffix in
               ('.py', '.pkl', '.log', '.out', '_sella.log'))


def resolve_job(db, job):
    """Prefer a current job; otherwise read its explicitly verified old name."""
    db = _raw(db)
    job = str(job)
    match = _JOB.search(job)
    if match is None or _has_job(db, job):
        return job
    rows = list(db.select(name='stereo_route/' + match.group(2)))
    if not rows:
        return job
    route = rows[-1].data
    legacy = job[:match.start(2)] + route['legacy'] + job[match.end(2):]
    if not route['prefix_alias'] and job not in route['jobs']:
        return job
    return legacy if _has_job(db, legacy) else job


class RoutingDatabase:
    """ASE database facade: reads honor verified aliases; writes stay literal."""
    def __init__(self, raw):
        self.raw = raw

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, 'raw'), name)

    def select(self, *args, **kwargs):
        if 'name' in kwargs:
            kwargs['name'] = resolve_job(self.raw, kwargs['name'])
        return self.raw.select(*args, **kwargs)

    def get(self, *args, **kwargs):
        if 'name' in kwargs:
            kwargs['name'] = resolve_job(self.raw, kwargs['name'])
        return self.raw.get(*args, **kwargs)


def connect(*args, **kwargs):
    return RoutingDatabase(ase_connect(*args, **kwargs))


def input_species(filename, *, legacy_evidence=False):
    from kinbot.stationary_pt import StationaryPoint
    data = json.loads(Path(filename).read_text())
    species = StationaryPoint('saved input', data.get('charge', 0), data.get('mult', 1),
        structure=data.get('structure'), smiles=data.get('smiles') or None)
    species.characterize()
    apply_input_reference(species, data)
    if (legacy_evidence and routing_name(species) != str(species.chemid)
            and not data.get('structure') and data.get('smiles')):
        from rdkit import Chem
        mol = Chem.MolFromSmiles(data['smiles'])
        if mol is None or any(str(info.specified) == 'Unspecified'
                              for info in Chem.FindPotentialStereo(mol)):
            from kinbot.stereo_routing import refuse_routing
            refuse_routing(f'{filename}: saved SMILES leaves configuration unspecified; '
                           'supply the original coordinate input', [species])
    return species


def apply_input_reference(species, parameters):
    """Restore the declared population even if its selected member is a mirror."""
    population = parameters.get('optical_population', 'specified')
    reference = parameters.get('stereo_reference')
    if reference is not None:
        observed = canonical_identity(species)
        allowed = {observed.get('id')}
        if population == 'racemic':
            allowed.add(observed.get('mirror_id'))
        if (reference.get('status') != 'assigned' or observed.get('status') != 'assigned'
                or reference.get('id') not in allowed):
            from kinbot.stereo_routing import refuse_routing
            refuse_routing('Serialized stereo reference disagrees with the input geometry/population', [species])
        # Recompute all reference fields instead of trusting serialized hashes
        # or a supplied mirror_id independently of the actual input geometry.
        if reference['id'] != observed['id']:
            import numpy as np
            observed = canonical_identity(species, np.asarray(species.geom) * [-1., 1., 1.])
        species.optical_reference = observed
    species.optical_population = population


def prepare_qc_routing(qc, species):
    """Reuse verified legacy results, otherwise start the configured job."""
    from kinbot.stereo_routing import StereoRoutingError
    try:
        _prepare_qc_routing(qc, species)
    except StereoRoutingError as error:
        key = routing_name(species)
        # With no prior alias the configured namespace is already an isolated
        # place for fresh work. An ambiguous old result need not block it.
        if (key != str(species.chemid) and hasattr(qc, 'db')
                and not list(_raw(qc.db).select(name='stereo_route/' + key))):
            logging.getLogger('KinBot').warning(
                '%s: leaving unverified legacy calculations untouched and using '
                'the configured job name. %s', key, error)
            return
        raise


def _prepare_qc_routing(qc, species):
    """Bind a stereo namespace only when the old input/result establishes it.

    A completed result alone licenses that well job, not all old conformer/TS
    children. An original saved input licenses the prefix. Conflicting old
    identities are left separate; unidentified old jobs require user resolution.
    """
    if getattr(species, 'wellorts', 0):
        return
    from kinbot.stereo_identity import optical_scope
    species.optical_population = getattr(qc, 'par', {}).get('optical_population', 'specified')
    optical_scope(species, species.optical_population)
    key = routing_name(species)
    legacy = str(species.chemid)
    if key == legacy or not hasattr(qc, 'db'):
        return
    db = _raw(qc.db)
    from kinbot.stereo_routing import _row_species, guard_well_job, refuse_routing
    routes = list(db.select(name='stereo_route/' + key))
    route = routes[-1] if routes else None
    reference = list(db.select(name=f'stereochemistry/{legacy}_well'))
    rows = list(db.select(name=legacy + '_well'))
    original = None
    prefix_alias = False
    if reference:
        original = _row_species(species, reference[-1], legacy + '_well')
        # A guard added after the result is evidence for that result only;
        # it must not masquerade as an original input for all earlier children.
        prefix_alias = not rows or reference[-1].id < rows[0].id
    supplied = getattr(qc, 'par', {}).get('stereo_legacy_inputs', {}).get(legacy)
    filename = Path(supplied) if supplied else Path(db.filename).parent / f'{legacy}.json'
    if filename.exists():
        saved = input_species(filename, legacy_evidence=True)
        if original is not None and not same_species(original, saved):
            refuse_routing(f'{legacy}: saved input and database reference disagree', [original, saved])
        original, prefix_alias = saved, True
    if rows and rows[-1].data.get('status') == 'normal' and not prefix_alias:
        observed = _row_species(species, rows[-1], legacy + '_well')
        if any(not state['observations']
               for state in observed.calculation_state_evidence.values()):
            refuse_routing(f'{legacy}: completed result does not record charge and multiplicity; '
                           'set stereo_legacy_inputs to the original KinBot JSON input',
                           [species, observed])
    evidence = original
    if evidence is None and rows and rows[-1].data.get('status') == 'normal':
        evidence = _row_species(species, rows[-1], legacy + '_well')
    if evidence is None:
        if _has_job(db, legacy + '_well') or getattr(qc, 'job_ids', {}).get(legacy + '_well'):
            refuse_routing(f'{legacy}: old job has no identifiable original input. '
                'Set stereo_legacy_inputs to the original KinBot JSON input and retry', [species])
        return
    identity = canonical_identity(evidence)
    declared = getattr(evidence, 'optical_reference', None) or identity
    expected = key.split('-s', 1)[1]
    if identity.get('status') != 'assigned':
        refuse_routing(f'{legacy}: legacy configuration cannot be assigned', [species, evidence])
    if declared.get('id') != expected:
        if route is not None:
            refuse_routing(f'{legacy}: previously verified legacy input has changed configuration',
                           [species, evidence])
        if not any(any(tag in graph for tag in ('@', '/', '\\'))
                   for graph in identity['canonical_graphs']):
            refuse_routing(f'{legacy}: original input does not resolve the requested configuration',
                           [species, evidence])
        return  # It belongs to another configured species; use the new name.
    # Verify explicit result state and input/output agreement before aliasing.
    require_cache_atom_order(species, evidence, legacy)
    guard_well_job(qc, evidence, evidence.geom, legacy + '_well')
    if route is not None and (route.data['prefix_alias'] or not prefix_alias):
        return
    from ase import Atoms
    db.write(Atoms(species.atom, positions=species.geom), name='stereo_route/' + key,
        data={'schema': 'kinbot.routing.v1', 'canonical_key': key, 'identity': identity,
              'legacy': legacy, 'prefix_alias': prefix_alias,
              'jobs': [key + '_well', 'stereochemistry/' + key + '_well']})


def routed_qc_job(qc, species, job):
    prepare_qc_routing(qc, species)
    source = resolve_job(qc.db, job) if hasattr(qc, 'db') else job
    if source != job and qc.check_qc(source) == 0:
        return job  # Missing/invalidated old calculation needs a fresh job.
    return source


def prepare_pes_directory(root, species):
    """Reuse a verified old PES directory through an explicit directory alias.

    Original inputs/results are retained. New output uses configured names.
    An existing directory without its original input is never guessed to match.
    """
    root = Path(root)
    key = routing_name(species)
    target = root / key
    legacy = root / str(species.chemid)
    if key != legacy.name and target.is_symlink():
        from kinbot.stereo_routing import refuse_routing
        filename = target / f'{species.chemid}.json'
        if (not filename.exists() or
                not same_species(species, input_species(filename, legacy_evidence=True), 'PES directory alias')):
            refuse_routing(f'{target}: directory alias does not identify the requested configuration',
                           [species])
    if key != legacy.name and not target.exists() and legacy.exists():
        filename = legacy / f'{species.chemid}.json'
        from kinbot.stereo_routing import refuse_routing
        if not filename.exists():
            logging.getLogger('KinBot').warning(
                '%s: original input is unavailable; leaving the old PES directory '
                'untouched and starting the configured directory %s.', legacy, target)
        else:
            from kinbot.stereo_routing import StereoRoutingError
            try:
                previous = input_species(filename, legacy_evidence=True)
            except StereoRoutingError as error:
                logging.getLogger('KinBot').warning(
                    '%s: using fresh configured directory %s because %s', legacy, target, error)
            else:
                if same_species(species, previous, 'legacy PES directory'):
                    target.symlink_to(legacy.name, target_is_directory=True)
    target.mkdir(exist_ok=True, parents=True)
    return target


def expand_pes_names(root, choices):
    """Resolve connectivity selectors to the configured directories on disk."""
    root = Path(root)
    names = []
    for choice in choices:
        candidates = sorted(path for path in root.iterdir()
                            if path.is_dir() and is_species_name(path.name)
                            and matches_name(path.name, [str(choice)]))
        if not candidates:
            names.append(str(choice))  # Preserve the existing missing-well error.
        for path in candidates:
            filename = path / f'{path.name}.json'
            if '-s' not in path.name and filename.exists():
                species = input_species(filename, legacy_evidence=True)
                path = prepare_pes_directory(root, species)
            if path.name not in names:
                names.append(path.name)
    return names


def configured_result_matches(db, name, row, population=None):
    """Validate a configured PES energy record, including its trusted labels."""
    match = _KEY.fullmatch(str(name))
    if match is None or match.group(2) is None:
        return True
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.stereo_routing import _row_species
    source = resolve_job(db, f'{name}_well')
    references = list(db.select(name=f'stereochemistry/{source}'))
    reference = references[-1] if references else None
    template = StationaryPoint('energy input',
        reference.data['input_charge'] if reference else 0,
        int(match.group(1)[-1]), atom=row.symbols, geom=row.positions)
    observed = _row_species(template, row, str(name), reference)
    identity = canonical_identity(observed)
    allowed = {match.group(2)}
    if reference is not None:
        scope = reference.data.get('chemical_context', {})
        requested = reference.data.get('identity', {})
        declared = scope.get('optical_reference') or requested
        if ((population or scope.get('optical_population')) == 'racemic'
                and declared.get('id') == match.group(2)
                and requested.get('id') in {declared.get('id'), declared.get('mirror_id')}):
            allowed.add(declared.get('mirror_id'))
    if population == 'racemic' and identity.get('mirror_id') == match.group(2):
        allowed.add(identity['id'])
    return identity.get('status') == 'assigned' and identity['id'] in allowed
