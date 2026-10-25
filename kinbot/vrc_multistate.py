"""State-resolved VRC orchestration using the unmodified ROTD_py integrator.

Each spin/root has its own native path, correction, cache and scalar flux.
Only physical electronic contributions are summed. Correction variants are
alternative potentials, never extra electronic states. Conserved fragment
vibrations are common to the surfaces and are convolved by MESS afterwards.
"""
from copy import copy, deepcopy
import json
from pathlib import Path
import re

import numpy as np

from kinbot.execution import atomic_json, file_sha256


_SAFE = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')
_OVERRIDES = {
    'vrc_tst_active_space', 'vrc_tst_state_model', 'vrc_tst_state_model_file',
    'vrc_tst_scan_method', 'vrc_tst_scan_basis', 'vrc_tst_scan_gradient',
    'vrc_tst_scan_shift', 'vrc_tst_scan_ipea', 'vrc_tst_scan_points',
    'vrc_tst_sample_method', 'vrc_tst_sample_basis',
    'vrc_tst_sample_extrapolation_power', 'vrc_tst_high_method',
    'vrc_tst_high_basis', 'vrc_tst_correction_model', 'vrc_tst_spin_splitting',
    'vrc_tst_high_extrapolation_power', 'vrc_tst_reference_review',
    'vrc_tst_state_average_roots',
    'vrc_tst_selected_roots',
    'vrc_tst_asymptote_distance',
    'vrc_tst_fragment_geometry', 'vrc_tst_pivot_level',
    'vrc_tst_pivot_backend', 'vrc_tst_pivot_method', 'vrc_tst_pivot_basis',
}


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'{name} must be a positive integer.')
    return value


class SurfaceModelRequired(ValueError):
    """A discovered channel needs reviewed electronic-state input."""


def extend_surface_model(previous, current):
    """Permit added channels while preserving every existing MR contract."""
    current = validate_surface_model(current)
    if previous is not None:
        previous = validate_surface_model(previous)
        channels = {tuple(sorted(c['bond'])): c for c in current['channels']}
        if any(channels.get(tuple(sorted(c['bond']))) != c for c in previous['channels']):
            raise ValueError('Existing VRC channel settings changed; a scientific migration is required.')
    return current


def validate_surface_model(payload, *, reviewed=True):
    """Validate declarative channel/state choices without inventing chemistry."""
    model = deepcopy(payload)
    if (not isinstance(model, dict) or model.get('schema') != 1
            or set(model) - {'schema', 'reviewed', 'review_note', 'channels'}
            or not isinstance(model.get('channels'), list) or not model['channels']):
        raise ValueError('VRC surface model requires schema 1 and channels.')
    if reviewed and (model.get('reviewed') is not True
                     or not str(model.get('review_note', '')).strip()):
        raise ValueError('VRC surface model requires reviewed=true and a chemistry review_note.')
    bonds = set()
    for channel in model['channels']:
        if (not isinstance(channel, dict)
                or set(channel) - {'bond', 'combination', 'equilibration', 'bottlenecks'}):
            raise ValueError('Unknown VRC surface channel controls.')
        bond = channel.get('bond')
        if (not isinstance(bond, list) or len(bond) != 2 or bond[0] == bond[1]
                or any(isinstance(index, bool) or not isinstance(index, int)
                       or index < 1 for index in bond)):
            raise ValueError('VRC surface channel needs two one-based bond indices.')
        key = tuple(sorted(bond))
        if key in bonds:
            raise ValueError('Duplicate VRC surface channel bond.')
        bonds.add(key)
        combination = channel.get('combination', 'single')
        if combination not in ('single', 'statistical', 'minimum'):
            raise ValueError('VRC combination must be single, statistical, or minimum. '
                             'Complex wells use the existing KinBot network machinery.')
        channel['combination'] = combination
        stages = channel.get('bottlenecks')
        required_names = {'single'} if combination == 'single' else {'inner', 'outer'}
        if (not isinstance(stages, list) or len(stages) != len(required_names)
                or {stage.get('name') for stage in stages if isinstance(stage, dict)} != required_names):
            raise ValueError('VRC bottlenecks require single, or one inner and one outer stage.')
        if combination == 'statistical' and channel.get('equilibration') != 'rapid':
            raise ValueError('Statistical inner/outer VRC requires explicit rapid equilibration.')
        for stage in stages:
            if set(stage) - {'name', 'coordinate', 'distances_angstrom', 'states', 'parameters'}:
                raise ValueError('Unknown VRC bottleneck controls.')
            if stage.get('coordinate') not in ('reactive', 'com'):
                raise ValueError('VRC coordinate must be reactive or com.')
            distances = stage.get('distances_angstrom')
            if (not isinstance(distances, list) or len(distances) < 3
                    or any(isinstance(value, bool) or not isinstance(value, (int, float))
                           or not np.isfinite(value) or value <= 0. for value in distances)
                    or np.any(np.diff(distances) <= 0.)):
                raise ValueError('Each VRC bottleneck needs at least three increasing distances.')
            overrides = stage.get('parameters', {})
            if not isinstance(overrides, dict) or set(overrides) - _OVERRIDES:
                raise ValueError('Unsupported VRC bottleneck parameter override.')
            stage['parameters'] = overrides
            states = stage.get('states')
            if not isinstance(states, list) or not states:
                raise ValueError('Each VRC bottleneck needs electronic states.')
            seen = set()
            for state in states:
                if (not isinstance(state, dict)
                        or set(state) - {'multiplicity', 'root', 'weight',
                                         'energy_offset_cm_inverse', 'energy_reference_note', 'electronic_levels'}):
                    raise ValueError('Unknown VRC electronic state controls.')
                mult = _positive_integer(state.get('multiplicity'), 'State multiplicity')
                root = _positive_integer(state.get('root'), 'State root')
                if (mult, root) in seen:
                    raise ValueError('Duplicate electronic state in one VRC bottleneck.')
                seen.add((mult, root))
                weight = state.setdefault('weight', 1.)
                offset = state.get('energy_offset_cm_inverse')
                if (isinstance(weight, bool) or not isinstance(weight, (int, float))
                        or not np.isfinite(weight) or weight <= 0.
                        or isinstance(offset, bool) or not isinstance(offset, (int, float))
                        or not np.isfinite(offset)):
                    raise ValueError('VRC state needs positive weight and explicit finite energy offset.')
                if reviewed and not str(state.get('energy_reference_note', '')).strip():
                    raise ValueError('Every VRC state needs an energy_reference_note.')
                if 'electronic_levels' in state:
                    levels = state['electronic_levels']
                    if not isinstance(levels, list) or not levels:
                        raise ValueError('Surface electronic_levels must be nonempty.')
                    for level in levels:
                        if (not isinstance(level, dict)
                                or set(level) != {'energy_cm_inverse', 'degeneracy'}
                                or isinstance(level['energy_cm_inverse'], bool)
                                or not isinstance(level['energy_cm_inverse'], (int, float))
                                or not np.isfinite(level['energy_cm_inverse'])
                                or level['energy_cm_inverse'] < 0.):
                            raise ValueError('Invalid surface electronic level energy.')
                        _positive_integer(level['degeneracy'], 'Electronic level degeneracy')
    return model


def resolve_state_model_paths(payload, base):
    """Resolve binary record paths relative to the containing single input."""
    from kinbot.molpro import resolve_orbital_record_path
    model = deepcopy(payload)
    if isinstance(model, dict) and model.get('schema') == 2:
        if not isinstance(model.get('active_space'), dict) or not isinstance(model.get('basis_models'), dict) or not model['basis_models']:
            raise ValueError('Inline state model requires active_space and nonempty basis_models.')
        if any(not isinstance(entry, dict) for entry in model['basis_models'].values()):
            raise ValueError('Every inline basis model must be an object.')
        model['basis_models'] = {basis: resolve_orbital_record_path(entry, base)
                                 for basis, entry in model.get('basis_models', {}).items()}
    elif isinstance(model, dict):
        model = resolve_orbital_record_path(model, base)
    return model


def normalize_inline_model_paths(parameters, base):
    name = parameters.get('vrc_tst_state_model_file')
    if name:
        parameters['vrc_tst_state_model_file'] = str((Path(base)/Path(name).expanduser()).resolve())
    models = parameters.get('vrc_tst_surface_models', {})
    if not isinstance(models, dict):
        raise ValueError('vrc_tst_surface_models must map species selectors to reviewed surface models.')
    for selector, model in models.items():
        if not isinstance(selector, str) or not selector or not isinstance(model, dict):
            raise ValueError('Species surface models require nonempty selectors and inline model objects.')
        holder = {'vrc_tst_surface_model': model}
        normalize_inline_model_paths(holder, base)
        parameters['vrc_tst_surface_models'][selector] = holder['vrc_tst_surface_model']
    binding = parameters.get('vrc_tst_reference_review')
    if isinstance(binding, dict) and binding.get('path'):
        binding['path'] = str((Path(base)/Path(binding['path']).expanduser()).resolve())
    if parameters.get('vrc_tst_state_model') is not None:
        parameters['vrc_tst_state_model'] = resolve_state_model_paths(parameters['vrc_tst_state_model'], base)
    surface = parameters.get('vrc_tst_surface_model')
    if surface is not None:
        surface = validate_surface_model(surface, reviewed=False)
        for channel in surface['channels']:
            for stage in channel['bottlenecks']:
                normalize_inline_model_paths(stage['parameters'], base)
        parameters['vrc_tst_surface_model'] = surface


def load_surface_model(parameters):
    name = parameters.get('vrc_tst_surface_model_file', '')
    inline = parameters.get('vrc_tst_surface_model')
    if inline is not None and name:
        raise ValueError('Specify vrc_tst_surface_model or vrc_tst_surface_model_file, not both.')
    if inline is None and not name:
        return None
    base = (Path(parameters['input_file']).expanduser().resolve().parent
            if parameters.get('input_file') else Path.cwd())
    if inline is not None:
        import hashlib
        model = validate_surface_model(inline)
        holder = {'vrc_tst_surface_model': model}
        normalize_inline_model_paths(holder, base)
        model = holder['vrc_tst_surface_model']
        digest = hashlib.sha256(json.dumps(model, sort_keys=True).encode()).hexdigest()
        return model, {'source': 'inline:vrc_tst_surface_model', 'sha256': digest}
    path = Path(name).expanduser()
    if not path.is_absolute():
        base = (Path(parameters['input_file']).expanduser().resolve().parent
                if parameters.get('input_file') else Path.cwd())
        path = base / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f'Reviewed VRC surface model unavailable: {path}. '
                                'Complete orbital/state reconnaissance and review first.')
    model = validate_surface_model(json.loads(path.read_text()))
    # Paths inside an external model belong to that model, not the launch cwd.
    for channel in model['channels']:
        for stage in channel['bottlenecks']:
            name = stage['parameters'].get('vrc_tst_state_model_file')
            if name:
                stage['parameters']['vrc_tst_state_model_file'] = str(
                    (path.parent / name).resolve())
            normalize_inline_model_paths(stage['parameters'], path.parent)
    return model, {'source': str(path), 'sha256': file_sha256(path)}


def channel_for_reaction(parameters, reaction):
    # Atom/bond indices are local to a reviewed parent, not transferable to
    # every isomer reached by PES exploration.
    from kinbot.species_routing import routing_name
    models = parameters.get('vrc_tst_surface_models', {})
    scope = parameters.get('vrc_tst_surface_model_scope')
    require_parent_space = False
    if models or scope:
        parent = getattr(reaction, 'species', None)
        if parent is None:
            raise ValueError('Scoped VRC surface models require the reaction parent identity.')
        route = routing_name(parent)
        override = models.get(route, models.get(getattr(parent, 'smiles', '')))
        if override is not None:
            require_parent_space = bool(scope and scope != route)
            parameters = dict(parameters, vrc_tst_surface_model=override, vrc_tst_surface_model_file='')
        elif scope not in (None, '', route):
            raise SurfaceModelRequired(f'{route}: VRC needs a reviewed surface model for this newly explored well; '
                             'set vrc_tst_surface_models[species] in the input. '
                             'The initial well orbital indices cannot be silently reused.')
    loaded = load_surface_model(parameters)
    if loaded is None:
        return None
    model, provenance = loaded
    bond = sorted(int(index) + 1 for index in reaction.scan_coo)
    matches = [channel for channel in model['channels'] if sorted(channel['bond']) == bond]
    if not matches:
        raise SurfaceModelRequired(f'Accepted VRC bond {bond} needs a reviewed surface channel.')
    if len(matches) != 1:
        raise ValueError(f'Accepted VRC bond {bond} has no unique reviewed surface channel.')
    if require_parent_space:
        for stage in matches[0]['bottlenecks']:
            settings = stage['parameters']
            if not (settings.get('vrc_tst_state_model') or settings.get('vrc_tst_state_model_file')):
                raise SurfaceModelRequired(f'{route}, bond {bond}: each new-parent bottleneck needs its own '
                    'reviewed vrc_tst_state_model and basis-specific orbital selection; '
                    'an initial-parent active space will not be inherited.')
    from kinbot.vrc_correction import coupled_multiplicities
    from kinbot.electronic_states import state_average_roots
    allowed = coupled_multiplicities(*(int(product.mult) for product in reaction.products))
    roots = state_average_roots(reaction.products, allowed, parameters)
    for stage in matches[0]['bottlenecks']:
        from kinbot.vrc_correction import root_overrides
        counts = {**roots, **root_overrides(stage['parameters'], 'vrc_tst_state_average_roots')}
        for state in stage['states']:
            if state['multiplicity'] not in allowed or state['root'] > counts[state['multiplicity']]:
                raise ValueError('VRC surface state is outside the declared separated-fragment terms.')
    return deepcopy(matches[0]), provenance


def state_parameters(parameters, stage, state, fragment_multiplicities):
    """Resolve a scalar run while keeping the real accepted reaction intact."""
    from kinbot.vrc_correction import coupled_multiplicities, correction_model, root_overrides
    par = deepcopy(parameters)
    par.update(deepcopy(stage['parameters']))
    par['vrc_tst_surface_model_file'] = ''  # scalar leaf, no recursive expansion
    par['vrc_tst_surface_model'] = None
    par['mult'] = state['multiplicity']
    spin = deepcopy(par.get('vrc_tst_spin_splitting') or {})
    overrides = stage['parameters'].get('vrc_tst_spin_splitting')
    if isinstance(overrides, dict):
        # A channel may change only the auxiliary reference multiplicity;
        # retain the case's CC/gap basis and extrapolation policy.
        inherited = deepcopy(parameters.get('vrc_tst_spin_splitting') or {})
        overrides = deepcopy(overrides)
        for settings in (inherited, overrides):
            for old, new in (('singlet_multiplicity', 'target_multiplicity'),
                             ('triplet_multiplicity', 'reference_multiplicity'),
                             ('triplet_reference', 'reference')):
                if old in settings:
                    if new in settings:
                        raise ValueError(f'Specify {old} or {new}, not both.')
                    settings[new] = settings.pop(old)
        for key, value in overrides.items():
            if key in ('reference', 'weighting') and isinstance(value, dict):
                inherited[key] = {**inherited.get(key, {}), **deepcopy(value)}
            else:
                inherited[key] = deepcopy(value)
        spin = inherited
    spin.pop('singlet_multiplicity', None)
    spin['target_multiplicity'] = state['multiplicity']
    spin['fragment_multiplicities'] = list(fragment_multiplicities)
    par['vrc_tst_spin_splitting'] = spin
    # Default to direct for a highest coupled-spin leaf only when no separate
    # energy reference was requested. Explicit references need not be
    # reactive states of the ground-fragment channel.
    if (correction_model(par) == 'spin_splitting'
            and 'reference_multiplicity' not in spin
            and 'triplet_multiplicity' not in spin
            and state['multiplicity'] == max(coupled_multiplicities(*fragment_multiplicities))):
        par['vrc_tst_correction_model'] = 'direct'
    if correction_model(par) == 'direct':
        par.pop('vrc_tst_spin_splitting', None)
    selected = {str(key): value for key, value in root_overrides(
        par, 'vrc_tst_selected_roots').items()}
    selected[str(state['multiplicity'])] = state['root']
    # Do not overwrite the reviewed reference root. When absent, resolution
    # preserves the external state model's root or the automatic ground root.
    par['vrc_tst_selected_roots'] = selected
    par['rotdpy_dist'] = list(stage['distances_angstrom'])
    par['_vrc_dividing_coordinate'] = stage['coordinate']
    return par


def _surface_job(vts, name, stage, state):
    """Own coordinates and parameters; retain only read-only QC services."""
    from kinbot.vrc_tst_scan import VTS
    from kinbot.molpro import vrc_multireference_energy_contract
    from kinbot.vrc_orbitals import copy_fragment_geometries
    reaction = vts.scan_reac[name]
    par = state_parameters(vts.par, stage, state,
                           [int(product.mult) for product in reaction.products])
    if vrc_multireference_energy_contract(par['vrc_tst_scan_method']) is None:
        raise ValueError('State-resolved VRC requires an explicit Molpro MR scan method.')
    leaf_name = f'{name}__{stage["name"]}_m{state["multiplicity"]}_r{state["root"]}'
    clone = copy(reaction)
    clone.products = copy_fragment_geometries(reaction.products)
    clone.species = copy(reaction.species)
    clone.species.mult = state['multiplicity']
    clone.instance_name = leaf_name
    child = VTS(vts.well, par, vts.qc)
    child._surface_name = leaf_name
    child.scan_reac[leaf_name] = clone

    def execute():
        child.do_scan([leaf_name])
        child.energies([leaf_name])
        if getattr(child, '_scan_review_requests', []):
            from kinbot.vrc_scan_recovery import ScanReviewRequired
            raise ScanReviewRequired(child._scan_review_requests)
        correction = Path('vrctst') / f'corr_{leaf_name}.json'
        payload = json.loads(correction.read_text())
        payload['surface_parameters'] = {key: deepcopy(value) for key, value in par.items()
            if key in _OVERRIDES or key in {'vrc_tst_selected_roots', 'mult',
                                           'rotdpy_dist', '_vrc_dividing_coordinate'}}
        payload['electronic_surface'] = deepcopy(state)
        atomic_json(correction, payload)
        from kinbot.vrc_scan_recovery import bind_correction_review
        bind_correction_review(correction)
        return {**state, 'reaction': leaf_name,
                'correction': str(correction.resolve()),
                'correction_sha256': file_sha256(correction),
                '_plot': getattr(child, '_vrc_plot', None)}
    return leaf_name, execute


def calculate_channel(vts, name, selection):
    """Compatibility entry point for a single, serial channel."""
    return calculate_channels(vts, {name: selection}, limit=1)[name]


def calculate_channels(vts, selections, *, limit=None):
    """Advance independent reviewed leaves; keep each orbital chain serial.

    The same leaf names, input specifications and graph bindings are retained.
    Shared Slurm admission, not the number of controller threads, limits QC.
    Only the caller assembles channel records and subsequently launches ROTD.
    """
    from kinbot import run_policy
    from kinbot.vrc_parallel import run, SurfaceInterrupted
    from kinbot.pes_phase import ExplorationReopened
    from kinbot.vrc_scan_recovery import ScanReviewRequired
    import logging
    if limit is None:
        limit = vts.par.get('vrc_tst_parallel_surfaces', 4)
        if type(limit) is not int or limit < 1:
            raise ValueError('vrc_tst_parallel_surfaces must be a positive integer.')
        if run_policy.find() is None:
            limit = 1  # Independent graphs otherwise have independent node caps.
    logger = logging.getLogger('KinBot')
    jobs, stages, paths, outcomes = [], {}, {}, {}
    for name, (channel, provenance) in selections.items():
        stages[name] = []
        for stage in channel['bottlenecks']:
            keys = []
            for state in stage['states']:
                key, execute = _surface_job(vts, name, stage, state)
                if key in outcomes:
                    raise ValueError(f'Duplicate VRC surface identity: {key}')
                outcomes[key] = {'channel': name, 'status': 'queued'}
                jobs.append((key, execute))
                keys.append(key)
            stages[name].append((stage, keys))
    if not jobs:
        return paths
    status_path = Path('vrctst/parallel_surfaces.json')
    def publish():
        atomic_json(status_path, {'schema': 1, 'max_parallel_surfaces': limit,
                                 'surfaces': outcomes})
    publish()
    values, failures, reviews = {}, {}, []
    logger.info('Advancing %d reviewed VRC surfaces with up to %d local controllers; '
                'Slurm jobs use existing admission limits.', len(jobs), limit)

    def started(key):
        outcomes[key].update(status='running')
        publish()

    def completed(key, value, error):
        name = outcomes[key]['channel']
        if error is not None:
            paused = isinstance(error, (ExplorationReopened, SurfaceInterrupted))
            held = isinstance(error, ScanReviewRequired)
            outcomes[key].update(status='paused' if paused else 'needs_review' if held else 'failed', error=str(error))
            if held:
                reviews.extend(dict(c, channel=name) for c in error.channels)
            else:
                failures[key] = error
            logger.warning('VRC surface %s: %s', key, error)
        else:
            plot = value.pop('_plot', None)
            if plot is not None:
                from kinbot.utils import create_matplotlib_graph
                create_matplotlib_graph(**plot)
            values[key] = value
            outcomes[key].update(status='complete')
            from kinbot.rotdpy_campaign import ready_leaf
            try:
                ready_leaf(vts, name, value)
            except Exception as launch_error:
                outcomes[key]['sampling_error'] = str(launch_error)
                logger.exception('Cannot start ready ROTD surface %s', key)
            if all(k in values for _, keys in stages[name] for k in keys):
                channel, provenance = selections[name]
                assembled = [{**deepcopy(stage), 'states': [values[k] for k in keys]}
                             for stage, keys in stages[name]]
                master = {'schema': 1, 'kind': 'multistate', 'reaction': name,
                          'surface_model': provenance, 'combination': channel['combination'],
                          'equilibration': channel.get('equilibration'), 'bottlenecks': assembled}
                path = Path('vrctst') / f'corr_{name}.json'
                atomic_json(path, master)
                paths[name] = path
                from kinbot.rotdpy_campaign import ready
                try:
                    ready(vts, [name])
                except Exception as launch_error:
                    # Preserve all other surface work; final channel staging
                    # will report this error again rather than omit a channel.
                    outcomes[key]['sampling_error'] = str(launch_error)
                    logger.exception('Cannot start ready ROTD channel %s', name)
        publish()

    run(jobs, limit=limit, completed=completed, started=started)
    if failures:
        detail = ' | '.join(f'{key}: {error}' for key, error in failures.items())
        raise RuntimeError('VRC surface failure(s); independent work retained: ' + detail) from next(iter(failures.values()))
    if reviews:
        raise ScanReviewRequired(reviews)
    return paths


def create_inputs(par, payload, reactant, products, barrier, folder, correction_file):
    """Expand scalar ROTD runs sequentially, preserving the global job cap."""
    from kinbot.pes import create_rotdpy_inputs
    created = []
    stages = []
    for stage in payload['bottlenecks']:
        leaves = []
        for state in stage['states']:
            correction = Path(state['correction'])
            if file_sha256(correction) != state['correction_sha256']:
                raise ValueError('State-specific VRC correction changed after assembly.')
            leaf = json.loads(correction.read_text())
            leaf_par = deepcopy(par)
            leaf_par.update(leaf['surface_parameters'])
            leaf_par['vrc_tst_surface_model_file'] = ''
            leaf_par['vrc_tst_surface_model'] = None
            leaf_par['vrc_tst_scan'] = {reactant: [state['reaction']]}
            leaf_par['vrc_tst_noscan'] = {}
            files = create_rotdpy_inputs(leaf_par,
                [[reactant, state['reaction'], products, barrier]], [],
                correction_root=correction.parent.parent)
            if len(files) != 1:
                raise RuntimeError('State-resolved ROTD input was not generated.')
            path = Path(files[0]).resolve()
            created.extend(files)
            leaves.append({**state, 'input': path.name, 'input_sha256': file_sha256(path)})
        stages.append({**stage, 'states': leaves})
    manifest = {**payload, 'bottlenecks': stages,
                'correction': str(Path(correction_file).resolve()),
                'correction_sha256': file_sha256(correction_file)}
    atomic_json(Path(folder) / f'{payload["reaction"]}.multistate.json', manifest)
    return created


def _table(path):
    values = np.loadtxt(path, ndmin=2)
    if (values.shape[1] != 2 or len(values) < 2 or not np.isfinite(values).all()
            or np.any(values[:, 1] < 0.) or np.any(np.diff(values[:, 0]) <= 0.)):
        raise ValueError(f'Invalid state number-of-states table: {path}')
    # Same explicit zero-flux convention as the scalar MESS handoff.
    if values[0, 1] == 0.:
        return values
    if values[0, 0] < 0.:
        raise ValueError(f'{path}: negative-energy table needs a zero-flux lower boundary.')
    if values[0, 0] == 0. and values[0, 1] > 0.:
        values[0, 0] = min(.01, values[1, 0] / 2.)
    if values[0, 0] > 0.:
        values = np.vstack(([0., 0.], values))
    return values


def assemble_states(entries, destination, *, baseline_cm_inverse):
    """Sum spin-weighted physical states on one aligned microcanonical grid.

    Positive weights are additional spatial/path populations, not spin.
    Linear interpolation is bounded by the shared available upper energy;
    it never holds a final sample constant or extrapolates an unsampled tail.
    """
    # This is MESS Union algebra before the common conserved-mode convolution,
    # not an average of canonical rate coefficients. Temperature dependence
    # follows from integrating the aligned N(E); no extra 1/n or spin factor.
    tables = []
    for entry in entries:
        levels = entry.get('electronic_levels', [
            {'energy_cm_inverse': 0., 'degeneracy': entry['multiplicity']}])
        for level in levels:
            table = _table(entry['source'])
            shift = entry['energy_offset_cm_inverse'] + level['energy_cm_inverse'] - baseline_cm_inverse
            if shift < -1.e-9:
                raise ValueError('State energy offset lies below the common baseline.')
            table[:, 0] += max(0., shift)
            tables.append((table, level['degeneracy'] * entry.get('weight', 1.)))
    if not tables:
        raise ValueError('No electronic state fluxes supplied.')
    upper = min(table[-1, 0] for table, _ in tables)
    grid = np.unique(np.concatenate(([0.], *[table[:, 0] for table, _ in tables])))
    grid = grid[grid <= upper]
    if len(grid) < 3:
        raise ValueError('Insufficient common state energy range after alignment.')
    total = sum(weight * np.interp(grid, table[:, 0], table[:, 1], left=0.)
                for table, weight in tables)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(destination, np.column_stack((grid, total)), fmt='%.14e')
    from kinbot.rotdpy import number_of_states_ground
    return {'source_sha256': file_sha256(destination),
            'ground_cm_inverse': number_of_states_ground(destination),
            'upper_energy_cm_inverse': float(upper), 'baseline_cm_inverse': baseline_cm_inverse,
            'spin_counting': 'explicit_levels_or_multiplicity_in_flux; MESS_electronic_degeneracy=1',
            'electronic_assembly': 'sum_of_shifted_microcanonical_fluxes',
            'conserved_modes': 'common_fragment_modes_convolved_once_after_sum',
            'temperature_weighting': 'MESS integrates the combined N(E); no state-average weights'}


def channel_sources(manifest_file, *, required_profile=None):
    """Revalidate every physical-state source without writing MESS tables."""
    from kinbot.rotdpy import number_of_states_file
    from kinbot.vrc_scan_recovery import require_ready_for_rates
    path = Path(manifest_file).resolve()
    payload = json.loads(path.read_text())
    if payload.get('kind') != 'multistate' or payload.get('schema') != 1:
        raise ValueError('Invalid multistate VRC manifest.')
    source = payload['surface_model']
    if file_sha256(source['source']) != source['sha256']:
        raise ValueError('Reviewed VRC surface model changed after calculation.')
    if file_sha256(payload['correction']) != payload['correction_sha256']:
        raise ValueError('Multistate VRC correction changed after generation.')
    baseline = min(state['energy_offset_cm_inverse']
                   for stage in payload['bottlenecks'] for state in stage['states'])
    stages = []
    for stage in payload['bottlenecks']:
        entries = []
        states_seen, inputs_seen = set(), set()
        for state in stage['states']:
            require_ready_for_rates(state['correction'])
            inp = path.parent / state['input']
            key = (state['multiplicity'], state['root'])
            if key in states_seen or inp.resolve() in inputs_seen:
                raise ValueError('Duplicate electronic surface or reused scalar input in VRC flux manifest.')
            states_seen.add(key)
            inputs_seen.add(inp.resolve())
            if (inp.parent.resolve() != path.parent or file_sha256(inp) != state['input_sha256']
                    or file_sha256(state['correction']) != state['correction_sha256']):
                raise ValueError('State input or correction changed after multistate assembly.')
            table, provenance = number_of_states_file(inp, -1, required_profile=required_profile)
            entries.append({**state, 'source': str(table), 'provenance': provenance})
        stages.append({**stage, 'states': entries})
    return {**payload, 'bottlenecks': stages, 'baseline_cm_inverse': baseline,
            'manifest_sha256': file_sha256(path)}


def stage_channel(manifest_file, destination, *, required_profile=None):
    """Revalidate every source and derive one MESS table per bottleneck."""
    payload = channel_sources(manifest_file, required_profile=required_profile)
    baseline = payload['baseline_cm_inverse']
    stages = []
    for stage in payload['bottlenecks']:
        entries = stage['states']
        destination = Path(destination)
        target = destination / f'rotd_{payload["reaction"]}_{stage["name"]}.out'
        if not _SAFE.fullmatch(target.name):
            raise ValueError('Unsafe multistate result name.')
        info = assemble_states(entries, target, baseline_cm_inverse=baseline)
        stages.append({'name': stage['name'], 'file': target.name, **info, 'states': entries})
    return {**payload, 'staged_bottlenecks': stages,
            'baseline_cm_inverse': baseline}


def render_mess(provenance, render_rrho):
    """Use MESS's existing two-transition-state statistical/dynamical model."""
    stages = {stage['name']: stage for stage in provenance['staged_bottlenecks']}
    if provenance['combination'] == 'single':
        return render_rrho(stages['single'])
    method = {'statistical': 'statistical', 'minimum': 'dynamical'}[provenance['combination']]
    inner = render_rrho(stages['inner'])
    outer = render_rrho(stages['outer']).replace('    RRHO\n', '    OuterRRHO\n', 1)
    return ('    Variational\n' + f'      2TSMethod {method}\n'
            + outer + inner + '    End ! Variational\n')
