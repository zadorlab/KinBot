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
    'vrc_tst_high_extrapolation_power',
    'vrc_tst_state_average_roots',
}


def _positive_integer(value, name):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f'{name} must be a positive integer.')
    return value


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


def load_surface_model(parameters):
    name = parameters.get('vrc_tst_surface_model_file', '')
    if not name:
        return None
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
    return model, {'source': str(path), 'sha256': file_sha256(path)}


def channel_for_reaction(parameters, reaction):
    loaded = load_surface_model(parameters)
    if loaded is None:
        return None
    model, provenance = loaded
    bond = sorted(int(index) + 1 for index in reaction.scan_coo)
    matches = [channel for channel in model['channels'] if sorted(channel['bond']) == bond]
    if len(matches) != 1:
        raise ValueError(f'Accepted VRC bond {bond} has no unique reviewed surface channel.')
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
    from kinbot.vrc_correction import coupled_multiplicities, correction_model
    par = deepcopy(parameters)
    par.update(deepcopy(stage['parameters']))
    par['vrc_tst_surface_model_file'] = ''  # scalar leaf, no recursive expansion
    par['mult'] = state['multiplicity']
    spin = deepcopy(par.get('vrc_tst_spin_splitting') or {})
    spin.pop('singlet_multiplicity', None)
    spin.pop('triplet_multiplicity', None)
    spin['target_multiplicity'] = state['multiplicity']
    spin['fragment_multiplicities'] = list(fragment_multiplicities)
    par['vrc_tst_spin_splitting'] = spin
    # The highest-spin leaf needs no low-spin-gap construction. Its trusted
    # direct MR level remains explicit and state-specific.
    if (correction_model(par) == 'spin_splitting'
            and state['multiplicity'] == max(coupled_multiplicities(*fragment_multiplicities))):
        par['vrc_tst_correction_model'] = 'direct'
    selected = {str(state['multiplicity']): state['root']}
    if correction_model(par) == 'spin_splitting':
        reference = spin.get('reference_multiplicity',
            max(coupled_multiplicities(*fragment_multiplicities)))
        selected[str(reference)] = 1
    par['vrc_tst_selected_roots'] = selected
    par['rotdpy_dist'] = list(stage['distances_angstrom'])
    par['_vrc_dividing_coordinate'] = stage['coordinate']
    return par


def calculate_channel(vts, name, selection):
    """Compute separate MR paths/corrections for all reviewed state leaves."""
    from kinbot.vrc_tst_scan import VTS
    from kinbot.molpro import vrc_multireference_energy_contract
    channel, provenance = selection
    reaction = vts.scan_reac[name]
    stages = []
    for stage in channel['bottlenecks']:
        leaves = []
        for state in stage['states']:
            par = state_parameters(vts.par, stage, state,
                                   [int(product.mult) for product in reaction.products])
            if vrc_multireference_energy_contract(par['vrc_tst_scan_method']) is None:
                raise ValueError('State-resolved VRC requires an explicit Molpro MR scan method.')
            leaf_name = f'{name}__{stage["name"]}_m{state["multiplicity"]}_r{state["root"]}'
            clone = copy(reaction)
            clone.species = copy(reaction.species)
            clone.species.mult = state['multiplicity']
            clone.instance_name = leaf_name
            child = VTS(vts.well, par, vts.qc)
            child.scan_reac[leaf_name] = clone
            child.do_scan([leaf_name])
            child.energies([leaf_name])
            correction = Path('vrctst') / f'corr_{leaf_name}.json'
            payload = json.loads(correction.read_text())
            # Persist only the resolved scalar controls needed by generation,
            # avoiding installation-specific transient parser objects.
            payload['surface_parameters'] = {key: deepcopy(value) for key, value in par.items()
                if key in _OVERRIDES or key in {'vrc_tst_selected_roots', 'mult',
                                               'rotdpy_dist', '_vrc_dividing_coordinate'}}
            payload['electronic_surface'] = deepcopy(state)
            atomic_json(correction, payload)
            leaves.append({**state, 'reaction': leaf_name,
                           'correction': str(correction.resolve()),
                           'correction_sha256': file_sha256(correction)})
        stages.append({**deepcopy(stage), 'states': leaves})
    master = {'schema': 1, 'kind': 'multistate', 'reaction': name,
              'surface_model': provenance, 'combination': channel['combination'],
              'equilibration': channel.get('equilibration'), 'bottlenecks': stages}
    path = Path('vrctst') / f'corr_{name}.json'
    atomic_json(path, master)
    return path


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
            'spin_counting': 'explicit_levels_or_multiplicity_in_flux; MESS_electronic_degeneracy=1'}


def stage_channel(manifest_file, destination, *, required_profile=None):
    """Revalidate every source and derive one MESS table per bottleneck."""
    from kinbot.rotdpy import number_of_states_file
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
        for state in stage['states']:
            inp = path.parent / state['input']
            if (inp.parent.resolve() != path.parent or file_sha256(inp) != state['input_sha256']
                    or file_sha256(state['correction']) != state['correction_sha256']):
                raise ValueError('State input or correction changed after multistate assembly.')
            table, provenance = number_of_states_file(inp, -1, required_profile=required_profile)
            entries.append({**state, 'source': str(table), 'provenance': provenance})
        destination = Path(destination)
        target = destination / f'rotd_{payload["reaction"]}_{stage["name"]}.out'
        if not _SAFE.fullmatch(target.name):
            raise ValueError('Unsafe multistate result name.')
        info = assemble_states(entries, target, baseline_cm_inverse=baseline)
        stages.append({'name': stage['name'], 'file': target.name, **info, 'states': entries})
    return {**payload, 'staged_bottlenecks': stages,
            'baseline_cm_inverse': baseline,
            'manifest_sha256': file_sha256(path)}


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
