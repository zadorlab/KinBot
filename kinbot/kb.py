from kinbot.species_routing import routing_name, apply_input_reference, matches_name
from contextlib import contextmanager
import fcntl
import os
import sys
import datetime
import copy
import logging
from pathlib import Path
from queue import Queue
import threading
import time
import traceback

from kinbot import license_message
from kinbot import postprocess
from kinbot.parameters import Parameters
from kinbot.mess import MESS
from kinbot.optimize import Optimize
from kinbot.reaction_finder import ReactionFinder
from kinbot.reaction_finder_bimol import ReactionFinderBimol
from kinbot.reaction_generator import ReactionGenerator
from kinbot.stationary_pt import StationaryPoint
from kinbot.vrc_tst_scan import VTS
from kinbot.electronic_states import apply_electronic_state_model
from kinbot.qc import QuantumChemistry
from kinbot.utils import make_dirs, clean_files
from kinbot.config_log import config_log
from kinbot.execution import atomic_json as _atomic_json
from kinbot.run_format import ensure_current_run
from kinbot.stereo_identity import require_supported_identity, log_input_stereochemistry, UnsupportedStereochemistry
from kinbot.stereo_routing import StereoRoutingError
from kinbot.theory import automatic_cbh_anl_requested


logger = logging.getLogger('KinBot')


def _configure_automatic_queue_admission(parameters, requested):
    """Scope and serialize scheduler admission for one automatic run.

    The ordinary reaction search and both post-reaction branches belong to
    one calculation tree.  Establish the shared scope before constructing
    any QC backend, including both members of a profiled L1/L2 router.
    Generated reference controllers receive their narrower scope from the
    automatic dispatcher and do not enter this function as requested runs.
    """
    if (not requested or parameters.get('queuing') == 'local'
            or int(parameters.get('queue_job_limit', -1)) <= 0):
        return None
    root = Path.cwd().resolve()
    if os.environ.get('KINBOT_PES_PHASE_ROOT') == str(root.parent):
        root = root.parent
    os.environ['KINBOT_QUEUE_ROOT'] = str(root)
    os.environ['KINBOT_QUEUE_ROOT_RECURSIVE'] = '1'
    os.environ['KINBOT_QUEUE_LOCK'] = str(root / '.kinbot_submission.lock')
    return root


@contextmanager
def _driver_lock():
    """Allow only one KinBot controller in a calculation directory."""
    ensure_current_run(create=True)
    from kinbot.run_format import MARKER
    with open(MARKER, 'r', encoding='utf-8') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(
                'Another KinBot driver is already managing this calculation '
                'directory.') from error
        yield


def main():
    with _driver_lock():
        try:
            from kinbot.driver_status import run_recorded
            return run_recorded(_main, 'kinbot_status.json')
        except (StereoRoutingError, UnsupportedStereochemistry) as error:
            logging.getLogger('KinBot').error('KinBot stopped: %s', error)
            raise SystemExit(1) from None


def _accepted_homolytic_channels(species):
    """Accepted homolytic or product-complex dissociation paths for VRC."""
    statuses = getattr(species, 'reac_ts_done', ())
    families = getattr(species, 'reac_type', ())
    selected = []
    for index, reaction in enumerate(getattr(species, 'reac_obj', ())):
        if (index < len(statuses) and statuses[index] == -1
                and ((index < len(families) and families[index] == 'hom_sci')
                     or getattr(reaction, 'do_vdW', False))
                and len(reaction.products) == 2):
            selected.append(reaction.instance_name)
    return selected


def _check_initial_minimum_frequencies(frequencies, species_name):
    """Apply the stable-well check without assuming two vibrational modes."""
    if frequencies and frequencies[0] <= 0:
        logger.warning('First frequency is %s for initial structure %s.',
                       frequencies[0], species_name)
        frequencies[0] *= -1.
    if len(frequencies) > 1 and frequencies[1] <= 0:
        logger.error('Second frequency is %s for initial structure %s.',
                     frequencies[1], species_name)
        return False
    return True


def _execute_vrc_branch(well, parameters, qc):
    logger.info('Starting the independent VRC-TST branch')
    previous_limit = qc.queue_job_limit
    previous_root = getattr(qc, 'queue_job_root', None)
    previous_recursive = getattr(qc, 'queue_job_root_recursive', False)
    qc.queue_job_limit = parameters.get(
        'vrc_tst_branch_slots', previous_limit)
    if 'vrc_tst_branch_slots' in parameters:
        qc.queue_job_root = str(Path.cwd().resolve())
        qc.queue_job_root_recursive = False
    try:
        # VRC must not overwrite the L2 geometries retained for live revisions.
        from kinbot.vrc_orbitals import copy_fragment_geometries
        vrc_well = copy.copy(well)
        vrc_well.reac_obj = []
        for reaction in well.reac_obj:
            clone = copy.copy(reaction)
            clone.products = copy_fragment_geometries(reaction.products)
            vrc_well.reac_obj.append(clone)
        vts = VTS(vrc_well, parameters, qc)
        from kinbot.vrc_multistate import SurfaceModelRequired
        from kinbot import live_control, pes_phase
        review = Path('vrc_model_review.json')
        try:
            vts.calculate_correction_potentials()
        except SurfaceModelRequired as error:
            _atomic_json(review, dict(status='needs_review', error=str(error),
                                     channels=getattr(vts, 'review_requests', [])))
            if not live_control.find() or not pes_phase.active_worker(Path.cwd()):
                raise
            logger.warning('%s Waiting for a live reviewed-model update.', error)
            while True:
                pes_phase.checkpoint()
                time.sleep(1)
        else:
            if review.exists():
                _atomic_json(review, dict(status='resolved'))
    finally:
        qc.queue_job_limit = previous_limit
        qc.queue_job_root = previous_root
        qc.queue_job_root_recursive = previous_recursive
    logger.info('Independent VRC-TST branch is complete')
    return {'status': 'complete'}


def _execute_automatic_branch(well, parameters, input_file, directory, plan):
    from kinbot.anl.automatic import (run_automatic_cbh_anl,
                                      validate_automatic_environment)
    logger.info('Starting the independent CBH/ANL composite branch')
    validate_automatic_environment(
        parameters, include_mess=False, include_vrc=False)
    result = run_automatic_cbh_anl(
        well, parameters, input_file, directory=directory, plan=plan)
    logger.info('Independent CBH/ANL composite branch is complete')
    return result


def _post_reaction_branches(well, parameters, qc, input_file,
                            automatic_requested):
    """Run VRC and composite work independently, then open the MESS gate.

    The automatic branch consumes a frozen JSON snapshot.  It never reads or
    mutates the live reaction objects used by VRC.  Their filesystem outputs
    are disjoint, and this controller is the sole writer of the join record.
    """
    vrc_requested = bool(parameters.get('vrc_tst_scan')
                         or parameters.get('vrc_tst_noscan'))
    if not automatic_requested and not vrc_requested:
        return None

    from kinbot import run_policy
    shared_pool = run_policy.find()
    total_slots = parameters.get('queue_job_limit', -1)
    if isinstance(total_slots, bool) or not isinstance(total_slots, int):
        raise ValueError('queue_job_limit must be an integer.')
    total_slots = max(1, total_slots)
    concurrent = automatic_requested and vrc_requested and total_slots > 1
    automatic_parameters = copy.deepcopy(parameters)
    vrc_parameters = copy.deepcopy(parameters)
    if concurrent:
        # Preserve the frozen plan/resource contract on restart. A live
        # shared policy overrides these legacy caps at admission time; its
        # cross-process lock replaces the static reservation in that mode.
        vrc_slots = max(1, total_slots // 2)
        automatic_slots = total_slots - vrc_slots
        automatic_parameters['queue_job_limit'] = automatic_slots
        vrc_parameters['vrc_tst_max_nodes'] = min(
            int(vrc_parameters['vrc_tst_max_nodes']), vrc_slots)
        vrc_parameters['rotdpy_max_jobs'] = min(
            int(vrc_parameters['rotdpy_max_jobs']), vrc_slots)
        vrc_parameters['vrc_tst_branch_slots'] = vrc_slots
        if not shared_pool:
            logger.info(
                'Running VRC-TST and CBH/ANL concurrently under the %d-slot cap: '
                '%d VRC slot(s), %d composite slot(s).',
                total_slots, vrc_slots, automatic_slots)
    else:
        vrc_slots = total_slots if vrc_requested else 0
        automatic_slots = total_slots if automatic_requested else 0

    if shared_pool:
        logger.info('Shared scheduling pool %s: live limit %d; no branch reservations.',
                    shared_pool.parent, run_policy.limit(total_slots))
    directory = parameters.get('anl_automatic_directory', 'anl_automatic')
    plan = None
    automatic_well = well
    preparation_error = None
    join_error = None
    try:
        from kinbot.anl.automatic import validate_mess_environment
        validate_mess_environment(parameters)
    except Exception as error:
        # Missing final solvers block only the join.  They cannot prevent the
        # independent VRC and composite branches from producing reusable
        # scientific results.
        join_error = error
    if automatic_requested:
        from kinbot.anl.automatic import prepare_automatic_cbh_anl
        try:
            from kinbot.pes_selection import selected_view
            automatic_well = selected_view(well, parameters)
            plan = prepare_automatic_cbh_anl(
                automatic_well, automatic_parameters, directory=directory)
        except Exception as error:  # preserve the independent VRC branch
            preparation_error = error

    state_path = Path('.kinbot_post_reaction_branches.json')
    state = {
        'schema': 1,
        'status': 'running',
        'mess_gate': 'waiting',
        'slot_cap': total_slots,
        'branches': {
            'vrc_tst': {
                'status': 'pending' if vrc_requested else 'not_requested',
                'slots': vrc_slots,
            },
            'cbh_anl': {
                'status': ('failed' if preparation_error else
                           ('pending' if automatic_requested else
                            'not_requested')),
                'slots': automatic_slots,
            },
        },
        'final_solver': {
            'status': 'failed' if join_error else 'ready',
        },
    }
    if shared_pool:
        state['resource_policy'] = str(shared_pool)
        state['slot_cap'] = 'live shared pool'
        for branch in state['branches'].values():
            branch['slots'] = 'shared'

    def failure_details(error):
        return {'error': str(error), 'error_type': type(error).__name__,
                'traceback': ''.join(traceback.format_exception(
                    type(error), error, error.__traceback__))}

    if preparation_error:
        state['branches']['cbh_anl'].update(failure_details(preparation_error))
    if join_error:
        state['final_solver'].update(failure_details(join_error))
    _atomic_json(state_path, state)

    results = {}
    errors = ([('cbh_anl', preparation_error)]
              if preparation_error is not None else [])
    requested_count = int(vrc_requested) + int(
        automatic_requested and preparation_error is None)
    state_lock = threading.Lock()

    def publish_branch(name, **update):
        # The composite thread can finish while VRC is still polling/running.
        # Publish that outcome immediately instead of hiding it until join().
        with state_lock:
            state['branches'][name].update(update)
            if update.get('status') == 'failed':
                state['mess_gate'] = 'blocked'
            _atomic_json(state_path, state)

    def record(name, function, *args):
        try:
            results[name] = function(*args)
        except Exception as error:
            errors.append((name, error))
            from kinbot.pes_phase import ExplorationReopened
            publish_branch(name, status=('paused' if isinstance(error, ExplorationReopened) else 'failed'),
                           **failure_details(error))
        else:
            publish_branch(name, status='complete')

    if concurrent and preparation_error is None:
        # Keep VRC in the main thread. Its ROTD_py wrapper installs signal
        # forwarding handlers so termination reaches its process group; Python
        # permits that only in the main thread. The immutable composite branch
        # is the sole background worker.
        state['branches']['vrc_tst']['status'] = 'running'
        state['branches']['cbh_anl']['status'] = 'running'
        _atomic_json(state_path, state)
        channel = Queue(maxsize=1)

        def run_composite():
            try:
                value = _execute_automatic_branch(
                    automatic_well, automatic_parameters, input_file, directory, plan)
            except BaseException as error:
                from kinbot.pes_phase import ExplorationReopened
                publish_branch('cbh_anl', status=('paused' if isinstance(error, ExplorationReopened) else 'failed'),
                               **failure_details(error))
                channel.put(('error', error))
            else:
                publish_branch('cbh_anl', status='complete')
                channel.put(('result', value))

        # A daemon worker is deliberate. On SIGTERM the main-thread VRC
        # wrapper terminates its ROTD_py process group and exits immediately;
        # Python must not then wait indefinitely for a composite polling
        # thread. Submitted composite jobs remain restartable from atomic
        # dispatcher state, and reference child PIDs are recorded separately.
        worker = threading.Thread(
            target=run_composite, name='kinbot-composite', daemon=True)
        worker.start()
        record('vrc_tst', _execute_vrc_branch,
               well, vrc_parameters, qc)
        worker.join()
        kind, payload = channel.get()
        if kind == 'error':
            errors.append(('cbh_anl', payload))
            from kinbot.pes_phase import ExplorationReopened
            state['branches']['cbh_anl'].update(
                status=('paused' if isinstance(payload, ExplorationReopened) else 'failed'), error=str(payload))
        else:
            results['cbh_anl'] = payload
            state['branches']['cbh_anl']['status'] = 'complete'
        _atomic_json(state_path, state)
    else:
        if vrc_requested:
            state['branches']['vrc_tst']['status'] = 'running'
            _atomic_json(state_path, state)
            record('vrc_tst', _execute_vrc_branch,
                   well, vrc_parameters, qc)
        if automatic_requested and preparation_error is None:
            state['branches']['cbh_anl']['status'] = 'running'
            _atomic_json(state_path, state)
            record('cbh_anl', _execute_automatic_branch,
                   automatic_well, automatic_parameters, input_file, directory, plan)

    from kinbot.pes_phase import ExplorationReopened, checkpoint
    if any(isinstance(error, ExplorationReopened) for _, error in errors):
        state.update(status='exploration_reopened', mess_gate='waiting')
        _atomic_json(state_path, state)
        raise ExplorationReopened('Postprocessing quiesced for network revision.')
    checkpoint()
    if join_error is not None:
        errors.append(('mess', join_error))
    if errors:
        state.update(status='failed', mess_gate='blocked')
        _atomic_json(state_path, state)
        names = ', '.join(name for name, _ in errors)
        raise RuntimeError(
            f'Post-reaction branch failure ({names}); MESS was not started.') \
            from errors[0][1]
    if requested_count:
        state.update(status='complete', mess_gate='ready')
        _atomic_json(state_path, state)
    return results.get('cbh_anl')


def _main():
    if sys.version_info.major < 3 or sys.version_info.minor < 8:
        print('KinBot only runs with python 3.8 or higher. You have python '
              f'{sys.version_info.major}.{sys.version_info.minor}. '
              'Please update python to run KinBot!')
        sys.exit(-1)

    try:
        input_file = sys.argv[1]
    except IndexError:
        print('To use KinBot, supply one argument being the input file!')
        sys.exit(-1)

    ensure_current_run(create=True)
    print(license_message.message)
    global logger
    logger = config_log('KinBot')
    # initialize the parameters for this run
    masterpar = Parameters(input_file, show_warnings=True)
    from kinbot.qc_scratch import configure as configure_qc_scratch
    configure_qc_scratch(masterpar.par)
    par = masterpar.par
    overrides = Path('.kinbot_search_options.json')
    if overrides.exists():
        import json
        from kinbot.network_revisit import settings
        from kinbot.vrc_multistate import normalize_inline_model_paths
        saved = json.loads(overrides.read_text())
        settings(saved.get('vdW_detection'))
        if set(saved) - {'vdW_detection', 'vrc_tst_surface_models'}:
            raise ValueError('Unknown persisted live search option.')
        normalize_inline_model_paths(saved, Path.cwd())
        par.update(saved)
    input_file = masterpar.input_file
    from kinbot.run_policy import initialize as initialize_resource_policy
    initialize_resource_policy(Path(input_file).resolve().parent, par['queue_live_policy'])
    from kinbot.live_control import initialize as initialize_live_control
    initialize_live_control(Path(input_file).resolve().parent, par['queue_live_control'])
    if (par.get('vrc_tst_surface_model') or par.get('vrc_tst_surface_model_file')):
        # Validate reviewed chemistry inputs before expensive discovery/QC.
        from kinbot.vrc_multistate import load_surface_model
        load_surface_model(par)
    automatic_run_requested = automatic_cbh_anl_requested(par)
    automatic_queue_root = _configure_automatic_queue_admission(
        par, automatic_run_requested)
    # set up the logging environment
    if par['verbose']:
        logger = config_log('KinBot', level='debug')

    # write the license message and the parameters to the log file
    logger.info('Input parameters')
    par_str = "\n\t".join([str(p) + ": " + str(par[p])for p in par])

    logger.info(par_str)
    # time stamp of the KinBot start
    logger.info('Starting KinBot')
    if automatic_queue_root is not None:
        from kinbot import run_policy
        resource_policy = run_policy.find()
        if resource_policy:
            logger.info('Shared scheduler admission: live limit %s under %s.',
                        run_policy.limit(par['queue_job_limit']), resource_policy.parent)
        else:
            logger.info(
                'Automatic scheduler admission: at most %s live jobs under %s.',
                par['queue_job_limit'], automatic_queue_root)

    make_dirs(par)
    from kinbot import pes_phase
    pes_phase.attach()
    rg = None
    auto_vrc = (par.get('rotdpy_run') and not par['vrc_tst_scan']
                and not par['vrc_tst_noscan'])

    if par['bimol'] == 0:
        # initialize the reactant
        well0 = StationaryPoint('well0',
                                par['charge'],
                                par['mult'],
                                smiles=par['smiles'],
                                structure=par['structure'])
        well0.short_name = 'w1'
        # write the initial reactant geometry to a file for visualization
        with open('initial_geometry.xyz', 'w') as geom_out:
            geom_out.write('{}\n\n'.format(well0.natom))
            for i, at in enumerate(well0.atom):
                x, y, z = well0.geom[i]
                geom_out.write('{} {:.6f} {:.6f} {:.6f}\n'.format(at, x, y, z))
            geom_out.write('\n\n')

        # characterize the initial reactant
        well0.characterize()
        apply_electronic_state_model(well0, par)
        require_supported_identity(well0)
        apply_input_reference(well0, par)
        log_input_stereochemistry(well0, par, logger)
        if par['cluster']:
            well0.make_hbonds()
            while 1:
                frags, maps = well0.start_multi_molecular(bond_mx=well0.bond)
                if len(frags) > 1:
                    p1, p2 = well0.make_extra_bond(frags[:2], maps[:2])
                    logger.info(f'Added extra bond between {p1} and {p2}')
                else:
                    break
#        well0.characterize()
        well0.name = routing_name(well0)
        start_name = well0.name
        if matches_name(well0.name, par['skip_chemids']):
            logger.info('This chemid is skipped, nothing to do here')
            logger.info('Finished KinBot at {}'.format(datetime.datetime.now()))
            print("Done!")
            return

        # initialize the qc instance
        qc = QuantumChemistry(par)
        if par['do_clean']:
            # Profiled theory can route L1 through a backend different from
            # the legacy top-level ``qc`` default.  Cleanup must follow the
            # backend that actually created the files; otherwise an L1
            # FairChem run can apply Gaussian rules to unrelated *.log files.
            clean_files(qc=qc.qc)

        # start the initial optimization of the reactant
        logger.info('Starting optimization of initial well...')
        qc.qc_opt(well0, well0.geom)
        err, well0.geom = qc.get_qc_geom(routing_name(well0) + '_well',
                                         well0.natom, wait=1)
        logger.debug(f'Initial well opt error {err}.')
        err, well0.freq = qc.get_qc_freq(routing_name(well0) + '_well',
                                         well0.natom, wait=1)
        logger.debug(f'Initial well freq error {err}, frequencies are {well0.freq}.')
        if err < 0:
            logger.error('Error with initial structure optimization.')
            return 1
        if not _check_initial_minimum_frequencies(well0.freq, well0.name):
            return 1

        # characterize again and look for differences
        initial_frequencies = list(well0.freq)
        starting_reference = getattr(well0, 'optical_reference', None)
        well0 = StationaryPoint('well0',
                                par['charge'],
                                par['mult'],
                                atom=copy.deepcopy(well0.atom),
                                geom=copy.deepcopy(well0.geom))
        # Rebuilding connectivity must not discard the checked QC frequencies.
        well0.freq = initial_frequencies
        well0.short_name = 'w1'
        well0.characterize()
        apply_input_reference(well0, dict(par, stereo_reference=starting_reference))
        if par['cluster']:
            well0.make_hbonds()
            while 1:
                frags, maps = well0.start_multi_molecular(bond_mx=well0.bond)
                if len(frags) > 1:
                    well0.make_extra_bond(frags[:2], maps[:2])
                else:
                    break
#        well0.characterize()
        well0.name = routing_name(well0)
        if well0.name != start_name and not par['cluster']:
            logger.error('The first well optimized to a structure different from the input.')
            return 1

        # do an MP2 optimization of the reactant,
        # to compare some scan barrier heigths to
        if qc.qc != 'nn_pes' and (par['families'] == ['all'] or \
                'birad_recombination_R' in par['families'] or \
                'r12_cycloaddition' in par['families'] or \
                'r14_birad_scission' in par['families'] or \
                'R_Addition_MultipleBond' in par['families'] or \
                (par['skip_families'] != ['none'] and
                ('birad_recombination_R' not in par['skip_families'] or
                'r12_cycloaddition' not in par['skip_families'] or
                'r14_birad_scission' not in par['skip_families'] or
                'R_Addition_MultipleBond' not in par['skip_families'])) or \
                par['reaction_search'] == 0):
            logger.debug('Starting MP2 optimization of initial well...')
            qc.qc_opt(well0, well0.geom, mp2=1)
            err, geom = qc.get_qc_geom(routing_name(well0) + '_well_mp2', well0.natom, 1)

        # comparison for barrierless scan
        if par['barrierless_saddle']:
            logger.debug('Optimization of initial well for barrierless at {}/{}'.
                    format(par['barrierless_saddle_method'], par['barrierless_saddle_basis']))
            qc.qc_opt(well0, well0.geom, bls=1)
            err, geom = qc.get_qc_geom(routing_name(well0) + '_well_bls', well0.natom, 1)
        
        # characterize again and look for differences
        well0.characterize()
        if par['cluster']:
            well0.make_hbonds()
            while 1:
                frags, maps = well0.start_multi_molecular(bond_mx=well0.bond)
                if len(frags) > 1:
                    well0.make_extra_bond(frags[:2], maps[:2])
                else:
                    break
#        well0.characterize()
        well0.name = routing_name(well0)

        err, well0.energy = qc.get_qc_energy(routing_name(well0) + '_well', 1)
        err, well0.zpe = qc.get_qc_zpe(routing_name(well0) + '_well', 1)
        # to save the starting energy to make thresholds consistent
        well0.start_energy = well0.energy
        well0.start_zpe = well0.zpe

        well_opt = Optimize(well0, par, qc, wait=1)
        well_opt.do_optimization()
        if well_opt.shigh == -999:
            logger.error('Error with high level optimization of initial structure.')
            return 1

        if par['reaction_search'] == 1:
            if par.get('profiled_l1_reuse_only'):
                logger.info('Reconstructing completed reaction network in '
                            'L1 reuse-only mode...')
            else:
                logger.info('Starting reaction search...')
            rf = ReactionFinder(well0, par, qc)
            rf.find_reactions()
            rg = ReactionGenerator(well0, par, qc, input_file)
            rg.generate()
            if (par.get('rotdpy_run') and not par['vrc_tst_scan']
                    and not par['vrc_tst_noscan']):
                selected_vrc = _accepted_homolytic_channels(well0)
                if selected_vrc:
                    par['vrc_tst_scan'] = {
                        routing_name(well0): selected_vrc}
                    logger.info(
                        'Automatically selected accepted homolytic channel(s) '
                        'for VRC-TST: %s', ', '.join(selected_vrc))

    # BIMOLECULAR REACTANTS
    elif par['bimol'] == 1:
        fragments = {'frag_a': None, 'frag_b': None} 
        # initialize the reacting fragments
        charge = 0
        for ii, frag in enumerate(fragments):
            fragments[frag] = StationaryPoint(frag,
                                              par['charge'][ii],
                                              par['mult'][ii],
                                              structure=par['structure'][ii])
            fragments[frag].characterize()
            apply_electronic_state_model(fragments[frag], par)
            require_supported_identity(fragments[frag])
            fragments[frag].name = routing_name(fragments[frag])
            charge += par['charge'][ii]
         
        # this is just formally a well
        well0 = StationaryPoint('bimolecular',
                                par['charge'],
                                par['mult'][2],
                                smiles=par['smiles'],
                                structure=par['structure'],
                                fragA=fragments['frag_a'],
                                fragB=fragments['frag_b'],
                                )
        well0.characterize()
        well0.short_name = 'w1'
        well0.chemid = '_'.join(sorted(str(frag.chemid) for frag in fragments.values()))
        well0.name = routing_name(well0)
 
        qc = QuantumChemistry(par)

        logger.info('\tStarting optimization of fragments...')
        for frag in fragments.values():
            qc.qc_opt(frag, frag.geom)
            err, frag.geom = qc.get_qc_geom(routing_name(frag) + '_well',
                                            frag.natom, wait=1)
            if err < 0:
                logger.error(f'Error with initial structure optimization of {frag.name}.')
                return 1
            err, frag.freq = qc.get_qc_freq(routing_name(frag) + '_well',
                                            frag.natom, wait=1)
            if frag.freq and frag.freq[0] <= 0. and frag.freq[0] >= -20.:
                logger.warning(f'Found imaginary frequency {frag.freq[0]} for {frag.name}. It is flipped.')
                frag.freq[0] *= -1.
            elif frag.freq and frag.freq[0] < -20.:
                logger.error(f'Found imaginary frequency {frag.freq[0]} for {frag.name}.')
                return 1
            if len(frag.freq) > 1 and frag.freq[1] <= 0:
                logger.error(f'Found two imaginary frequencies {frag.freq[1]} for {frag.name}.')
                return 1
            # characterize again and look for differences
            frag.characterize()
            if frag.name != routing_name(frag):
                logger.error(f'Reactant {frag} optimized to a structure different from the input.')
                return 1

        well0.energy = 0.
        well0.zpe = 0.
        well0.fragA = fragments['frag_a']  # update
        well0.fragB = fragments['frag_b']  # update
        well0.structure = [well0.fragA.structure, well0.fragB.structure]  # update, used in reactions
        for frag in fragments.values():
            err, frag.energy = qc.get_qc_energy(routing_name(frag) + '_well', 1)
            err, frag.zpe = qc.get_qc_zpe(routing_name(frag) + '_well', 1)
            well0.energy += frag.energy
            well0.zpe += frag.zpe

            frag_opt = Optimize(frag, par, qc, wait=1)
            frag_opt.do_optimization()
            if frag_opt.shigh == -999:
                logger.error(f'Error with high level optimization for {frag.name}.')
                return 1

        well0.fragA = fragments['frag_a']  # update
        well0.fragB = fragments['frag_b']  # update
        well0.energy = 0.
        well0.zpe = 0.
        for frag in fragments.values():
            well0.energy += frag.energy
            well0.zpe += frag.zpe

        if par['reaction_search'] == 1:
            logger.info('\tStarting bimolecular reaction search...')
            rf = ReactionFinderBimol(well0, par, qc)
            rf.find_reactions()
            rg = ReactionGenerator(well0, par, qc, input_file)
            rg.generate()

    while True:
        try:
            pes_phase.wait_for_exploration(rg)
            if auto_vrc:
                par['vrc_tst_scan'] = {routing_name(well0): _accepted_homolytic_channels(well0)}
            automatic = _post_reaction_branches(
                well0, par, qc, input_file, automatic_run_requested)
            from kinbot.pes_selection import selected_view
            downstream_well = selected_view(well0, par) if automatic_run_requested else well0
            if automatic is not None:
                par['anl_handoff_manifest'] = automatic['handoff_manifest']
                par['cbh_thermochemistry_models'] = automatic[
                    'thermochemistry_models']
                logger.info('Automatic CBH/ANL network is complete')

            if par.get('thermochemistry_model_output'):
                model_writer = MESS(par, well0)
                output = model_writer.write_thermochemistry_model(
                    well0, par['thermochemistry_model_output'])
                logger.info('Wrote accepted L2 thermochemistry model to %s', output)

            if par.get('anl_handoff_manifest'):
                from kinbot.anl.handoff import apply_handoff
                logger.info('Applying completed ANL/CBH results to the KinBot network')
                apply_handoff(downstream_well, par['anl_handoff_manifest'])

            if par['me'] > 0:  # it will be 2 for kinbots when the mess file is needed but not run
                mess = MESS(par, downstream_well)
                mess.write_input(qc)
                if par.get('run_messpf'):
                    logger.info('Starting MESS partition-function calculations')
                    mess.run_partition_functions()

            # Preserve the completed search even if the separate rate calculation fails.
            postprocess.create_summary_file(downstream_well, qc, par)
            postprocess.createPESViewerInput(downstream_well, qc, par)
            postprocess.creatMLInput(downstream_well, qc, par)

            if par['me'] == 1 and par['me_code'] == 'mess':
                logger.info('Starting Master Equation calculations')
                mess.run()

            if automatic is not None:
                from kinbot.anl.automatic import finalize_automatic_run
                report = finalize_automatic_run(automatic, par)
                logger.info('Wrote automatic CBH/ANL final report to %s', report)

            pes_phase.assembled()
            break
        except pes_phase.ExplorationReopened:
            if automatic_run_requested:
                from kinbot.anl.handoff import detach_handoff
                detach_handoff(well0)
                par.pop('anl_handoff_manifest', None)
                par.pop('cbh_thermochemistry_models', None)
            logger.info('Reopening PES exploration; existing QC jobs remain owned and reusable.')

    clean_files(diagnostic=True, qc=qc.qc)
    logger.info('KinBot finished.')
    try:
        print("Done!")
    except OSError:
        pass
    return


if __name__ == "__main__":
    raise SystemExit(main())
