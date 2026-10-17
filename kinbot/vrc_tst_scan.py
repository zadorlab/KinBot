from kinbot.species_routing import configured_selection
from kinbot.species_routing import routing_key, routing_name, same_species, matches_name
from ase.db import connect
from operator import ne
from typing import Any
import numpy as np
import time
import logging
import copy
import os
import json
import shutil
import hashlib
from pathlib import Path
import re

from subprocess import run
from ase import Atoms
from ase.io import read
from ase.units import Hartree
from kinbot.utils import reorder_coord
from kinbot.stationary_pt import StationaryPoint
from kinbot import geometry
from kinbot.molpro import (Molpro, method_with_active_space,
                           vrc_multireference_energy_contract)
from kinbot.vrc_correction import (
    active_space_settings, asymptote_distance,
    assemble_sampling_reference, assemble_spin_splitting_reference,
    correction_components,
    correction_model, coupled_multiplicities, expected_levels,
    resolve_vrc_level, resolved_state_model, sampling_settings,
    spin_splitting_settings,
)
from kinbot.utils import create_matplotlib_graph, NpEncoder
from kinbot import constants
from kinbot.electronic_states import apply_electronic_state_model, state_average_roots

logger = logging.getLogger('KinBot')


def vrc_dispatch_resources(parameters):
    """One memory/rank contract shared by relaxed scans and correction jobs."""
    memory = parameters.get('vrc_tst_memory_mb')
    if memory is not None and (type(memory) is not int or memory < 1):
        raise ValueError('vrc_tst_memory_mb must be a positive integer or null.')
    return {'cores': 'auto', 'memory_mb': 'node' if memory is None else memory,
            'walltime': parameters['vrc_tst_walltime'],
            'max_cores': parameters['single_point_ppn'],
            'min_stack_mw': parameters['vrc_tst_min_stack_mw'],
            'partition': parameters['queue_name']}


def fragment_routing_state(species, fallback=None):
    """Serialize a configured fragment, with legacy connectivity fallback."""
    try:
        name = routing_name(species)
    except AttributeError:
        if fallback is None:
            raise
        name = str(fallback)
    return {'routing_key': name, 'charge': int(getattr(species, 'charge', 0)),
            'multiplicity': int(species.mult),
            'optical_population': getattr(species, 'optical_population', 'specified'),
            'stereo_reference': getattr(species, 'optical_reference', None)}


def _native_state_integer(value, label):
    """Return a JSON-native integer without truncating malformed state data."""
    if (isinstance(value, (bool, np.bool_))
            or not isinstance(value, (int, float, np.integer, np.floating))):
        raise ValueError(f'{label} must be an integer.')
    numeric = float(value)
    if not np.isfinite(numeric) or not numeric.is_integer():
        raise ValueError(f'{label} must be an integer.')
    return int(numeric)


def _drive_vrc_dispatch(run_dir, task_ids, *, interval=20., max_attempts=2):
    """Drive VRC tasks and retry one transient native/scheduler failure.

    A failed process launch, node loss, or scheduler disappearance is safe to
    resubmit because the dispatcher archives the complete prior task directory
    and binds the replacement to the same immutable input. Parser, provenance,
    constraint, and result-validation failures are never retried here.
    """
    from kinbot.anl.dispatch import (
        SchedulerUnavailableError, _load, drive_run,
        preflight_until_scheduler_ready, retry_failed)

    task_ids = tuple(task_ids)
    while True:
        returncode, summary = drive_run(run_dir, interval=interval)
        if not returncode:
            return returncode, summary
        _, _, state = _load(run_dir)
        retry = []
        for task_id in task_ids:
            entry = state['tasks'].get(task_id, {})
            if entry.get('status') != 'failed':
                continue
            execution_path = Path(run_dir) / 'tasks' / task_id / 'execution.json'
            try:
                execution = json.loads(execution_path.read_text())
            except (OSError, json.JSONDecodeError):
                execution = {}
            phase = entry.get('failure_phase') or execution.get(
                'failure_phase')
            attempt = int(entry.get('attempt', 1))
            if phase in ('execution', 'scheduler') and attempt < max_attempts:
                retry.append((task_id, phase, attempt))
        if not retry:
            return returncode, summary
        for task_id, phase, attempt in retry:
            while True:
                try:
                    archive = retry_failed(run_dir, task_id)
                    break
                except SchedulerUnavailableError as error:
                    logger.warning(
                        '%s: scheduler unavailable while preparing transient '
                        'retry; retaining failed attempt: %s', task_id, error)
                    time.sleep(float(interval))
            logger.warning(
                '%s: archived transient %s failure from attempt %d at %s; '
                'retrying once with the same input.',
                task_id, phase, attempt, archive)
        preflight_until_scheduler_ready(run_dir, interval=interval)


class VTS:
    """
    Class to run the scans for VRC-TST to
    eventually generate correction potentials
    """
    def __init__(self, well, par, qc):
        # instance of the VTS object
        self.well = well
        self.par = par
        self.qc = qc

        # reaction object for reactions scanned
        self.scan_reac = {}
        self.prod_chemids = []

    @staticmethod
    def _state_average_roots(products, multiplicities, parameters=None):
        """Infer C1 roots from the separated ground electronic terms."""
        return state_average_roots(products, multiplicities, parameters)

    def _spin_context(self, reaction):
        products = list(reaction.products)
        parent = reaction.species
        active_space = active_space_settings(
            self.par, getattr(parent, 'atom', None),
            getattr(parent, 'charge', 0))
        if (active_space is not None
                and active_space.get('selection')
                == 'two_doublet_radical_pair'):
            fragment_multiplicities = [int(product.mult)
                                       for product in products]
            if len(products) != 2 or fragment_multiplicities != [2, 2]:
                raise ValueError(
                    'vrc_tst_active_space=radical_pair is valid only for an '
                    'accepted channel containing two doublet fragments; use '
                    'full_valence or an explicit orbital model here.')
        settings = spin_splitting_settings(
            self.par, parent_multiplicity=int(reaction.species.mult),
            fragment_multiplicities=[int(product.mult)
                                     for product in products])
        roots = self._state_average_roots(
            products, (settings['target_multiplicity'],
                       settings['reference_multiplicity']), self.par)
        return settings, roots

    def calculate_correction_potentials(self):
        """
        The main driver for scanning the potential.
        """
        # Automatic homolytic-channel selection happens after KinBot's
        # initial ``make_dirs`` call.  VRC execution must therefore own its
        # directory precondition instead of relying on the input mapping to
        # have been populated during startup.
        for directory in ('vrctst/molpro', 'perm/vrctst/molpro'):
            Path(directory).mkdir(parents=True, exist_ok=True)
        self.gaussian_tools = {}
        if str(self.par.get('pp_orient', 'homo')).casefold() == 'homo':
            # HOMO-oriented pivot points require both native Gaussian
            # conversion utilities. Resolve them before submitting fragment
            # work so a missing artifact cannot surface later in ROTD_py.
            from kinbot.anl.site import _executable
            for program in ('formchk', 'cubegen'):
                executable = _executable(program, 'gaussian')
                if executable is None:
                    raise FileNotFoundError(
                        f'HOMO-oriented VRC-TST requires {program!r}; '
                        f'configure KINBOT_{program.upper()}_EXECUTABLE, '
                        'KINBOT_GAUSSIAN_ROOT, or PATH.')
                self.gaussian_tools[program] = str(executable)
        if self.par.get('rotdpy_run'):
            from kinbot.rotdpy import ensure_available
            ensure_available()
        # structure of par['vrc_tst_scan']:
        # {chemid1: ["reaction_name1","reaction_name2], chemid2: [...]}
        for option, noscan in (('vrc_tst_scan', False), ('vrc_tst_noscan', True)):
            reactions = self.configured_reactions(
                configured_selection(self.par[option], self.well, []))
            if not reactions:
                continue
            self.opt_products(reactions)
            self.save_products(reactions)
            self.find_scan_coos(reactions)
            invalid = [name for name in reactions
                       if not isinstance(self.scan_reac[name].scan_coo,
                                         (list, tuple, np.ndarray))
                       or len(self.scan_reac[name].scan_coo) != 2]
            if invalid:
                raise ValueError(
                    'VRC scan coordinate is unavailable for: '
                    + ', '.join(invalid))
            self.find_equiv(reactions)
            from kinbot.vrc_multistate import channel_for_reaction, calculate_channel
            scalar = []
            for name in reactions:
                selection = channel_for_reaction(self.par, self.scan_reac[name])
                if selection is None:
                    scalar.append(name)
                else:
                    if noscan:
                        raise ValueError('State-resolved VRC requires correction scans.')
                    calculate_channel(self, name, selection)
            if scalar:
                self.do_scan(scalar, noscan=noscan)
                self.energies(scalar, noscan=noscan)
        if self.scan_reac:
            # A direct KinBot run has all information needed here.  Waiting
            # for a later PES aggregation step left single-well runs without
            # the rotdPy input that the user requested.
            from kinbot.pes import create_rotdpy_inputs
            barrierless = []
            reactant = routing_name(self.well)
            for name, reaction in self.scan_reac.items():
                barrierless.append([
                    reactant, name,
                    [routing_name(product) for product in reaction.products],
                    0.0])
            # Full-PES aggregation owns sampling after all worker corrections
            # are ready. A worker must not run a duplicate ROTD campaign.
            rotd_parameters = dict(self.par)
            if self.par.get('pes'):
                rotd_parameters['rotdpy_run'] = 0
            create_rotdpy_inputs(rotd_parameters, barrierless, [],
                                 correction_root='.')
        return

    def configured_reactions(self, reactions):
        return [reaction.instance_name for reaction in self.well.reac_obj
                if matches_name(reaction.instance_name, reactions)]

    def opt_products(self, reactions):
        """
        reac: name of reaction for which scan is needed
        Find the products based on the reaction name.
        Submit the optimizations for the fragments.
        Does not wait for the jobs to finish.
        """

        prod_chemid = []  # to avoid double submission
        for reac in reactions:
            hit = False  # reaction found in previously explored set
            for ro in self.well.reac_obj:
                if ro.instance_name == reac:
                    hit = True
                    if len(ro.products) != 2:
                        logger.warning(
                            f'There are {len(ro.products)} products for \
                                {reac}, cannot do scan if not bimolecular.')
                    else:
                        # everything is there about the original reactions
                        self.scan_reac[reac] = ro
                        self.scan_reac[reac].usym = [[], []]
                        for prod in self.scan_reac[reac].products:
                            apply_electronic_state_model(prod, self.par)
                            prod.name = routing_name(prod)
                            if routing_key(prod) in prod_chemid:
                                continue
                            else:  # submit at vrc_tst level
                                prod_chemid.append(routing_key(prod))
                                self.qc.qc_vts_frag(prod)
            if not hit:
                logger.warning(f'Reaction {reac} requested \
                               to scan is not found in reaction list.')

        return

    def save_products(self, reactions):
        """
        Wait until all optimizations are done and save the geometries.
        Not everything is checked since no big changes are expected.
        """

        for reac in reactions:
            for prod in self.scan_reac[reac].products:
                # read geom
                job = f'vrctst/{routing_name(prod)}_vts'
                status, geom, atoms = self.qc.get_qc_geom(
                    job,
                    prod.natom,
                    wait=1,
                    allow_error=1,
                    reorder=True)
                geometry_array = np.asarray(geom, dtype=float)
                if (status != 0 or geometry_array.shape != (prod.natom, 3)
                        or not np.isfinite(geometry_array).all()
                        or len(atoms) != prod.natom):
                    raise RuntimeError(
                        f'{job}: VRC fragment optimization did not produce a '
                        'finite, complete geometry.')
                # Publish the fragment geometry only after the complete native
                # result has passed validation.
                prod.geom = geometry_array
                prod.atom = atoms

                if str(self.par.get('pp_orient', 'homo')).casefold() != 'homo':
                    continue
                chk = Path(f'{job}.chk')
                fchk = Path(f'{job}.fchk')
                cube = Path(f'{job}.cube')
                if not fchk.is_file():
                    if not chk.is_file():
                        raise FileNotFoundError(
                            f'{job}: Gaussian checkpoint is unavailable for '
                            'HOMO-oriented VRC pivot points.')
                    result = run(
                        [self.gaussian_tools['formchk'], str(chk), str(fchk)],
                        capture_output=True, text=True, check=False)
                    if (result.returncode or not fchk.is_file()
                            or not fchk.stat().st_size):
                        detail = result.stderr.strip() or result.stdout.strip()
                        raise RuntimeError(
                            f'{job}: formchk failed to create {fchk}.'
                            + (f' {detail}' if detail else ''))
                if not cube.is_file():
                    orbital = 'AMO=HOMO' if prod.mult != 1 else 'MO=HOMO'
                    result = run(
                        [self.gaussian_tools['cubegen'], '1', orbital,
                         str(fchk), str(cube), '-4'],
                        capture_output=True, text=True, check=False)
                    if (result.returncode or not cube.is_file()
                            or not cube.stat().st_size):
                        detail = result.stderr.strip() or result.stdout.strip()
                        raise RuntimeError(
                            f'{job}: cubegen failed to create {cube}.'
                            + (f' {detail}' if detail else ''))
        return

    def find_scan_coos(self, reactions):
        """
        Find bond to be scanned
        This is done in two ways
        - If it's a real well, then we use the bond of hom_sci,
          or the breaking ts bond (only 1 is allowed)
        - If it's a vdw well, then find the closest atoms between
          the two parts using the core function
        """

        for reac in reactions:
            # find the fragments from the IRCs
            self.scan_reac[reac].parts, self.scan_reac[reac].maps = \
                self.scan_reac[reac].irc_prod.start_multi_molecular()
            for part in self.scan_reac[reac].parts:
                apply_electronic_state_model(part, self.par)
                part.name = routing_name(part)
            self.scan_reac[reac].parts[0].characterize()
            self.scan_reac[reac].parts[1].characterize()
            self.match_order(reac)
            self.scan_reac[reac].irc_prod.find_bond()

            if 'hom_sci' in reac:
                ww = self.scan_reac[reac].instance_name.split('_')
                self.scan_reac[reac].scan_coo = [int(ww[-2]) - 1,
                                                 int(ww[-1]) - 1]
                logger.info(f'Bond to be scanned for {reac} is \
                            {np.array(self.scan_reac[reac].scan_coo)+1}')
            elif self.scan_reac[reac].do_vdW:
                self.scan_reac[reac].scan_coo = [None, None]
                self.scan_reac[reac].scan_coo[0], \
                    self.scan_reac[reac].scan_coo[1] = \
                    self.scan_reac[reac].irc_prod.make_extra_bond(
                        self.scan_reac[reac].parts, self.scan_reac[reac].maps)
            else:
                # adding up bond breaks only
                nreacbond = np.sum(np.array([b for bb in
                                             self.scan_reac[reac].ts.reac_bond
                                             for b in bb if b < 0])) / 2
                if nreacbond == 0:
                    logger.warning(f'No bond change detected,\
                                   unable to determine scan coo for {reac}')
                    self.scan_reac[reac].scan_coo = False
                elif nreacbond < -1:
                    logger.warning(f'More than one bonds changed in {reac},\
                                   unable to determine scan coo for {reac}')
                    self.scan_reac[reac].scan_coo = False
                else:
                    hit = False
                    for i, bb in enumerate(self.scan_reac[reac].ts.reac_bond):
                        for j, b in bb:
                            if b < 0:
                                logger.info(f'Bond to be scanned for\
                                            {reac} is {i+1}-{j+1}')
                                self.scan_reac[reac].scan_coo = [i, j]
                                hit = True
                                break
                        if hit:
                            break

        return

    def explicit(self, prod, atomid, equiv, mapping, unique=None):
        '''
        Add user defined reaction centers to the equivalent list to forge
        communication between various entrances
        Now also adds equivalency based on resonance
        stabilized centers - this is automatic
        prod: product st_pt object
        atomid: the scan point's chemid in that product
        equiv: the list of equivalent atoms to be appended
        mapping: each element of this list tells which atom the fragment
        corresponds to in the scan object
        '''
        try:
            for ai in configured_selection(self.par['vrc_tst_scan_reac_cent'], prod, []):
                if ai != atomid:
                    for ii, aii in enumerate(prod.atomid):
                        if aii == ai:
                            equiv.append(mapping[ii])
                            if unique is not None:
                                unique.append([mapping[ii]])
        except KeyError:
            pass

        # look over resonances and automatically add them
        for rad in prod.rads:  # rad is 1 at the radical
            if any(rad == 1):
                for atm in np.where(rad == 1)[0]:
                    if mapping[atm] not in equiv:
                        equiv.append(mapping[list(rad).index(1)])
                        if unique is not None:
                            unique.append([list(rad).index(1)])

        return

    def match_order(self, reac) -> None:
        '''
        Keeps the object to be scanned intact, but rearranges the products
        and the atoms in the products so that:
        1. fragment0 of scan is product0
        2. in both fragments the order of the atoms,
           labeled by atomid, is identical
        It might still scramble chirality,
        but that is perhaps a very rare edge case?
        '''
        # match product ordering with scanned fragments' order
        if not same_species(self.scan_reac[reac].parts[0], self.scan_reac[reac].products[0]):
            self.scan_reac[reac].products = \
                list(reversed(self.scan_reac[reac].products))
        if not same_species(self.scan_reac[reac].parts[0], self.scan_reac[reac].products[0]):
            logger.warning('IRC prod and prod are not the same!')
        if not same_species(self.scan_reac[reac].parts[1], self.scan_reac[reac].products[1]):
            logger.warning('IRC prod and prod are not the same!')
        # match atom ordering of individual products
        # with scanned fragment's atom order
        for fi, frag in enumerate(self.scan_reac[reac].parts):
            self.scan_reac[reac].products[fi].reset_order()
            # This reorders the map so that if fits products and not parts
            self.scan_reac[reac].maps[fi] = reorder_coord(
                mol_A=self.scan_reac[reac].products[fi],
                mol_B=frag,
                map_B=self.scan_reac[reac].maps[fi])
        return

    def do_scan(self, reactions, noscan=False):
        """
        There are two different scans to be made.
        1. relaxed scan: all degrees of freedom are optimized
           except the B-C distance that is scanned.
           Final geometry is saved in each step.
           The user can request that the RMSD < then a
           threshold during optimization..
           Also, if the optimization crashed, the last valid point is taken.
        2. frozen scan: no degrees of freedom are optimized.
           The frozen fragments are oriented in a way the minimizes
           RMSD relative to the relaxed scan.
        Following this, a high-level calculation is done on both set of
        geometries, just single pt.
        All of these calculations are invoked from a single template per point,
        and the points are run in sequence.
        """

        scan_method = str(self.par['vrc_tst_scan_method']).strip()
        if vrc_multireference_energy_contract(scan_method) is not None:
            return self._do_molpro_multireference_scan(
                reactions, noscan=noscan)
        if any(token in scan_method.casefold()
               for token in ('caspt2', 'mrci', 'casscf')):
            raise ValueError(
                f'{scan_method} cannot provide this relaxed VRC path. '
                'Supported Molpro scan methods are CASPT2(ne,no), '
                'CASPT2-F12[A/B](ne,no), MRCI+Q(ne,no), and '
                'MRCI-F12[A/B]+Q(ne,no).')

        status = ['ready'] * len(reactions)
        step = np.zeros(len(reactions), dtype=int)
        geoms = []
        for reac in reactions:
            if not self.scan_reac[reac].do_vdW:
                geoms.append(self.scan_reac[reac].species.geom)
            else:
                geoms.append(self.scan_reac[reac].irc_prod.geom)
        jobs = [''] * len(reactions)
        step0_geoms = [np.array(geoms[ri]) for ri in range(len(reactions))]

        while 1:
            for ri, reac in enumerate(reactions):
                if status[ri] == 'ready' and (step[ri] < \
                   len(self.par['vrc_tst_scan_points']) + 1 or
                   noscan and step[ri] < 1):
                    # shift geometries along the bond to next desired distance
                    # scanning between atoms A and B
                    pos_A = geoms[ri][self.scan_reac[reac].scan_coo[0]]
                    pos_B = geoms[ri][self.scan_reac[reac].scan_coo[1]]
                    dist_AB = np.linalg.norm(pos_B - pos_A)
                    vec_AB = geometry.unit_vector(
                        np.array(pos_B) - np.array(pos_A))
                    # shift so that atom A is at origin
                    geoms[ri] = list(np.array(geoms[ri]) - np.array(pos_A))
                    # stretch frag B along B-A vector
                    if step[ri] < len(self.par['vrc_tst_scan_points']) and\
                       not noscan:
                        shift = vec_AB * (
                            self.par['vrc_tst_scan_points'][step[ri]] -
                            dist_AB)
                        asymptote = False
                    else:
                        shift = vec_AB * (asymptote_distance(self.par) - dist_AB)
                        asymptote = True
                    for mi in self.scan_reac[reac].maps[1]:
                        geoms[ri][mi] = [gi + shift[i]
                                         for i, gi in enumerate(geoms[ri][mi])]
                    new_distAB = np.linalg.norm(
                        geoms[ri][self.scan_reac[reac].scan_coo[0]] -
                        geoms[ri][self.scan_reac[reac].scan_coo[1]])
                    # Temporary fix in case shift is in wrong direction
                    if step[ri] != 0 and step[ri] < len(self.par['vrc_tst_scan_points']):
                        if (new_distAB < dist_AB and\
                           self.par['vrc_tst_scan_points'][step[ri]] > self.par['vrc_tst_scan_points'][step[ri] -1]) or\
                           (new_distAB > dist_AB and\
                           self.par['vrc_tst_scan_points'][step[ri]] < self.par['vrc_tst_scan_points'][step[ri] -1]):
                            for mi in self.scan_reac[reac].maps[1]:
                                geoms[ri][mi] = [gi - 2*shift[i]
                                                for i, gi in enumerate(geoms[ri][mi])]

                    jobs[ri] = self.qc.qc_vts(self.scan_reac[reac],
                                              geoms[ri],
                                              step[ri],
                                              self.scan_reac[reac].equiv,
                                              asymptote,
                                              # needed for alignment of
                                              # rigid fragments later
                                              step0_geoms[ri].tolist()
                                              )
                    logger.info(f'\trunning {jobs[ri]}')
                    status[ri] = 'running'
                elif status[ri] == 'running':
                    geometry_status, geom = \
                        self.qc.get_qc_geom(jobs[ri],
                                            self.scan_reac[reac].species.natom,
                                            allow_error=1)
                    qcst = self.qc.check_qc(jobs[ri])
                    if qcst in ['normal', 'error']:
                        geom = np.asarray(geom, dtype=float)
                        expected = (self.scan_reac[reac].species.natom, 3)
                        if (geometry_status < 0 or geom.shape != expected
                                or not np.isfinite(geom).all()):
                            raise RuntimeError(
                                f'{jobs[ri]}: Gaussian VRC scan point ended '
                                'without a finite geometry.')
                        status[ri] = 'ready'
                        step[ri] += 1
                        geoms[ri] = copy.deepcopy(geom)
                if (step[ri] == len(self.par['vrc_tst_scan_points']) + 1 or
                      (noscan and step[ri] == 1)):
                    status[ri] = 'done'
            if len([st for st in status if st == 'done']) == len(reactions):
                break
            else:
                time.sleep(1)

        return (jobs)

    @staticmethod
    def _shift_scan_geometry(reaction, geom, distance):
        """Move fragment B so the selected atom pair has ``distance``."""
        geom = np.asarray(geom, dtype=float).copy()
        first, second = reaction.scan_coo
        vector = geom[second] - geom[first]
        current = float(np.linalg.norm(vector))
        if not np.isfinite(current) or current <= 1.e-8:
            raise ValueError('VRC scan starts from a zero or nonfinite distance.')
        geom -= geom[first]
        shift = vector / current * (float(distance) - current)
        geom[np.asarray(reaction.maps[1], dtype=int)] += shift
        observed = float(np.linalg.norm(geom[first] - geom[second]))
        if not np.isclose(observed, distance, atol=1.e-8, rtol=0.):
            # Product ordering can put the second scanned atom in the first
            # map.  The opposite translation is then the physically correct
            # one and is deterministic.
            geom[np.asarray(reaction.maps[1], dtype=int)] -= 2. * shift
            observed = float(np.linalg.norm(geom[first] - geom[second]))
        if not np.isclose(observed, distance, atol=1.e-8, rtol=0.):
            raise RuntimeError('Unable to construct the requested VRC distance.')
        return geom

    def _molpro_multireference_scan_task(
            self, reaction, geom, distance, logical_name, *,
            state_model_override=None):
        """Run one constrained Molpro multireference scan point."""
        from kinbot.anl.dispatch import (
            _load, preflight_until_scheduler_ready, prepare)

        constraints = {'fixed_bonds': [list(map(int, reaction.scan_coo))],
                       'minimum_bonds': []}
        for first in reaction.equiv[0]:
            for second in reaction.equiv[1]:
                pair = [int(first), int(second)]
                if pair == constraints['fixed_bonds'][0] or \
                        pair[::-1] == constraints['fixed_bonds'][0]:
                    continue
                constraints['minimum_bonds'].append(
                    [pair[0], pair[1], float(distance)])
        method = str(self.par['vrc_tst_scan_method']).strip()
        basis = str(self.par['vrc_tst_scan_basis']).strip()
        settings = roots = None
        if correction_model(self.par) == 'spin_splitting':
            settings, roots = self._spin_context(reaction)
        elif len(reaction.products) == 2:
            roots = self._state_average_roots(reaction.products,
                [int(reaction.species.mult)], self.par)
        active_space = active_space_settings(
            self.par, reaction.species.atom, reaction.species.charge)
        level = resolve_vrc_level(
            {'method': method, 'basis': basis}, active_space)
        method, basis = level['method'], level['basis']
        state_model = (copy.deepcopy(state_model_override)
                       if state_model_override is not None else
                       self.par.get('vrc_tst_state_model')
                       if active_space is None else
                       resolved_state_model(
                           self.par, reaction.species.atom,
                           reaction.species.charge, settings, roots=roots,
                           basis=basis))
        molecule = {
            'symbols': list(map(str, reaction.species.atom)),
            'positions': np.asarray(geom, dtype=float).tolist(),
            'charge': _native_state_integer(
                reaction.species.charge, 'Molecular charge'),
            'multiplicity': _native_state_integer(
                reaction.species.mult, 'Molecular multiplicity'),
        }
        task = {
            'id': 'relaxed_scan', 'kind': 'ase_optimize',
            'geometry_from': 'initial', 'geometry_output': 'final.xyz',
            'profile': {
                'calculator': 'molpro', 'method': method, 'basis': basis,
                'command': 'molpro', 'optimizer': 'sella',
                'frequency_mode': 'ase_forces',
                'calculator_kwargs': {
                    'caspt2_shift': self.par.get(
                        'vrc_tst_scan_shift', 0.1),
                    'caspt2_ipea': self.par.get(
                        'vrc_tst_scan_ipea', 0.25),
                    'state_model': state_model,
                    'gradient_mode': self.par.get(
                        'vrc_tst_scan_gradient', 'auto'),
                    'scf_convergence': self.par.get(
                        'vrc_tst_scf_convergence', 1.e-10),
                    'correlation_convergence': self.par.get(
                        'vrc_tst_correlation_convergence', 1.e-8),
                },
                'label': f'RS2-{method.upper()}/{basis}',
            },
            'optimizer': {
                'fmax': 0.03, 'steps': 100, 'order': 0,
                'sella_kwargs': {'internal': False},
                'constraints': constraints,
            },
            'resources': vrc_dispatch_resources(self.par),
        }
        fingerprint = hashlib.sha256(json.dumps({
            'molecule': molecule, 'task': task,
        }, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        safe_name = re.sub(r'[^A-Za-z0-9_.-]', '_', logical_name)
        run_dir = (Path('vrctst/molpro_scan').resolve() /
                   f'{safe_name}_{fingerprint[:16]}')
        spec = {
            'schema': 1, 'name': f'vrc-caspt2-{safe_name}',
            'molecule': molecule,
            'limits': {'max_nodes': 1}, 'tasks': [task],
            'intent': {
                'cap_walltime_to_explicit_partition': True,
                'vrc_scan_fingerprint': fingerprint,
                'vrc_scan_distance_angstrom': float(distance),
            },
        }
        if not run_dir.exists():
            prepare(spec, run_dir)
        else:
            _, stored, _ = _load(run_dir)
            if stored.get('intent', {}).get(
                    'vrc_scan_fingerprint') != fingerprint:
                raise RuntimeError(
                    f'VRC scan provenance collision at {run_dir}.')
        preflight_until_scheduler_ready(run_dir, interval=20.)
        _, _, state = _load(run_dir)
        if state['tasks']['relaxed_scan']['status'] != 'complete':
            returncode, _ = _drive_vrc_dispatch(
                run_dir, ('relaxed_scan',), interval=20.)
            if returncode:
                _, _, state = _load(run_dir)
                error = state['tasks']['relaxed_scan'].get(
                    'error', 'unknown task failure')
                raise RuntimeError(
                    f'Molpro multireference VRC scan failed at {distance:g} '
                    f'angstrom; inspect {run_dir}: {error}')
        task_dir = run_dir / 'tasks' / 'relaxed_scan'
        execution = json.loads((task_dir / 'execution.json').read_text())
        if execution.get('status') != 'executed':
            raise RuntimeError(f'Incomplete VRC scan execution at {run_dir}.')
        atoms = read(task_dir / 'final.xyz')
        observed = float(np.linalg.norm(
            atoms.positions[reaction.scan_coo[0]] -
            atoms.positions[reaction.scan_coo[1]]))
        if not np.isclose(observed, distance, atol=0.01, rtol=0.):
            raise RuntimeError(
                f'Molpro VRC point lost its distance constraint: '
                f'{observed:.6f} versus {distance:.6f} angstrom.')
        # A fixed association coordinate must not be satisfied by breaking an
        # internal fragment bond.  Check the original connectivity after each
        # optimization before publishing the geometry for correction jobs.
        broken = self._broken_scan_fragment_bonds(reaction, atoms.positions)
        if broken:
            detail = ', '.join(
                f'{first + 1}-{second + 1}: {initial:.4f}->{final:.4f} A'
                for first, second, initial, final in broken)
            raise RuntimeError(
                'Molpro VRC relaxed scan broke internal fragment bond(s): '
                + detail)
        energy_ev = float(execution['details']['energy_ev'])
        recorded_model = execution.get('details', {}).get(
            'orbital_checkpoint_model')
        if isinstance(recorded_model, dict):
            selection = recorded_model.get('orbital_selection', {})
            checkpoint = Path(str(selection.get('path', '')))
            if (selection.get('mode') != 'molpro_record'
                    or not checkpoint.is_file()
                    or hashlib.sha256(checkpoint.read_bytes()).hexdigest()
                    != selection.get('sha256')):
                raise RuntimeError(
                    f'Molpro VRC point has an invalid saved orbital record: '
                    f'{run_dir}.')
            return atoms, energy_ev, run_dir, execution, recorded_model
        generated = execution.get('details', {}).get('generated_files', [])
        checkpoints = [task_dir / Path(name).name for name in generated
                       if str(name).endswith('.wfu')]
        checkpoints = [path for path in checkpoints
                       if path.is_file() and path.stat().st_size]
        if not checkpoints:
            raise RuntimeError(
                f'Molpro VRC point did not publish an orbital checkpoint: '
                f'{run_dir}.')
        checkpoint = checkpoints[-1]
        next_state_model = copy.deepcopy(state_model)
        if next_state_model is None:
            raise RuntimeError(
                'A multireference scan completed without an electronic-state '
                'model; orbital continuity cannot be certified.')
        next_state_model['orbital_selection'] = {
            'mode': 'molpro_record',
            'path': str(checkpoint.resolve()),
            'record': '2140.2',
            'sha256': hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
        }
        return atoms, energy_ev, run_dir, execution, next_state_model

    @staticmethod
    def _broken_scan_fragment_bonds(reaction, positions):
        """Return fragment bonds stretched beyond the isolated geometry."""
        positions = np.asarray(positions, dtype=float)
        fragment_by_atom = {
            int(atom): fragment
            for fragment, mapping in enumerate(reaction.maps)
            for atom in mapping
        }
        broken = []
        for first, second in reaction.irc_prod.bondlist:
            first, second = int(first), int(second)
            fragment = fragment_by_atom.get(first)
            if fragment is None or fragment != fragment_by_atom.get(second):
                continue
            mapping = list(map(int, reaction.maps[fragment]))
            product = np.asarray(reaction.products[fragment].geom, dtype=float)
            initial = float(np.linalg.norm(
                product[mapping.index(first)] - product[mapping.index(second)]))
            final = float(np.linalg.norm(
                positions[first] - positions[second]))
            if initial <= 1.e-8 or final > 1.2 * initial:
                broken.append((first, second, initial, final))
        return broken

    @staticmethod
    def _frozen_scan_geometry(reaction, relaxed):
        """Orient optimized isolated fragments onto a relaxed VRC point."""
        import rmsd

        maps = [list(map(int, mapping)) for mapping in reaction.maps]
        scan = list(map(int, reaction.scan_coo))
        relaxed = np.asarray(relaxed, dtype=float)
        result = np.zeros_like(relaxed)
        for fragment_index in range(2):
            mapping = maps[fragment_index]
            relaxed_part = relaxed[mapping]
            product = reaction.products[fragment_index]
            frozen = np.asarray(product.geom, dtype=float).copy()
            reactive = scan[fragment_index]
            if reactive not in mapping:
                reactive = scan[1 - fragment_index]
            pivot = mapping.index(reactive)
            weights = np.asarray(product.bond01[pivot], dtype=float).copy()
            weights[weights == 0.] = 0.02
            weights[pivot] = 1.
            fitted, _ = rmsd.kabsch_weighted_fit(
                frozen, relaxed_part, W=weights)
            fitted -= fitted[pivot] - relaxed_part[pivot]
            result[mapping] = fitted
        if not np.isfinite(result).all():
            raise RuntimeError('Frozen VRC fragment alignment is nonfinite.')
        return result

    def _publish_scan_geometry(self, reaction, name, atoms, energy_ev,
                               run_dir, *, frozen):
        """Publish a verified scan geometry under the legacy database key."""
        atoms = atoms.copy()
        charges = np.zeros(len(atoms))
        charges[0] = int(reaction.species.charge)
        atoms.set_initial_charges(charges)
        connect('kinbot.db').write(atoms, name=f'vrctst/{name}', data={
            'energy': float(energy_ev), 'status': 'normal',
            'vrc_scan_backend': 'molpro',
            'vrc_scan_method': str(self.par['vrc_tst_scan_method']),
            'vrc_scan_basis': str(self.par['vrc_tst_scan_basis']),
            'vrc_scan_role': 'frozen' if frozen else 'relaxed',
            'vrc_scan_dispatch': str(run_dir),
        })
        path = Path('vrctst') / f'{name}.log'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            'Molpro multireference VRC scan geometry accepted\n'
            f'dispatch={run_dir}\ndone\n')

    def _do_molpro_multireference_scan(self, reactions, noscan=False):
        """Compute a relaxed VRC path with the requested Molpro MR method."""
        Path('vrctst/molpro_scan').mkdir(parents=True, exist_ok=True)
        completed = []
        for reac in reactions:
            reaction = self.scan_reac[reac]
            settings = roots = None
            if correction_model(self.par) == 'spin_splitting':
                settings, roots = self._spin_context(reaction)
            elif len(reaction.products) == 2:
                roots = self._state_average_roots(reaction.products,
                    [int(reaction.species.mult)], self.par)
            active_space = active_space_settings(
                self.par, reaction.species.atom, reaction.species.charge)
            state_model = (self.par.get('vrc_tst_state_model')
                           if active_space is None else
                           resolved_state_model(
                               self.par, reaction.species.atom,
                               reaction.species.charge, settings,
                               roots=roots,
                               basis=self.par['vrc_tst_scan_basis']))
            geom = np.asarray(
                reaction.irc_prod.geom if reaction.do_vdW
                else reaction.species.geom, dtype=float)
            separation = asymptote_distance(self.par)
            distances = ([separation] if noscan else
                         [*map(float, self.par['vrc_tst_scan_points']), separation])
            for step, distance in enumerate(distances):
                geom = self._shift_scan_geometry(reaction, geom, distance)
                base = (f'{reac}_vts_pt_asymptote'
                        if step == len(distances) - 1
                        else f'{reac}_vts_pt{step:02d}')
                if step == len(distances) - 1:
                    # At the asymptote there is no intermolecular geometry to
                    # relax.  The correction energies are evaluated later at
                    # this explicitly separated structure.
                    relaxed = Atoms(reaction.species.atom, positions=geom)
                    energy_ev = 0.
                    run_dir = Path('vrctst/molpro_scan/asymptote_identity')
                else:
                    kwargs = ({} if state_model is None else
                              {'state_model_override': state_model})
                    result = self._molpro_multireference_scan_task(
                        reaction, geom, distance, base, **kwargs)
                    if len(result) == 5:
                        relaxed, energy_ev, run_dir, _, state_model = result
                    else:
                        # Preserve the established callback contract for
                        # custom direct scans which do not use an orbital
                        # checkpoint model.
                        relaxed, energy_ev, run_dir, _ = result
                    geom = relaxed.positions.copy()
                self._publish_scan_geometry(
                    reaction, base, relaxed, energy_ev, run_dir, frozen=False)
                frozen_geom = self._frozen_scan_geometry(
                    reaction, relaxed.positions)
                frozen = Atoms(reaction.species.atom, positions=frozen_geom)
                self._publish_scan_geometry(
                    reaction, base + '_fr', frozen, energy_ev, run_dir,
                    frozen=True)
                logger.info(
                    'Accepted %s at RS2-%s/%s (%.3f angstrom)', base,
                    self.par['vrc_tst_scan_method'],
                    self.par['vrc_tst_scan_basis'], distance)
                completed.append(str(run_dir))
        return completed

    def energies(self, reactions, noscan=False):
        """Run and read the VRC Molpro correction calculations.

        These calculations used to be written into a batch file and left for
        a manual submission that the direct KinBot workflow never resumed.
        Stage them through the restartable exclusive-node dispatcher instead,
        then create the correction records before returning.
        """
        db = connect('kinbot.db')
        records = {}
        pending = {}
        model = correction_model(self.par)
        components_by_reaction = {}
        context_by_reaction = {}
        for reac in reactions:
            reaction = self.scan_reac[reac]
            settings = roots = None
            if model == 'spin_splitting':
                fragment_mults = [int(product.mult)
                                  for product in reaction.products]
                parent_mult = int(reaction.species.mult)
                if len(fragment_mults) != 2:
                    raise ValueError(
                        'The spin-splitting VRC correction requires exactly '
                        'two separated fragments.')
                settings, roots = self._spin_context(reaction)
                declared = settings.get('fragment_multiplicities')
                if declared is not None and fragment_mults != declared:
                    raise ValueError(
                        'The spin-splitting fragment multiplicities differ '
                        'from the accepted reaction products.')
                allowed = coupled_multiplicities(*fragment_mults)
                if parent_mult != settings['target_multiplicity']:
                    raise ValueError(
                        'The spin-splitting target multiplicity differs from '
                        'the associated molecular state.')
                if (parent_mult not in allowed
                        or settings['reference_multiplicity'] not in allowed):
                    raise ValueError(
                        'The requested target/reference spin surfaces are '
                        'not allowed by the fragment multiplicities.')
            elif len(reaction.products) == 2:
                roots = self._state_average_roots(reaction.products,
                    [int(reaction.species.mult)], self.par)
            context_by_reaction[reac] = (settings, roots)
            parent = reaction.species
            components = correction_components(
                self.par, symbols=getattr(parent, 'atom', None),
                charge=getattr(parent, 'charge', 0), settings=settings,
                state_roots=roots)
            components_by_reaction[reac] = components
            records[reac] = []
            ndist = 1 if noscan else len(self.par['vrc_tst_scan_points']) + 1
            for step in range(ndist):
                base = (f'{reac}_vts_pt{str(step).zfill(2)}'
                        if step < ndist - 1
                        else f'{reac}_vts_pt_asymptote')
                for component in components:
                    job = base + component['suffix']
                    geometry_job = (base + '_fr'
                                    if (component['geometry'] == 'frozen'
                                        or step == ndist - 1) else base)
                    rows = list(db.select(name=f'vrctst/{geometry_job}',
                                          sort='-id', limit=1))
                    if not rows:
                        raise RuntimeError('Missing accepted VRC geometry for '
                                           f'{geometry_job}.')
                    last_row = rows[0]
                    scan_spec = StationaryPoint.from_ase_atoms(
                        last_row.toatoms(),
                        charge=int(getattr(
                            self.scan_reac[reac].species, 'charge', 0)),
                        mult=int(self.scan_reac[reac].species.mult))
                    scan_spec.characterize()

                    molp = Molpro(scan_spec, self.par)
                    molp.create_molpro_input(
                        name=job, VTS=True,
                        sample=component['role'].startswith('sampling'),
                        vrc_level=component)
                    if self._vrc_output_matches_input(job):
                        e_stat, e = molp.get_molpro_energy(
                            key=self.par['vrc_tst_scan_molpro_key'],
                            name=f'{job}', VTS=True)
                    else:
                        e_stat, e = 0, -1.
                    logger.debug(f'{job}, {e_stat}, {e}')
                    if not e_stat:
                        pending[job] = scan_spec
                    records[reac].append((component['role'], job, molp))

        if pending:
            self._dispatch_molpro_corrections(pending)

        for reac in reactions:
            reaction = self.scan_reac[reac]
            settings, roots = context_by_reaction[reac]
            components = components_by_reaction[reac]
            parent = reaction.species
            active_space = active_space_settings(
                self.par, getattr(parent, 'atom', None),
                getattr(parent, 'charge', 0))
            scan_method = self.par.get('vrc_tst_scan_method', 'unspecified')
            if active_space is not None:
                scan_method = method_with_active_space(
                    scan_method, active_space)
            scan_basis = self.par.get('vrc_tst_scan_basis', 'unspecified')
            scan_contract = vrc_multireference_energy_contract(scan_method)
            energies = {component['role']: [] for component in components}
            native_components = {component['role']: []
                                 for component in components}
            for role, job, molp in records[reac]:
                e_stat, energy = molp.get_molpro_energy(
                    key=self.par['vrc_tst_scan_molpro_key'], name=job,
                    VTS=True)
                if not e_stat or not self._vrc_output_matches_input(job):
                    raise RuntimeError(f'VRC Molpro result {job} has no '
                                       f'{self.par["vrc_tst_scan_molpro_key"]} energy.')
                energies[role].append(energy)
                directory = Path('vrctst/molpro')
                native_components[role].append({
                    'job': job,
                    'input_sha256': hashlib.sha256(
                        (directory / f'{job}.inp').read_bytes()).hexdigest(),
                    'output_sha256': hashlib.sha256(
                        (directory / f'{job}.out').read_bytes()).hexdigest(),
                    'energy_hartree': energy,
                    'energy_contract':
                        vrc_multireference_energy_contract(
                            next(item['method'] for item in components
                                 if item['role'] == role)),
                })

            dist = ([asymptote_distance(self.par)] if noscan
                    else self.par['vrc_tst_scan_points'] + [asymptote_distance(self.par)])
            e_samp, sampling_cbs = assemble_sampling_reference(
                dist, energies, sampling_settings(
                    self.par, active_space=active_space))
            spin_splitting = None
            if model == 'direct':
                from kinbot.vrc_correction import assemble_trusted_reference
                e_high, trusted_cbs = assemble_trusted_reference(
                    dist, energies, self.par, active_space=active_space)
            else:
                e_high, spin_splitting = assemble_spin_splitting_reference(
                    dist, energies, settings)
            ens = []  # energies for sample and high in kcal/mol
            asyms = []  # asymptotic energies in hartree
            for curve in (e_samp, e_high):
                if len(curve) != len(dist):
                    raise RuntimeError(f'Incomplete VRC energy series for {reac}.')
                asyms.append(curve[-1])
                ens.append(list((np.array(curve) - curve[-1])
                                * constants.AUtoKCAL))

            # Create scan references between all equivalent atoms:
            scan_ref = []

            for i in self.scan_reac[reac].usym[0][0]:
                for j in self.scan_reac[reac].usym[1][0]:
                    scan_ref.append([i, j])

            # Create list of reactive atoms (fragment indexed)
            ra: list[list[int]] = [[], []]
            for i in range(2):
                for j in self.scan_reac[reac].equiv[i]:
                    ra[i].append(
                        np.where(self.scan_reac[reac].maps[i] == j)[0][0])

            if not noscan:
                create_matplotlib_graph(x=dist,
                                        data=ens,
                                        name=f'{reac}',
                                        x_label=f"{reac}",
                                        y_label="Energy (kcal/mol)",
                                        data_legends=['sample', 'high'])

            smallest = np.linalg.norm(
                self.well.geom[self.scan_reac[reac].equiv[0][0]] -
                self.well.geom[self.scan_reac[reac].equiv[1][0]])

            corr: dict[str, Any] = {
                'dist': dist,
                'e_samp': ens[0],
                'e_high': ens[1],
                'scan_ref': scan_ref,
                'ra': ra,
                'smallest': smallest,
                'unique': self.scan_reac[reac].usym,
                'e_inf_samp': asyms[0],
                'e_inf_high': asyms[1],
                'correction_model': model,
                'levels': expected_levels(
                    self.par, symbols=getattr(parent, 'atom', None),
                    charge=getattr(parent, 'charge', 0), settings=settings,
                    state_roots=roots),
                'native_components': native_components,
                'vrc_molecule': {
                    'symbols': list(map(str, getattr(parent, 'atom', scan_spec.atom))),
                    'charge': int(getattr(parent, 'charge', 0)),
                    'multiplicity': int(parent.mult),
                },
                'sampling_atom_order': [int(index)
                                        for mapping in reaction.maps
                                        for index in mapping],
                'state_average_roots': roots,
                'relaxed_scan_level': {
                    'backend': ('molpro' if scan_contract is not None
                        else getattr(self.qc, 'qc', 'configured')),
                    'method': scan_method,
                    'basis': scan_basis,
                    'gradient': (
                        'Molpro RS2 analytic CASPT2'
                        if (scan_contract is not None and
                            scan_contract['method'] == 'CASPT2')
                        else ('Molpro numerical multireference gradient'
                              if scan_contract is not None
                              else 'configured legacy scan')),
                },
                'frags_atom': [list(self.scan_reac[reac].products[0].atom),
                               list(self.scan_reac[reac].products[1].atom)],
                'frags_geom': [self.scan_reac[reac].products[0].geom,
                               self.scan_reac[reac].products[1].geom],
                'frags_mult': [self.scan_reac[reac].products[0].mult,
                                self.scan_reac[reac].products[1].mult],
                'frags_routing': [fragment_routing_state(product,
                                                          f'fragment_{index}')
                                  for index, product in enumerate(
                                      self.scan_reac[reac].products)]
                }
            if spin_splitting is not None:
                corr['spin_splitting'] = spin_splitting
            if sampling_cbs is not None:
                corr['sampling_cbs'] = sampling_cbs
            if model == 'direct' and trusted_cbs is not None:
                corr['trusted_cbs'] = trusted_cbs

            with open(f'vrctst/corr_{reac}.json', 'w',
                      encoding='utf-8') as f:
                json.dump(corr, f, ensure_ascii=False, indent=4,
                          cls=NpEncoder)
        return

    @staticmethod
    def _vrc_output_matches_input(job):
        """Require a VRC output to be bound to its exact generated input."""
        directory = Path('vrctst/molpro')
        input_path = directory / f'{job}.inp'
        output_path = directory / f'{job}.out'
        fingerprint = directory / f'{job}.input.sha256'
        if not input_path.is_file() or not output_path.is_file() \
                or not fingerprint.is_file():
            return False
        expected = hashlib.sha256(input_path.read_bytes()).hexdigest()
        return fingerprint.read_text().strip() == expected

    def _dispatch_molpro_corrections(self, pending):
        """Execute missing VRC correction points with the ANL dispatcher."""
        from kinbot.anl.dispatch import (
            _load, preflight_until_scheduler_ready, prepare)

        first = next(iter(pending.values()))
        molecule = {
            'symbols': list(first.atom),
            'positions': np.asarray(first.geom).tolist(),
            # StationaryPoint.from_ase_atoms obtains charge by summing ASE's
            # floating-point initial charges.  Dispatcher specifications are
            # JSON data and require native integers for molecular state.
            'charge': _native_state_integer(first.charge,
                                            'Molecular charge'),
            'multiplicity': _native_state_integer(
                first.mult, 'Molecular multiplicity'),
        }
        tasks = []
        for job in sorted(pending):
            input_path = Path('vrctst/molpro') / f'{job}.inp'
            tasks.append({
                'id': job, 'kind': 'external', 'backend': 'molpro',
                'geometry_from': 'initial',
                'resources': vrc_dispatch_resources(self.par),
                'input_name': f'{job}.inp',
                'input_template': input_path.read_text(),
                'command': ['molpro', '-n', '{cores}', '-m',
                            '{molpro_stack_mw}', '{input}'],
                'stdout': 'launcher.stdout', 'stderr': 'launcher.stderr',
                'required_outputs': [f'{job}.out'],
                'success_marker': {
                    'file': f'{job}.out',
                    'contains': 'Molpro calculation terminated'},
            })
        spec = {
            'schema': 1, 'name': 'kinbot-vrc-correction-potentials',
            'molecule': molecule,
            'limits': {'max_nodes': self.par['vrc_tst_max_nodes']},
            'tasks': tasks,
            'intent': {
                # Keep an explicitly selected short validation partition.
                # Its finite limit is part of the execution policy rather
                # than a reason to route work to another partition.
                'cap_walltime_to_explicit_partition': True,
            },
        }
        base_dir = Path('vrctst/molpro').resolve()
        run_dir = base_dir / 'dispatch'
        if not run_dir.exists():
            prepare(spec, run_dir)
        else:
            _, stored, _ = _load(run_dir)
            wanted = {task['id']: task['input_template'] for task in tasks}
            present = {task['id']: task['input_template']
                       for task in stored['tasks']}
            if present != wanted:
                digest = hashlib.sha256(json.dumps(
                    wanted, sort_keys=True).encode()).hexdigest()[:12]
                run_dir = base_dir / f'dispatch_{digest}'
                if not run_dir.exists():
                    prepare(spec, run_dir)
        (base_dir / 'current_dispatch.json').write_text(json.dumps({
            'schema': 1, 'run_dir': str(run_dir),
        }, indent=2) + '\n')
        preflight_until_scheduler_ready(run_dir, interval=20.)
        returncode, _ = _drive_vrc_dispatch(
            run_dir, pending, interval=20.)
        if returncode:
            _, _, state = _load(run_dir)
            failures = []
            for job in sorted(pending):
                entry = state['tasks'].get(job, {})
                if entry.get('status') != 'failed':
                    continue
                outcome = run_dir / 'tasks' / job / 'execution.json'
                try:
                    error = json.loads(outcome.read_text()).get('error')
                except (OSError, json.JSONDecodeError):
                    error = entry.get('error')
                failures.append(f'{job}: {error or "unknown task failure"}')
            detail = ' | '.join(failures)
            raise RuntimeError('One or more VRC Molpro correction jobs failed; '
                               f'inspect {run_dir}.'
                               + (f' {detail}' if detail else ''))
        for job in pending:
            source = run_dir / 'tasks' / job / f'{job}.out'
            destination = base_dir / source.name
            shutil.copyfile(source, destination)
            input_path = base_dir / f'{job}.inp'
            (base_dir / f'{job}.input.sha256').write_text(
                hashlib.sha256(input_path.read_bytes()).hexdigest() + '\n')

    def find_equiv(self, reactions):
        for reac in reactions:
            # determine equivalent atoms
            equiv_A = []
            equiv_B = []
            if self.scan_reac[reac].scan_coo[0] in \
                self.scan_reac[reac].maps[0]:
                # find index of self.scan_reac[reac].scan_coo[0]
                # in prod0 and give its atomid
                index_A = np.where(
                    self.scan_reac[reac].maps[0] ==
                    self.scan_reac[reac].scan_coo[0])[0][0]
                self.scan_reac[reac].usym[0].append([index_A])
                atomid_A = self.scan_reac[reac].\
                    products[0].atomid[index_A]
                for ii, mi in enumerate(
                    self.scan_reac[reac].maps[0]):
                    if self.scan_reac[reac].products[0].\
                        atomid[ii] == atomid_A:
                        equiv_A.append(mi)
                self.explicit(
                    prod=self.scan_reac[reac].products[0],
                    atomid=atomid_A,
                    equiv=equiv_A,
                    mapping=self.scan_reac[reac].maps[0],
                    unique=self.scan_reac[reac].usym[0])
                index_B = np.where(self.scan_reac[reac].maps[1] ==
                                    self.scan_reac[reac].scan_coo[1]
                                    )[0][0]
                self.scan_reac[reac].usym[1].append([index_B])
                atomid_B = self.scan_reac[reac].\
                    products[1].atomid[index_B]
                for ii, mi in enumerate(
                    self.scan_reac[reac].maps[1]):
                    if self.scan_reac[reac].products[1].\
                        atomid[ii] == atomid_B:
                        equiv_B.append(mi)
                self.explicit(
                    prod=self.scan_reac[reac].products[1],
                    atomid=atomid_B,
                    equiv=equiv_B,
                    mapping=self.scan_reac[reac].maps[1],
                    unique=self.scan_reac[reac].usym[1])
                # Complete self.scan_reac[reac].usym to contain equivalent atoms
                for fnum, frag_ra in enumerate(self.scan_reac[reac].usym):
                    # ura is a list with the index of a single unique atom
                    for ura in frag_ra:
                        uaid = self.scan_reac[reac].products[fnum].\
                                atomid[ura[0]]
                        for ii, aid in enumerate(
                            self.scan_reac[reac].products[fnum].atomid):
                            if aid == uaid and \
                                ii not in ura:
                                ura.append(ii)
            else:
                index_A = np.where(
                    self.scan_reac[reac].maps[1] ==
                    self.scan_reac[reac].scan_coo[0])[0][0]
                self.scan_reac[reac].usym[0].append([index_A])
                atomid_A = self.scan_reac[reac].products[1].\
                    atomid[index_A]
                for ii, mi in enumerate(
                    self.scan_reac[reac].maps[1]):
                    if self.scan_reac[reac].products[1].\
                        atomid[ii] == atomid_A:
                        equiv_A.append(mi)
                self.explicit(
                    prod=self.scan_reac[reac].products[1],
                    atomid=atomid_A,
                    equiv=equiv_A,
                    mapping=self.scan_reac[reac].maps[1],
                    unique=self.scan_reac[reac].usym[0])
                index_B = np.where(
                    self.scan_reac[reac].maps[0] ==
                    self.scan_reac[reac].scan_coo[1])[0][0]
                self.scan_reac[reac].usym[1].append([index_B])
                atomid_B = self.scan_reac[reac].\
                    products[0].atomid[index_B]
                for ii, mi in enumerate(
                    self.scan_reac[reac].maps[0]):
                    if self.scan_reac[reac].\
                        products[0].atomid[ii] == atomid_B:
                        equiv_B.append(mi)
                self.explicit(
                    prod=self.scan_reac[reac].products[0],
                    atomid=atomid_B,
                    equiv=equiv_B,
                    mapping=self.scan_reac[reac].maps[0],
                    unique=self.scan_reac[reac].usym[1])
                # Complete self.scan_reac[reac].usym to contain equivalent atoms
                for fnum, frag_ra in enumerate(self.scan_reac[reac].usym):
                    # ura is a list with the index of a single unique atom
                    if fnum == 0:
                        idx = 1
                    elif fnum == 1:
                        idx = 0
                    for ura in frag_ra:
                        uaid = self.scan_reac[reac].products[idx].\
                                atomid[ura[0]]
                        for ii, aid in enumerate(
                            self.scan_reac[reac].products[idx].atomid):
                            if aid == uaid and \
                                ii not in ura:
                                ura.append(ii)
            self.scan_reac[reac].equiv = [equiv_A, equiv_B]
