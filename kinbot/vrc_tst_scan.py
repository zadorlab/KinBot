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
import subprocess
import sys
import hashlib
from pathlib import Path

from shutil import which
from subprocess import Popen, PIPE
from kinbot.utils import reorder_coord
from kinbot.stationary_pt import StationaryPoint
from kinbot import geometry
from kinbot.molpro import Molpro
from kinbot.utils import create_matplotlib_graph, NpEncoder
from kinbot import constants

logger = logging.getLogger('KinBot')


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

    def calculate_correction_potentials(self):
        """
        The main driver for scanning the potential.
        """
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
            self.find_equiv(reactions)
            self.do_scan(reactions, noscan=noscan)
            self.energies(reactions, noscan=noscan)
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
            create_rotdpy_inputs(self.par, barrierless, [],
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
                # update geom in self.scan_reac[reac] values,
                # this already saves it
                prod.geom = geom
                prod.atom = atoms
                if status == 0:
                    # Check if commands are available
                    commands = ['formchk', 'cubegen']
                    missing_com = False
                    for com in commands:
                        if which(com) is None:
                            logger.warning(f'Make sure the command {com} is installed.')
                            missing_com = True
                    if not missing_com:
                        # Create the formchk
                        if not os.path.isfile(f'{job}.fchk'):
                            if os.path.isfile(f'{job}.chk'):
                                command = ['formchk', f'{job}.chk', f'{job}.fchk']
                                process = Popen(
                                    args=command,
                                    shell=False,
                                    stdout=PIPE,
                                    stdin=PIPE,
                                    stderr=PIPE)
                                out, err = process.communicate()
                                out = out.decode()
                            else:
                                logger.warning(f"Could not create formatted checkpoint file for {job} as checkpoint file doesn't exist.")

                        # Create the cubegen file
                        if not os.path.isfile(f'{job}.cube'):
                            if os.path.isfile(f'{job}.fchk'):
                                if prod.mult != 1:
                                    command = ['cubegen', '1', 'AMO=HOMO', f'{job}.fchk', f'{job}.cube', '-4']
                                else:
                                    command = ['cubegen', '1', 'MO=HOMO', f'{job}.fchk', f'{job}.cube', '-4']
                                process = Popen(
                                    command,
                                    shell=False,
                                    stdout=PIPE,
                                    stdin=PIPE,
                                    stderr=PIPE)
                                out, err = process.communicate()
                                out = out.decode()
                            else:
                                logger.warning(f"Cannnot create cube file for {job} as formated checkpoint file doesn't exist.")
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
                        shift = vec_AB * (30. - dist_AB)
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
                    _, geom = \
                        self.qc.get_qc_geom(jobs[ri],
                                            self.scan_reac[reac].species.natom,
                                            allow_error=1)
                    qcst = self.qc.check_qc(jobs[ri])
                    if qcst in ['normal', 'error']:
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
        for reac in reactions:
            records[reac] = []
            ndist = 1 if noscan else len(self.par['vrc_tst_scan_points']) + 1
            for step in range(ndist):
                for sample in [True, False]:
                    if sample:
                        if step < ndist - 1:
                            job = f'{reac}_vts_pt{str(step).zfill(2)}_fr'
                        else:
                            job = f'{reac}_vts_pt_asymptote_fr'
                    else:
                        if step < ndist - 1:
                            job = f'{reac}_vts_pt{str(step).zfill(2)}'
                        else:
                            # Both theories use the frozen asymptotic geometry,
                            # but need distinct inputs and native outputs.
                            job = f'{reac}_vts_pt_asymptote'
                    geometry_job = (job if sample or not job.endswith(
                        '_asymptote') else job + '_fr')
                    rows = list(db.select(name=f'vrctst/{geometry_job}',
                                          sort='-id', limit=1))
                    if not rows:
                        raise RuntimeError('Missing accepted VRC geometry for '
                                           f'{geometry_job}.')
                    last_row = rows[0]
                    scan_spec = StationaryPoint.from_ase_atoms(
                        last_row.toatoms())
                    scan_spec.characterize()

                    molp = Molpro(scan_spec, self.par)
                    molp.create_molpro_input(name=job, VTS=True, sample=sample)
                    if self._vrc_output_matches_input(job):
                        e_stat, e = molp.get_molpro_energy(
                            key=self.par['vrc_tst_scan_molpro_key'],
                            name=f'{job}', VTS=True)
                    else:
                        e_stat, e = 0, -1.
                    logger.debug(f'{job}, {e_stat}, {e}')
                    if not e_stat:
                        pending[job] = scan_spec
                    records[reac].append((sample, job, molp))

        if pending:
            self._dispatch_molpro_corrections(pending)

        for reac in reactions:
            e_samp = []
            e_high = []
            for sample, job, molp in records[reac]:
                e_stat, energy = molp.get_molpro_energy(
                    key=self.par['vrc_tst_scan_molpro_key'], name=job,
                    VTS=True)
                if not e_stat or not self._vrc_output_matches_input(job):
                    raise RuntimeError(f'VRC Molpro result {job} has no '
                                       f'{self.par["vrc_tst_scan_molpro_key"]} energy.')
                (e_samp if sample else e_high).append(energy)

            dist = ([30] if noscan
                    else self.par['vrc_tst_scan_points'] + [30])
            ens = []  # energies for sample and high in kcal/mol
            asyms = []  # asymptotic energies in hartree
            for energies in (e_samp, e_high):
                if len(energies) != len(dist):
                    raise RuntimeError(f'Incomplete VRC energy series for {reac}.')
                asyms.append(energies[-1])
                ens.append(list((np.array(energies) - energies[-1])
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
                'levels': {
                    'sampling': {
                        'method': self.par['vrc_tst_sample_method'],
                        'basis': self.par['vrc_tst_sample_basis'],
                    },
                    'trusted_correction': {
                        'method': self.par['vrc_tst_high_method'],
                        'basis': self.par['vrc_tst_high_basis'],
                    },
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
        from kinbot.anl.dispatch import _load, preflight, prepare

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
                'resources': {
                    'cores': 'auto', 'memory_mb': 'node',
                    'walltime': self.par['vrc_tst_walltime'],
                    'max_cores': self.par['single_point_ppn'],
                    'min_stack_mw': self.par['vrc_tst_min_stack_mw'],
                    'partition': self.par['queue_name'],
                },
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
        preflight(run_dir)
        result = subprocess.run(
            [sys.executable, '-m', 'kinbot.anl.dispatch', 'drive',
             str(run_dir), '--interval', '20'], check=False)
        if result.returncode:
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
