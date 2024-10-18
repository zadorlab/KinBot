from operator import ne
from typing import Any
import numpy as np
import time
import logging
import copy
import os
import stat
import json

from ase.db import connect
from shutil import which
from subprocess import Popen, PIPE
from kinbot.utils import reorder_coord
from kinbot.stationary_pt import StationaryPoint
from kinbot import geometry
from kinbot.molpro import Molpro
from kinbot.utils import queue_command, create_matplotlib_graph, NpEncoder
from kinbot import constants

logger = logging.getLogger('KinBot')


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
        # structure of par['vrc_tst_scan']:
        # {chemid1: ["reaction_name1","reaction_name2], chemid2: [...]}
        for chemid, reactions in self.par['vrc_tst_scan'].items():
            # select the part of the input relevant for this kb
            if chemid == str(self.well.chemid):
                self.opt_products(reactions)
                self.save_products(reactions)
                self.find_scan_coos(reactions)
                jobs = self.do_scan(reactions)
                self.energies(reactions)
        return

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
                        for prod in self.scan_reac[reac].products:
                            if prod.chemid in prod_chemid:
                                continue
                            else:  # submit at vrc_tst level
                                prod_chemid.append(prod.chemid)
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
                job = f'vrctst/{str(prod.chemid)}_vts'
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
            for ai in self.par['vrc_tst_scan_reac_cent'][str(prod.chemid)]:
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
        if self.scan_reac[reac].parts[0].chemid !=\
           self.scan_reac[reac].products[0].chemid:
            self.scan_reac[reac].products = \
                list(reversed(self.scan_reac[reac].products))
        if self.scan_reac[reac].parts[0].chemid !=\
           self.scan_reac[reac].products[0].chemid:
            logger.warning('IRC prod and prod are not the same!')
        if self.scan_reac[reac].parts[1].chemid !=\
           self.scan_reac[reac].products[1].chemid:
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

    def do_scan(self, reactions):
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
            # Stores index of radicals and their equivalent
            self.scan_reac[reac].usym = [[], []]
            if not self.scan_reac[reac].do_vdW:
                geoms.append(self.scan_reac[reac].species.geom)
            else:
                geoms.append(self.scan_reac[reac].irc_prod.geom)
        jobs = [''] * len(reactions)
        equiv = []
        step0_geoms = [np.array(geoms[ri]) for ri in range(len(reactions))]

        while 1:
            for ri, reac in enumerate(reactions):
                if status[ri] == 'ready' and step[ri] < \
                   len(self.par['vrc_tst_scan_points']) + 1:
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
                    if step < len(self.par['vrc_tst_scan_points']):
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
                    if step[ri] == 0:
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
                        equiv.append([equiv_A, equiv_B])
                        self.scan_reac[reac].equiv = [equiv_A, equiv_B]

                    jobs[ri] = self.qc.qc_vts(self.scan_reac[reac],
                                              geoms[ri],
                                              step[ri],
                                              equiv[ri],
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
                elif step[ri] == len(self.par['vrc_tst_scan_points']) + 1:
                    status[ri] = 'done'
            if len([st for st in status if st == 'done']) == len(reactions):
                break
            else:
                time.sleep(1)

        return (jobs)

    def energies(self, reactions):
        '''
        Create and submit molpro calculations
        '''
        cmd, ext = queue_command(self.par['queuing'])
        batch_submit = ''
        db = connect('kinbot.db')
        for reac in reactions:
            all_done = True
            e_samp = []
            e_high = []
            for step in range(len(self.par['vrc_tst_scan_points']) + 1):
                for sample in [True, False]:
                    if sample:
                        if step < len(self.par['vrc_tst_scan_points']):
                            job = f'{reac}_vts_pt{str(step).zfill(2)}_fr'
                        else:
                            job = f'{reac}_vts_pt_asymptote_fr'
                    else:
                        if step < len(self.par['vrc_tst_scan_points']):
                            job = f'{reac}_vts_pt{str(step).zfill(2)}'
                        else:
                            # take the geometry from the _fr case as well here
                            job = f'{reac}_vts_pt_asymptote_fr'
                    *_, last_row = db.select(name=f'vrctst/{job}', sort='-1')
                    scan_spec = StationaryPoint.from_ase_atoms(
                        last_row.toatoms())
                    scan_spec.characterize()

                    molp = Molpro(scan_spec, self.par)
                    if not sample and step == len(
                       self.par['vrc_tst_scan_points']):
                        job = job[:-3]  # the actual job name to run
                    molp.create_molpro_input(name=job, VTS=True, sample=sample)
                    molp.create_molpro_submit(name=job, VTS=True)
                    e_stat, e = \
                        molp.get_molpro_energy(
                            key=self.par['vrc_tst_scan_molpro_key'],
                            name=f'{job}',
                            VTS=True)
                    logger.debug(f'{job}, {e_stat}, {e}')
                    if not e_stat:
                        batch_submit += f'{cmd} {job}.{ext}\n'
                        all_done = False
                    elif sample:
                        e_samp.append(e)
                    else:
                        e_high.append(e)

            if all_done:
                dist = self.par['vrc_tst_scan_points'] + [30]
                ens = []  # energies for sample and high in kcal/mol
                asyms = []  # asymptotic energies in hartree
                for sample in [True, False]:
                    if sample:
                        eee = copy.copy(e_samp)
                    else:
                        eee = copy.copy(e_high)
                    asyms.append(eee[-1])
                    eee = list((np.array(eee) - eee[-1]) * constants.AUtoKCAL)
                    ens.append(eee)

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

                # TODO instead of writing files, create and save png
                # TODO simple text file with 3 columns: R, e_samp, e_high
                create_matplotlib_graph(x=dist,
                                        data=ens,
                                        name=f'{reac}',
                                        x_label=f"{reac}",
                                        y_label="Energy (kcal/mol)",
                                        data_legends=['sample', 'high'],
                                        )

                smallest = np.linalg.norm(
                    self.well.geom[self.scan_reac[reac].equiv[0][0]] -
                    self.well.geom[self.scan_reac[reac].equiv[1][0]])

                # write small file with correction data
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
                    'frags_atom': [self.scan_reac[reac].products[0].atom,
                                   self.scan_reac[reac].products[1].atom],
                    'frags_geom': [self.scan_reac[reac].products[0].geom,
                                   self.scan_reac[reac].products[1].geom],
                    'frags_mult': [self.scan_reac[reac].products[0].mult,
                                   self.scan_reac[reac].products[1].mult]
                    }

                with open(f'vrctst/corr_{reac}.json',
                          'w',
                          encoding='utf-8') as f:
                    json.dump(corr,
                              f,
                              ensure_ascii=False,
                              indent=4,
                              cls=NpEncoder)

        batch = f'vrctst/molpro/batch_vts_{self.par["queuing"]}.sub'
        if self.par['queuing'] != 'local' and batch_submit != '':
            with open(batch, 'w') as f:
                f.write(batch_submit)
            os.chmod(batch, stat.S_IRWXU)  # read, write, execute by owner
        return
