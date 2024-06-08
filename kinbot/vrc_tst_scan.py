import numpy as np
import time
import os
import logging
from shutil import copyfile

from kinbot import kb_path
from kinbot.stationary_pt import StationaryPoint

logger = logging.getLogger('KinBot')


class VTS:
    """
    Class to run the scans for VRC-TST to eventually generate correction potentials
    """
    def __init__(self, well, par, qc):
        # instance of the VTS object
        self.well = well
        self.par = par
        self.qc = qc

        self.scan_reac = {}  # reaction object for reactions scanned 
        self.scan_coo = {}
        self.prod_chemids = []

    def calculate_correction_potentials(self):
        """
        The main driver for scanning the potential.
        """
        # structure of par['vrc_tst_scan']: {chemid1: ["reaction_name1","reaction_name2], chemid2: [...]} 
        for chemid, reactions in self.par['vrc_tst_scan'].items():
            if chemid == str(self.well.chemid):  # select the part of the input relevant for this kb
                self.opt_products(reactions)
                self.save_products(reactions)
                self.find_scan_coos(reactions)
                self.do_scan(reactions)
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
                        logger.warning(f'There are {len(ro.products)} products for {reac}, cannot do scan if not bimolecular.')
                    else:
                        self.scan_reac[reac] = ro  # everything is there about the original reactions 
                        for prod in self.scan_reac[reac].products:
                            if prod.chemid in prod_chemid:
                                continue
                            else:  # submit at vrc_tst level, we have the geometry and everything else, mult etc.
                                prod_chemid.append(prod.chemid)    
                                self.qc.qc_vts(prod, prod.geom)
            if not hit:
                logger.warning(f'Reaction {reac} requested to scan is not found in reaction list.')

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
                status, geom = self.qc.get_qc_geom(job, prod.natom, wait=1, allow_error=1)
                # update geom in self.scan_reac[reac] values, this already saves it
                prod.geom = geom

        return

    def find_scan_coos(self, reactions):
        """
        Find bond to be scanned
        This is done in two ways
        - If it's a real well, then we use the bond of 
        - If it's a vdw well, then find the closest atoms between the two parts using the core function
        """

        parts, maps = self.well.start_multi_molecular()
        if len(parts) == 1:
            self.well.find_bond()
            for reac in reactions:
                if 'hom_sci' in reac:
                    ww = self.scan_reac[reac].instance_name.split('_')
                    self.scan_coo[reac] = [int(ww[-2]), int(ww[-1])]
                    logger.info(f'Bond to be scanned for {reac} is {self.scan_coo[reac]}')
                else:
                    # adding up bond breaks only
                    nreacbond = np.sum(np.array([b for b in bb for bb in reac.ts.reac_bond if b < 0])) / 2
                    if nreacbond == 0:
                        logger.warning(f'No bond change detected, unable to determine scan coo for {reac}')
                        self.scan_coo[reac] = False
                    elif nreacbond < -1:
                        logger.warning(f'More than one bonds changed in {reac}, unable to determine scan coo for {reac}')
                        self.scan_coo[reac] = False
                    else:
                        hit = False
                        for i, bb in enumerate(self.scan_reac[reac].ts.reac_bond):
                            for j, b in bb:
                                if b < 0:
                                    logger.info(f'Bond to be scanned for {reac} is {i+1}-{j+1}')
                                    self.scan_coo[reac] = [i, j]
                                    hit = True
                                    break
                            if hit:
                                break

        elif len(parts) == 2:  # assuming that there is just a single scan for a vdW, between the closest atoms
            self.scan_coo[reac] = well.add_extra_bond()
            logger.info(f'Bond to be scanned for {reac} is {self.scan_coo[reac]}')
        else:
            self.scan_coo[reac] = False
            logging.warning(f'Well {well.chemid} has more than two parts. Cannot do VTS.')

        return

    def do_scan(self, reactions):
        """
        There are two different scans to be made.
        1. relaxed scan: all degrees of freedom are optimized except the B-C distance that is scanned.
            In these scans the interfragmental degrees of freedom are saved. These are:
            - distance B-C
            - angles A-B-C and B-C-D
            - dihedrals A-B-C-D, X-A-B-C (if exists) and B-C-D-Y (if exists)
            A, D, X, and Y are selected as neighbors of B, C, A and D, respectively
            The user can request that certain angles and bond lengths do not change more than a threshold between
            consecutive steps to prevent sudden changes. In this case the last geometry obeying these constraints is picked
            and the corresponding energy. Also, if the optimization crashed, the last valid point is taken.
        2. rigid scan: no degrees of freedom are optimized. The cartesian definition is converted into
            a z-matrix one, elements are adjusted based on saved values of internal coordinates of separated fragmens,
            converted back to cartesian, and a single pt calculation is done.
        Following this, a high-level calculation is done on both set of geometries, just single pt.
        All of these calculations are invoked from a single template per point, and the points are run in sequence.
        """

        status = ['ready'] * len(reactions)
        step = np.zeros(len(reactions))
        geoms = [self.scan_reac[reac].species.geom for reac in reactions]
        jobs = [''] * len(reactions)

        while 1:
            for ri, reac in enumerate(reactions):
                if status[ri] == 'ready' and step[ri] < len(self.par['vrc_tst_scan_points']):
                    jobs[ri] = self.qc.qc_vts(self.scan_reac[reac].species, geoms[ri], scan_coo=self.scan_coo[reac], step=step[ri])
                    status[ri] = 'running'
                elif status[ri] == 'running':
                    qcst, geom = self.qc.get_qc_geom(jobs[ri], self.scan_reac[reac].species.natom, allow_error=1)
                    if qcst in ['normal', 'error']:
                        status[ri] = 'ready'
                        step[ri] += 1
                elif step[ri] == len(self.par['vrc_tst_scan_points']):
                    status[ri] = 'done'
            if len([st for st in status if st == 'done']) == len(reactions):
                break
            else:
                time.sleep(1)

        return
