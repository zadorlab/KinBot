import numpy as np
import time
import logging
import copy

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
                                self.qc.qc_vts_frag(prod)
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
        - If it's a real well, then we use the bond of hom_sci, or the breaking ts bond (only 1 is allowed)
        - If it's a vdw well, then find the closest atoms between the two parts using the core function
        """

        wellparts, wellmaps = self.well.start_multi_molecular()
        if len(wellparts) == 1:
            self.well.find_bond()
            for reac in reactions:
                # find the fragments from the IRCs
                self.scan_reac[reac].parts, self.scan_reac[reac].maps = self.scan_reac[reac].irc_prod.start_multi_molecular()
                self.scan_reac[reac].parts[0].characterize()
                self.scan_reac[reac].parts[1].characterize()
                self.match_order(reac)

                if 'hom_sci' in reac:
                    ww = self.scan_reac[reac].instance_name.split('_')
                    self.scan_reac[reac].scan_coo = [int(ww[-2]) - 1, int(ww[-1]) - 1]
                    logger.info(f'Bond to be scanned for {reac} is {np.array(self.scan_reac[reac].scan_coo)+1}')
                else:
                    # adding up bond breaks only
                    nreacbond = np.sum(np.array([b for b in bb for bb in reac.ts.reac_bond if b < 0])) / 2
                    if nreacbond == 0:
                        logger.warning(f'No bond change detected, unable to determine scan coo for {reac}')
                        self.scan_reac[reac].scan_coo = False
                    elif nreacbond < -1:
                        logger.warning(f'More than one bonds changed in {reac}, unable to determine scan coo for {reac}')
                        self.scan_reac[reac].scan_coo = False
                    else:
                        hit = False
                        for i, bb in enumerate(self.scan_reac[reac].ts.reac_bond):
                            for j, b in bb:
                                if b < 0:
                                    logger.info(f'Bond to be scanned for {reac} is {i+1}-{j+1}')
                                    self.scan_reac[reac].scan_coo = [i, j]
                                    hit = True
                                    break
                            if hit:
                                break

        elif len(wellparts) == 2:  # assuming that there is just a single scan for a vdW, between the closest atoms
            # save the fragments based on the original split
            for reac in reactions:
                self.scan_reac[reac].parts, self.scan_reac[reac].maps = wellparts, wellmaps
                self.scan_reac[reac].scan_coo = well.add_extra_bond()
                self.match_order(reac)
                logger.info(f'Bond to be scanned for {reac} is {self.scan_reac[reac].scan_coo}')
        else:
            self.scan_reac[reac].scan_coo = False
            logging.warning(f'Well {well.chemid} has more than two parts. Cannot do VTS.')

        return

    def match_order(self, reac): 
        '''
        Keeps the object to be scanned intact, but rearranges the products and the atoms in the products so that:
        1. fragment0 of scan is product0
        2. in both fragments the order of the atoms, labeled by chemid, is identical
        It might still scramble chirality, but that is perhaps a very rare edge case?
        '''
        # match product ordering with scanned fragments' order 
        if self.scan_reac[reac].parts[0].chemid != self.scan_reac[reac].products[0].chemid:
            self.scan_reac[reac].products = list(reversed(self.scan_reac[reac].products))
        # match atom ordering of individual products with scanned fragment's atom order
        for fi, frag in enumerate(self.scan_reac[reac].parts):
            for ii, aid in enumerate(frag.atomid):
                if aid != self.scan_reac[reac].products[fi].atomid[ii]: 
                    # swap to something that works
                    pos = self.scan_reac[reac].products[fi].atomid[ii:].index(aid)
                    # need to swap atom and geom
                    self.scan_reac[reac].products[fi].geom[ii],  self.scan_reac[reac].products[fi].geom[pos] = \
                        self.scan_reac[reac].products[fi].geom[pos],  self.scan_reac[reac].products[fi].geom[ii]
                    self.scan_reac[reac].products[fi].atom[ii],  self.scan_reac[reac].products[fi].atom[pos] = \
                        self.scan_reac[reac].products[fi].atom[pos],  self.scan_reac[reac].products[fi].atom[ii]
                    self.scan_reac[reac].products[fi].characterize()
        return

    def do_scan(self, reactions):
        """
        There are two different scans to be made.
        1. relaxed scan: all degrees of freedom are optimized except the B-C distance that is scanned.
           Final geometry is saved in each step.. 
           The user can request that the RMSD < then a threshold during optimization..
           Also, if the optimization crashed, the last valid point is taken.
        2. frozen scan: no degrees of freedom are optimized.
           The frozen fragments are oriented in a way the minimizes RMSD relative to the relaxed scan.
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
                    jobs[ri] = self.qc.qc_vts(self.scan_reac[reac], 
                                              geoms[ri], 
                                              step=step[ri],
                                              )
                    status[ri] = 'running'
                elif status[ri] == 'running':
                    qcst, geom = self.qc.get_qc_geom(jobs[ri], self.scan_reac[reac].species.natom, allow_error=1)
                    if qcst in ['normal', 'error']:
                        status[ri] = 'ready'
                        step[ri] += 1
                        geoms[ri] = copy.deepcopy(geom)
                elif step[ri] == len(self.par['vrc_tst_scan_points']):
                    status[ri] = 'done'
            if len([st for st in status if st == 'done']) == len(reactions):
                break
            else:
                time.sleep(1)

        return
