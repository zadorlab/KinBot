# -*- coding: utf-8 -*-
import copy
import logging
from kinbot import constants
from kinbot.optimize import Optimize
from kinbot.stationary_pt import StationaryPoint


class HomolyticScission:
    """
    Class with the information of one homolytic scission reaction
    """
    def __init__(self, species, par, qc, bond):
        self.species = species
        self.par = par
        self.qc = qc
        self.bond = bond  # bond that is being broken
        self.products = []  # list of stationary points with the products
        self.prod_opt = []  # list of optimize object for each product
        self.status = 0

    def create_geometries(self):
        """
        Cut the species in two parts and save the geometries of both parts into
        new stationary point objects
        """
        # copy the reactant stationary point
        atom = copy.deepcopy(self.species.atom)
        geom = copy.deepcopy(self.species.geom)
        temp = StationaryPoint('temp', self.species.charge,
                               self.species.mult, atom=atom, geom=geom)
        temp.characterize()
        # set the bond order of the breaking bond to 0
        temp.bond[self.bond[0]][self.bond[1]] = 0
        temp.bond[self.bond[1]][self.bond[0]] = 0
        self.products, maps = temp.start_multi_molecular()


class HomolyticScissions:
    """
    Class to find all the potential homolytic scission reactions
    """
    def __init__(self, species, par, qc):
        self.species = species
        self.qc = qc
        self.par = par
        # list of homolytic scission objects
        self.hss = []

    def find_homolytic_scissions(self):
        """
        Enumerate all unique homolytic scission reactions
        """
        # set of unique bonds to break
        bonds = []
        for i in range(len(self.species.atom)-1):
            for j in range(i+1, len(self.species.atom)):
                # only consider a bond of which both
                # atoms are not in the same cycle
                cycle = []
                for cyc in self.species.cycle_chain:
                    if i in cyc:
                        cycle.extend(cyc)
                if j not in cycle:
                    # only consider single bonds for homolytic scissions
                    if self.species.bond[i][j] == 1:
                        # check if the bond is present in the 'homolytic_bonds' list, or if the list is empty
                        # add is a boolean that tells if the bond needs to be added to the homolytic scissions
                        # to be considered (i.e. the 'bonds' list)
                        add = 0
                        # get the elements of the current bond under consideration, and put them
                        # in alphabetical order
                        bond_elements = sorted([self.species.atom[i], self.species.atom[j]])
                        # check all bonds to be considered, supplied by the user
                        for user_bond in self.par['homolytic_bonds']:
                            if sorted(user_bond) == bond_elements:
                                # if there is a match between the current bond and a user-defined
                                # bond, put the add boolean to 1
                                add = 1
                        # add the bond if the 'add' boolean is 1 or if the user defined list is empty
                        # (in the latter case, all single bonds that are not in a cycle are considered. 
                        if len(self.par['homolytic_bonds']) == 0 or add:
                            # check if a bond with identical atomids
                            # has been added to the bonds list yet
                            new = 1
                            for bi in bonds:
                                if sorted([self.species.atomid[at] for at in bi]) == sorted([self.species.atomid[i], self.species.atomid[j]]):
                                    new = 0
                            if new:
                                bonds.append([i, j])
                                hs = HomolyticScission(self.species, self.par,
                                                       self.qc, [i, j])
                                hs.create_geometries()
                                self.hss.append(hs)
        
        # optimize the products of the hss
        while not all([hs.status < 0 for hs in self.hss]):
            for index, hs in enumerate(self.hss):
                if hs.status == 0:
                    # do the initial optimization
                    for prod in hs.products:
                        hs.qc.qc_opt(prod, prod.geom)
                    hs.status = 1
                if hs.status == 1:
                    # wait for the optimization to finish
                    err = 0
                    for prod in hs.products:
                        e, prod.geom = hs.qc.get_qc_geom(str(prod.chemid) + '_well', prod.natom)
                        if e < 0:
                            # optimizatin failed
                            logging.info("HS optimization failed for {}".format(prod.chemid))
                            hs.status = -999
                            err = -1
                        elif e != 0:
                            err = -1
                        else:
                            e2, prod.energy = hs.qc.get_qc_energy(str(prod.chemid) + '_well')
                            e2, prod.zpe = hs.qc.get_qc_zpe(str(prod.chemid) + '_well')
                    if err == 0:
                        hs.status = 2
                if hs.status == 2:
                    # Do the product conf search, high level opt and HIR
                    err = 0
                    for i, prod in enumerate(hs.products):
                        chemid = prod.chemid
                        prod_opt = Optimize(prod, self.par, self.qc)
                        if str(chemid) != str(prod_opt.species.chemid):
                            logging.info("HS product {} changed to {} during optimization.".format(chemid, prod_opt.species.chemid))
                            hs.qc.qc_opt(prod_opt.species, prod_opt.species.geom)
                            j = i - 1
                            # wait for new opt to finish
                            er, prod.geom = hs.qc.get_qc_geom(str(prod_opt.species.chemid) + '_well', prod_opt.species.natom)
                            if er < 0:
                                # optimization failed
                                logging.info("HS optimization failed for {}".format(prod_opt.species.chemid))
                                err = -1
                            elif er != 0:
                                err = -1
                            else:
                                er2, prod.energy = hs.qc.get_qc_energy(str(prod_opt.species.chemid) + '_well')
                                er2, prod.zpe = hs.qc.get_qc_zpe(str(prod_opt.species.chemid) + '_well')
                            hs.products.pop(i)
                            hs.products.insert(j, prod_opt.species)
                        prod_opt.do_optimization()
                        hs.prod_opt.append(prod_opt)
                    if err == 0:
                        hs.status = 3
                if hs.status == 3:
                    # check up on the optimization
                    opts_done = 1
                    fails = 0
                    for pr_opt in hs.prod_opt:
                        if not pr_opt.shir == 1:
                            opts_done = 0
                            pr_opt.do_optimization()
                        if pr_opt.shigh == -999:
                            fails = 1
                    if fails:
                        hs.status = -999
                    elif opts_done:
                        # check if the energy is higher
                        # than the barrier threshold
                        species_zeroenergy = self.species.energy + self.species.zpe
                        prod_zeroenergy = 0.
                        for pr_opt in hs.prod_opt:
                            prod_zeroenergy += pr_opt.species.energy + pr_opt.species.zpe
                        barrier = (prod_zeroenergy - species_zeroenergy) * constants.AUtoKCAL
                        prod_name = ' '.join(sorted([str(prod.species.chemid) for prod in hs.prod_opt]))
                        if barrier > self.par['barrier_threshold']:
                            logging.info("Energy of HS product {} is above the barrier threshold ({:.3} kcal/mol)".format(prod_name, barrier))
                            hs.status = -999
                        else:
                            hs.status = -1
                            name = '_'.join(sorted([str(prod.species.chemid) for prod in hs.prod_opt]))
                            logging.info('Homolytic scission (asymptote {:.2f} kcal/mol) lead to products {}'.format(barrier, name))
