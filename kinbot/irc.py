import numpy as np
import os
import logging
from shutil import copyfile

from kinbot import kb_path
from kinbot.stationary_pt import StationaryPoint

logger = logging.getLogger('KinBot')


class IRC:
    """
    Class to run the IRC's for one specific reaction
    """
    def __init__(self, rxn, par):
        # instance of the reac_family object
        # the family this reaction belongs to
        self.rxn = rxn
        self.par = par

    def irc2stationary_pt(self):
        """
        Read the irc files
        There are three possible scenarios:
        1. One of the ircs leads the initial well and
           the other to another well or bimolecular product
        2. Neither of the ircs lead to the inital well,
           transition state structure is not the one
           kinbot was looking for
        3. Both the ircs lead to the initial well,
           KinBot found either an identical reaction
           or the ts is not correct
        """
        instance_name = self.rxn.instance_name

        directions = ['Forward', 'Reverse']

        ini_well_hits = 0
        prod_hit = None
        st_pts = [None, None]
        for i, direction in enumerate(directions):
            irc_name = '{}_IRC_{}_prod'.format(instance_name, direction[0])
            err, geom = self.rxn.qc.get_qc_geom(irc_name,
                                                self.rxn.species.natom,
                                                allow_error=1)
            if err == -1:
                return 0
            if self.problem_in_geom(geom):
                # this happens seldom that all the atoms are
                # very close to one another (problem in Gaussian)
                logger.warning('Problem with product geometry for {}'.format(instance_name))
                return 0

            temp = StationaryPoint(irc_name,
                                   self.rxn.species.charge,
                                   self.rxn.species.mult,
                                   atom=self.rxn.species.atom,
                                   geom=geom)
            temp.characterize()

            st_pts[i] = temp
            if self.par['bimol']:
                fragments, maps = temp.start_multi_molecular()
                if len(fragments) == 1 or len(fragments) > 2:
                    prod_hit = i
                elif len(fragments) == 2:
                    fragments[0].characterize()
                    fragments[1].characterize()
                    if (fragments[0].chemid == self.rxn.species.fragA.chemid and fragments[1].chemid == self.rxn.species.fragB.chemid) or \
                       (fragments[1].chemid == self.rxn.species.fragA.chemid and fragments[0].chemid == self.rxn.species.fragB.chemid):
                        ini_well_hits += 1
                    else:
                        prod_hit = i
            else:
                if temp.chemid == self.rxn.species.chemid and all(temp.chiral[at] == self.rxn.species.chiral[at] for at in range(self.rxn.species.natom)):
                    ini_well_hits += 1
                else:
                    prod_hit = i  # this leaves the possibility of a chirality changing reaction

        if ini_well_hits == 0:
            if self.par['bimol']:
                logger.info('\tNeither IRC leads to the initial reactants for {}'.format(instance_name))
            else:
                logger.info('\tNeither IRC leads to the initial well for {}'.format(instance_name))
            return 0
        elif ini_well_hits == 2:
            if self.par['bimol']:
                logger.info('\tBoth IRCs lead to the initial reactants, identical reaction found: {}'.format(instance_name))
            else:
                logger.info('\tBoth IRCs lead to the initial well, identical reaction found: {}'.format(instance_name))
            return 0
        else:
            # ircs OK: well and product found
            logger.info('\tIRCs successful for {}'.format(instance_name))
            return st_pts[prod_hit]

    def problem_in_geom(self, geom):
        # check if interatomic distances are closer than 0.3 Angstrom
        for i in range(len(geom)):
            for j in range(i + 1, len(geom)):
                dist = np.linalg.norm(geom[i] - geom[j])
                if dist < 0.3:
                    return 1
        return 0

    def check_irc(self):
        instance_name = self.rxn.instance_name
        directions = ['Forward', 'Reverse']
        status = [-1, -1]
        for i, direction in enumerate(directions):
            irc_name = '{}_IRC_{}'.format(instance_name, direction[0])
            status[i] = self.rxn.qc.check_qc(irc_name)
        return status

    def do_irc_calculations(self):
        """
        Carry out the IRC calculation.
        """
        instance_name = self.rxn.instance_name
        err, geom = self.rxn.qc.get_qc_geom(instance_name,
                                            self.rxn.species.natom)
        directions = ['Forward', 'Reverse']
        for direction in directions:
            irc_name = '{}_IRC_{}'.format(instance_name, direction[0])

            # This boolean is false if the checkpoint file is available
            # and true if no checkpoint file is found.
            # In the latter case, the geometry needs to be supplied to
            # the gaussian calculation and the keywords
            # geom(AllCheck,NoKeepConstants) guess=Read need to be removed
            start_from_geometry = 0
            if self.rxn.qc.qc == 'gauss':
                code = 'gaussian'  # Sella
                Code = 'Gaussian'  # Sella
                # copy the chk file
                if os.path.exists(instance_name + '.chk'):
                    copyfile(instance_name + '.chk', irc_name + '.chk')
                else:
                    start_from_geometry = 1
            elif self.rxn.qc.qc == 'nwchem':
                code = 'nwchem'  # Sella
                Code = 'NWChem'  # Sella
                if direction == 'Reverse':
                    direction = 'Backward'
            elif self.rxn.qc.qc == 'qchem':
                code = 'qchem'  # Sella
                Code = 'QChem'  # Sella
            elif self.rxn.qc.qc == 'nn_pes':
                code = 'nn_pes'
                Code = 'Nn_surr'
            else:
                raise ValueError(f'Unexpected code name: {self.rxn.qc.qc}.')

            kwargs = self.rxn.qc.get_qc_arguments(irc_name,
                                                  self.rxn.species.mult,
                                                  self.rxn.species.charge,
                                                  irc=direction.lower(),
                                                  start_from_geom=start_from_geometry)
            prod_kwargs = self.rxn.qc.get_qc_arguments(irc_name + '_prod', self.rxn.species.mult, self.rxn.species.charge)
            if self.rxn.qc.qc == 'gauss':
                #prod_kwargs['opt'] = 'CalcFC, Tight'
                prod_kwargs['opt'] = 'CalcFC'
            if self.rxn.par['calc_kwargs']:
                kwargs = self.rxn.qc.merge_kwargs(kwargs)
                prod_kwargs = self.rxn.qc.merge_kwargs(prod_kwargs)
            if self.rxn.qc.use_sella:
                kwargs.pop('irc', None)
                kwargs.pop('geom', None)
                kwargs.pop('guess', None)
                prod_kwargs.pop('opt', None)
                prod_kwargs.pop('freq', None)
                template_file = f'{kb_path}/tpl/ase_sella_irc.tpl.py'
            else:
                template_file = f'{kb_path}/tpl/ase_{self.rxn.qc.qc}_irc.tpl.py'
            template = open(template_file, 'r').read()
            template = template.format(label=irc_name,
                                       kwargs=kwargs,
                                       prod_kwargs=prod_kwargs,
                                       atom=list(self.rxn.species.atom),
                                       geom=list([list(gi) for gi in geom]),
                                       ppn=self.rxn.qc.ppn,
                                       qc_command=self.par['qc_command'],
                                       working_dir=os.getcwd(),
                                       code=code,
                                       Code=Code,
                                       sella_kwargs=self.par['sella_kwargs']  # Sella
            )

            with open('{}.py'.format(irc_name), 'w') as f:
                f.write(template)

            self.rxn.qc.submit_qc(irc_name, singlejob=0)

        return 0
