import numpy as np
import os
import logging
from dataclasses import dataclass
from shutil import copyfile

from kinbot.utils import plain_symbols, plain_geometry
from kinbot import kb_path
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_identity import require_supported_identity, identity_matches

logger = logging.getLogger('KinBot')


@dataclass
class BimolecularIRCResult:
    """Optimized endpoints connected by a bimolecular transition state."""

    reactant: StationaryPoint
    product: StationaryPoint
    reactant_direction: str
    product_direction: str
    product_matches_expected: bool


class IRC:
    """
    Class to run the IRC's for one specific reaction
    """
    def __init__(self, rxn, par):
        # instance of the reac_family object
        # the family this reaction belongs to
        self.rxn = rxn
        self.par = par

    @property
    def qc(self):
        # A final bimolecular IRC validates the selected TS on its own surface.
        if self.par['bimol'] and getattr(self.rxn, 'final_irc_name', None):
            source = getattr(getattr(self.rxn, 'ts', None), 'source_job', None)
            if source:
                return self.rxn.qc.for_job(source)
            level = 'l2' if self.par['high_level'] else 'l1'
            return self.rxn.qc.for_level(level)
        return self.rxn.qc.for_level('l1')

    def _job_name(self):
        if self.par['bimol']:
            return getattr(self.rxn, 'final_irc_name', self.rxn.instance_name)
        return self.rxn.instance_name

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
        instance_name = self._job_name()

        directions = ['Forward', 'Reverse']

        ini_well_hits = 0
        ini_well_hit = None
        prod_hit = None
        product_matches_expected = False
        st_pts = [None, None]
        reference_identity = require_supported_identity(self.rxn.species)
        for i, direction in enumerate(directions):
            irc_name = '{}_IRC_{}_prod'.format(instance_name, direction[0])
            err, geom = self.qc.get_qc_geom(irc_name,
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
                from kinbot.bimol import _expected_graphs

                search = self.rxn.get_bimol_search_data()
                reactant_graph, product_graph = _expected_graphs(
                    self.rxn.species, search['formed_bonds'],
                    search['broken_bonds'])
                graph = {
                    (first, second)
                    for first in range(temp.natom)
                    for second in range(first + 1, temp.natom)
                    if temp.bond01[first][second]
                }
                if graph == reactant_graph:
                    if not identity_matches(
                            reference_identity, require_supported_identity(temp),
                            self.par.get('optical_population', 'specified')):
                        logger.info(
                            '\tIRC reactant configuration differs from the requested '
                            'population for %s.', instance_name)
                        return 0
                    ini_well_hits += 1
                    ini_well_hit = i
                else:
                    prod_hit = i
                    product_matches_expected = graph == product_graph
            else:
                matches = identity_matches(reference_identity, require_supported_identity(temp),
                                           self.par.get('optical_population', 'specified'))
                if matches:
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
            if self.par['bimol']:
                if not product_matches_expected:
                    logger.info(
                        '\tIRC for %s found an unexpected product; retaining '
                        'the saddle and characterizing the observed pathway.',
                        instance_name)
                return BimolecularIRCResult(
                    reactant=st_pts[ini_well_hit],
                    product=st_pts[prod_hit],
                    reactant_direction=directions[ini_well_hit][0],
                    product_direction=directions[prod_hit][0],
                    product_matches_expected=product_matches_expected,
                )
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
        instance_name = self._job_name()
        directions = ['Forward', 'Reverse']
        status = [-1, -1]
        for i, direction in enumerate(directions):
            irc_name = '{}_IRC_{}'.format(instance_name, direction[0])
            status[i] = self.qc.check_qc(irc_name)
            if (self.par['bimol'] and status[i] == 'normal'
                    and getattr(self.qc, 'qc', None) == 'gauss'
                    and not getattr(self.qc, 'use_sella', False)):
                # Old workers could mark failed corrector trials as normal.
                # Require the same native path evidence on saved jobs too.
                from kinbot import reader_gauss
                try:
                    _, accepted, _, _ = reader_gauss.read_all_irc_geoms(irc_name + '.log')
                except OSError:
                    accepted = None
                if accepted is None:
                    logger.warning('%s: no accepted Gaussian IRC points; saved endpoint '
                                   'does not validate this path.', irc_name)
                    status[i] = 'error'
        return status

    def do_irc_calculations(self):
        """
        Carry out the IRC calculation.
        """
        instance_name = self._job_name()
        if self.par['bimol'] and getattr(self.rxn, 'ts', None) is not None:
            geom = np.asarray(self.rxn.ts.geom)
        else:
            err, geom = self.qc.get_qc_geom(
                self.rxn.instance_name, self.rxn.species.natom)
        directions = ['Forward', 'Reverse']
        for direction in directions:
            irc_name = '{}_IRC_{}'.format(instance_name, direction[0])
            # This boolean is false if the checkpoint file is available
            # and true if no checkpoint file is found.
            # In the latter case, the geometry needs to be supplied to
            # the gaussian calculation and the keywords
            # geom(AllCheck,NoKeepConstants) guess=Read need to be removed
            start_from_geometry = 0
            if self.qc.qc == 'gauss':
                code = 'gaussian'  # Sella
                Code = 'Gaussian'  # Sella
                # copy the chk file
                if self.par['bimol'] and getattr(self.rxn, 'final_irc_name', None):
                    start_from_geometry = 1
                elif os.path.exists(instance_name + '.chk'):
                    copyfile(instance_name + '.chk', irc_name + '.chk')
                else:
                    start_from_geometry = 1
            elif self.qc.qc == 'qchem':
                code = 'qchem'  # Sella
                Code = 'QChem'  # Sella
            elif self.qc.qc == 'orca':
                code = 'orca'
                Code = 'ORCA'
            elif self.qc.qc == 'fc':
                code = 'fairchem'
                Code = 'Fairchem'
            else:
                raise ValueError(f'Unexpected code name: {self.qc.qc}.')

            kwargs = self.qc.get_qc_arguments(irc_name,
                                                  self.rxn.species.mult,
                                                  self.rxn.species.charge,
                                                  self.rxn.species.nel,
                                                  irc=direction.lower(),
                                                  start_from_geom=start_from_geometry)
            if self.par['bimol'] and start_from_geometry and self.qc.qc == 'gauss':
                # The final bimolecular IRC uses the selected TS coordinates,
                # not the discovery checkpoint and its force constants.
                kwargs['irc'] = ','.join(option for option in kwargs['irc'].split(',')
                                         if option.upper() != 'RCFC')
            prod_kwargs = self.qc.get_qc_arguments(irc_name + '_prod', self.rxn.species.mult, self.rxn.species.charge, self.rxn.species.nel)
            if self.qc.qc == 'gauss':
                #prod_kwargs['opt'] = 'CalcFC, Tight'
                prod_kwargs['opt'] = 'CalcFC'
            if self.qc.par['calc_kwargs']:
                kwargs = self.qc.merge_kwargs(kwargs)
                prod_kwargs = self.qc.merge_kwargs(prod_kwargs)
            if self.qc.use_sella:
                kwargs.pop('irc', None)
                kwargs.pop('geom', None)
                kwargs.pop('guess', None)
                prod_kwargs.pop('opt', None)
                prod_kwargs.pop('freq', None)
                if self.par['bimol']:
                    template_file = (
                        f'{kb_path}/tpl/ase_sella_bimol_irc.tpl.py')
                else:
                    template_file = f'{kb_path}/tpl/ase_sella_irc.tpl.py'
            else:
                template_file = f'{kb_path}/tpl/ase_{self.qc.qc}_irc.tpl.py'
            template = open(template_file, 'r').read()
            template = template.format(label=irc_name,
                                       frequency_mode=self.qc.frequency_mode,
                                       optimizer=self.qc.optimizer,
                                       profiled=self.qc.theory_level is not None,
                                       kwargs=kwargs,
                                       prod_kwargs=prod_kwargs,
                                       atom=plain_symbols(self.rxn.species.atom),
                                       geom=plain_geometry(geom),
                                       ppn=min(self.rxn.species.nel, self.qc.ppn),
                                       qc_command=self.qc.par['qc_command'],
                                       working_dir=os.getcwd(),
                                       code=code,
                                       Code=Code,
                                       sella_kwargs=self.qc.par['sella_kwargs'],  # Sella
                                       fmax=self.qc.par['sella_fmax'],
                                       steps=self.qc.par['sella_steps'],
                                       irc_maxpoints=(self.par['irc_maxpoints']
                                                      if self.par['bimol'] else 100),
                                       fc_model_path=self.qc.par['fc_model_path'],
                                       fc_task_name=self.qc.par['fc_task_name'],
                                       fc_device=self.qc.par['fc_device'],
                                       fc_compile=self.qc.par.get('fc_compile', False),
                                       bimol=bool(self.par['bimol']),
                                       )

            with open('{}.py'.format(irc_name), 'w') as f:
                f.write(template)

            self.qc.submit_qc(irc_name, min(self.rxn.species.nel, self.qc.ppn), singlejob=0)

        return 0
