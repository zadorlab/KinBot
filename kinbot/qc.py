import os
import sys
import subprocess
import logging
import re
import time
from datetime import datetime
import copy

import numpy as np
from ase.db import connect

from kinbot import kb_path
from kinbot import constants
from kinbot import geometry
from kinbot import exceptions

logger = logging.getLogger('KinBot')


class QuantumChemistry:
    '''
    This class provides the link between KinBot and the qc code
    It choses the write options, submits the jobs and checks
    the jobs for success or failure
    '''

    def __init__(self, par):
        self.par = par
        self.qc = par['qc'].lower()
        self.method = par['method']
        self.basis = par['basis']
        self.scan_method = par['scan_method']
        self.scan_basis = par['scan_basis']
        self.bls_method = par['barrierless_saddle_method']
        self.bls_basis = par['barrierless_saddle_basis']
        self.high_level_method = par['high_level_method']
        self.high_level_basis = par['high_level_basis']
        self.bls_high_level_method = par['barrierless_saddle_method_high']
        self.bls_high_level_basis = par['barrierless_saddle_basis_high']
        self.vts_method = par['vrc_tst_scan_method']
        self.vts_basis = par['vrc_tst_scan_basis']
        self.single_point_template = par['single_point_template']
        self.single_point_key = par['single_point_key']
        self.integral = par['integral']
        self.opt = par['opt']
        self.ppn = par['ppn']
        self.queuing = par['queuing']
        self.queue_name = par['queue_name']
        self.slurm_feature = par['slurm_feature']
        self.zf = par['zf']
        self.db = connect('kinbot.db')
        self.job_ids = {}
        self.irc_maxpoints = par['irc_maxpoints']
        self.irc_stepsize = par['irc_stepsize']
        self.qc_command = par['qc_command']
        if par['slurm_feature'] == '':
            self.slurm_feature = ''
        else:
            self.slurm_feature = '#SBATCH -C ' + par['slurm_feature']
        self.queue_job_limit = par['queue_job_limit']
        self.username = par['username']
        self.use_sella = par['use_sella']
        if not self.use_sella and self.qc.lower() == 'nn_pes':
            logger.warning('NNPES needs Sella optimizer. Turning "use_sella" on.')
            self.use_sella = True

    def get_qc_arguments(self, job, mult, charge, ts=0, step=0, max_step=0,
                         irc=None, scan=0, high_level=0, hir=0,
                         start_from_geom=0, rigid=0, aie=0, L3=False, vts=0):
        '''
        Method to get the argument to pass to ase, which are then passed to the qc codes.
        Job: name of the job
        mult: multiplicity of the job
        charge: charge of the job
        ts: 1 for transition state searches, 0 for wells and products
        step: Step number in the transition state search
        max_step: total number of necessary steps to get to the final transition state structure
            This is different for every reaction family
        irc: direction of the irc, None if this is not an irc job
        scan: this calculation is part of a scan of a bond length
        '''

        if self.qc == 'gauss' or (self.qc == 'nn_pes' and step < max_step):
            # arguments for Gaussian
            kwargs = {
                'method': self.method,
                'basis': self.basis,
                'nprocshared': self.ppn,
                'mem': '700MW',
                'chk': job,
                'label': job,
                'Symm': 'None',
                'mult': mult,
                'charge': charge,
                'scf': 'xqc',
                'pop': 'None',
            }
            if self.par['guessmix'] == 1 or \
                'barrierless_saddle' in job or \
                'bls' in job or \
                (mult == 1 and 'R_Addition_MultipleBond' in job):
                kwargs['guess'] = 'Mix,Always'
            if ts:
                # arguments for transition state searches
                kwargs['method'] = 'am1'
                kwargs['basis'] = ''

                if step == 0:
                    if not self.par['bimol']:
                        kwargs['opt'] = 'ModRedun,Loose,CalcFC'
                    else:
                        kwargs['opt'] = 'ModRedun,Loose,CalcFC'
                        kwargs['method'] = self.method
                        kwargs['basis'] = self.basis
                elif step < max_step:
                    kwargs['opt'] = 'ModRedun,Loose,CalcFC'
                    kwargs['guess'] = 'Read'
                    if self.par['guessmix'] == 1 or \
                        'barrierless_saddle' in job or \
                        (mult == 1 and 'R_Addition_MultipleBond' in job):
                        kwargs['guess'] = 'Read,Mix,Always'
                    if self.par['bimol']:
                        kwargs['method'] = self.method
                        kwargs['basis'] = self.basis
                        kwargs['opt'] = 'QST3,AddRedundant'
                        kwargs.pop('guess', None)
                else:
                    kwargs['method'] = self.method
                    kwargs['basis'] = self.basis
                    if self.par['calcall_ts'] == 1:
                        kwargs['opt'] = 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999'
                        # not sending the frequency calculation for CalcAll
                    else:
                        kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,MaxCycle=999'
                        kwargs['freq'] = 'freq'
            elif vts:
                kwargs['method'] = self.vts_method
                kwargs['basis'] = self.vts_basis
                kwargs['opt'] = 'ModRedun,CalcFC,MaxCycles=999'
                kwargs['guess'] = 'Mix, Always'
                kwargs.pop('freq', None)
                kwargs.pop('chk', None)
            else:
                kwargs['freq'] = 'freq'  # ? 
            if (scan or 'R_Addition_MultipleBond' in job) and not vts:
                kwargs['method'] = self.scan_method 
                kwargs['basis'] = self.scan_basis
            if 'barrierless_saddle' in job or 'bls' in job:
                kwargs['method'] = self.bls_method
                kwargs['basis'] = self.bls_basis
            if irc is not None:
                # arguments for the irc calculations
                if start_from_geom == 0:
                    kwargs['geom'] = 'AllCheck,NoKeepConstants'
                    if self.par['guessmix'] == 1 or \
                        'barrierless_saddle' in job or \
                        (mult == 1 and 'R_Addition_MultipleBond' in job):
                        kwargs['guess'] = 'Read,Mix'  # Always is illegal here
                    else:
                        kwargs['guess'] = 'Read'
                    kwargs['irc'] = 'RCFC,{},MaxPoints={},StepSize={}'.format(irc, self.irc_maxpoints, self.irc_stepsize)
                else:
                    kwargs['irc'] = 'RCFC,CalcFC,{},MaxPoints={},StepSize={}'.format(irc, self.irc_maxpoints, self.irc_stepsize)
                kwargs.pop('freq', None)
            if high_level:
                kwargs['method'] = self.high_level_method
                kwargs['basis'] = self.high_level_basis
                if len(self.opt) > 0:
                    if self.par['calcall_ts'] == 1 and ts == 1:
                        kwargs['opt'] = 'NoFreeze,TS,CalcAll,NoEigentest,' \
                                        'MaxCycle=999,{}'.format(self.opt)  
                    else:
                        kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,' \
                                        'MaxCycle=999,{}'.format(self.opt)  # to overwrite possible CalcAll
                        kwargs['freq'] = 'freq'
                else:
                    if self.par['calcall_ts'] == 1 and ts == 1:
                        kwargs['opt'] = 'NoFreeze,TS,CalcAll,NoEigentest,' \
                                        'MaxCycle=999'  
                    else:
                        kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,' \
                                        'MaxCycle=999'  # to overwrite possible CalcAll
                        kwargs['freq'] = 'freq'
                if len(self.integral) > 0:
                    kwargs['integral'] = self.integral
                if 'barrierless_saddle' in job:  # completely overwrite normal settings
                    kwargs['method'] = self.bls_high_level_method
                    kwargs['basis'] = self.bls_high_level_basis
                    kwargs['opt'] = 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999'  # to overwrite possible CalcAll
                    kwargs.pop('freq', None)
            if hir:
                kwargs['opt'] = 'ModRedun,CalcFC'
                #if (not ts) or (ts and (not self.par['calcall_ts'])):
                kwargs.pop('freq', None)
                if ts:
                    if self.par['hir_maxcycle'] is None:
                        kwargs['opt'] = 'ModRedun,CalcFC,TS,NoEigentest'
                    else:
                        kwargs['opt'] = 'ModRedun,CalcFC,TS,NoEigentest,MaxCycle={}'.format(self.par['hir_maxcycle'])
                if rigid == 1:
                    kwargs.pop('freq', None)
                    kwargs.pop('opt', None)
            if 'hom_sci' in job and not vts:
                kwargs.pop('opt', None)
            if aie:
                kwargs = {
                    'method': 'cbs-qb3',
                    'basis': None,
                    'nprocshared': self.ppn,
                    'mem': '700MW',
                    'label': job,
                    'Symm': 'None',
                    'mult': mult,
                    'charge': charge,
                    'scf': 'xqc'
                }

            return kwargs

        elif self.qc == 'nwchem':
            # arguments for NWChem
            odft = mult > 1
            kwargs = {
                'xc': self.method,
                'basis': self.basis,
                'scratch_dir': 'scratch',
                'permanent_dir': './perm',
                'label': job,
                'mult': mult,
                'charge': charge,
                'odft': odft,
                'task': 'optimize',
                'driver': 'tight',
                'maxiter': 100,
                'clear': 1,
                'inhess': 1
            }
            if ts:
                # arguments for transition state searches
                if step == max_step:
                    kwargs['task'] = 'saddle'
                else:
                    kwargs['driver'] = 'default'
            if irc is not None:
                # arguments for the irc calculations
                irc_kwargs = {
                    'task': 'mepgs',
                    'mepgs': 1,
                    'maxmep': 30,
                    'maxiter': 20,
                    'inhess': 2,
                    'xyz': 1,
                    'evib': 0.0005,
                    'stride': 0.1,
                    'direction': irc
                }
                kwargs.update(irc_kwargs)
            return kwargs

        elif self.qc == 'qchem':
            # arguments for QChem
            kwargs = {
                'label': job,
                'jobtype': 'opt',
                'nt': self.ppn,
                'method': self.method,
                'basis': self.basis,
                'unrestricted': 'True',
                'scf_algorithm': 'diis_gdm',
                'max_scf_cycles': '100',
                'geom_opt_max_cycles': '500',
                'multiplicity': mult,
                'charge': charge

            }
            if ts:
                if 0 < step < max_step:
                    kwargs['scf_guess'] = 'read'
                if step >= max_step:
                    kwargs['jobtype'] = 'ts'
                elif not scan and 'R_Addition_MultipleBond' not in job:
                        kwargs['method'] = 'b3lyp'
                        kwargs['basis'] = 'sto-3g'
                        kwargs['scf_convergence'] = '4'
                        kwargs['geom_opt_tol_gradient'] = '1500'
                        kwargs['geom_opt_tol_displacement'] = '6000'
                        kwargs['geom_opt_tol_energy'] = '500'

            if scan or 'R_Addition_MultipleBond' in job:
                kwargs['method'] = self.scan_method
                kwargs['basis'] = self.scan_basis
            if 'barrierless_saddle' in job or 'bls' in job:
                kwargs['method'] = self.bls_method
                kwargs['basis'] = self.bls_basis
            if vts:
                kwargs['method'] = self.vts_method
                kwargs['basis'] = self.vts_basis
            if high_level:
                kwargs['method'] = self.high_level_method
                kwargs['basis'] = self.high_level_basis
                kwargs['vibman_print'] = '4'
            if irc is not None:
                kwargs['jobtype'] = 'freq'
                kwargs['xc_grid'] = '3'
                kwargs['vibman_print'] = '4'
                kwargs['rpath_max_stepsize'] = str(self.irc_stepsize * 10)
                kwargs['rpath_max_cycles'] = str(self.irc_maxpoints)
                if irc == 'forward':
                    kwargs['rpath_direction'] = '-1'
                else:
                    kwargs['rpath_direction'] = '1'
        elif self.qc == 'nn_pes':
            if self.par['nn_model']:
                kwargs = {'fname': self.par['nn_model']}
            else:
                kwargs = {}

        return kwargs

    def qc_hir(self, species, geom, rot_index, ang_index, fix, rigid):
        '''Creates a constrained geometry optimization input and runs it.
        wellorts: 0 for wells and 1 for saddle points
        rot_index: index of the rotor in the molecule
        ang_index: index for the current size of the angle
        fix: four atoms of the dihedral that is currently fixed
        '''
        from kinbot.geometry import calc_dihedral
        if species.wellorts:
            job = 'hir/' + species.name + '_hir_' + str(rot_index) + '_' + str(ang_index).zfill(2)
        else:
            job = 'hir/' + str(species.chemid) + '_hir_' + str(rot_index) + '_' + str(ang_index).zfill(2)

        kwargs = self.get_qc_arguments(job, species.mult, species.charge, 
                                       ts=species.wellorts, step=1, max_step=1,
                                       high_level=1, hir=1, rigid=rigid)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        if self.qc == 'gauss': 
            code = 'gaussian'
            Code = 'Gaussian'
            if not self.use_sella:
                kwargs['addsec'] = f'{" ".join(str(f) for f in fix[0])} F'
                kwargs.pop('chk', None)
        elif self.qc == 'qchem':
            code = 'qchem'
            Code = 'QChem'
            if not self.use_sella:
                dihedral = calc_dihedral(geom[fix[0][0]-1], geom[fix[0][1]-1], 
                                         geom[fix[0][2]-1], geom[fix[0][3]-1])[0]
                kwargs['addsec'] = '\n$opt\nCONSTRAINT\ntors ' \
                                   f'{" ".join(str(f) for f in fix[0])} ' \
                                   f'{dihedral}\nENDCONSTRAINT\n$end'
        elif self.qc == 'nwchem':
            code = 'nwchem'
            Code = 'NWChem'
        elif self.qc == 'nn_pes':
            code = 'nn_pes'
            Code = 'Nn_surr'
        else:
            raise ValueError(f'Unexpected value for qc parameter: {self.qc}')
        if self.use_sella:
            kwargs.pop('chk', None)
            kwargs.pop('opt', None)
            kwargs.pop('freq', None)
            template_file = f'{kb_path}/tpl/ase_sella_hir.tpl.py'
        else:
            template_file = f'{kb_path}/tpl/ase_{self.qc}_hir.tpl.py'
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(species.atom),
                                   geom=list([list(gi) for gi in geom]),
                                   ppn=self.ppn,
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   fix=fix[0],
                                   code=code,  # Sella
                                   Code=Code,  # Sella
                                   order=species.wellorts,  # Sella
                                   sella_kwargs=self.par['sella_kwargs']  # Sella
                                   )

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job)

        return 0

    def qc_ring_conf(self, species, geom, fix, change, conf_idx, dih_idx):
        '''
        Creates a constrained geometry optimization input for the
        conformational search of cyclic structures and runs it.
        Makes use of Sella.

        qc: 'gauss' or 'nwchem' or 'qchem'
        wellorts: 0 for wells and 1 for saddle points
        conf_idx: index of ring conformer in the conformer search
        dih_idx: index of the dih for this conformer
        '''
        if species.wellorts:
            job = 'conf/' + species.name + '_r' + str(conf_idx).zfill(self.zf) \
                  + '_' + str(dih_idx).zfill(self.zf)
        else:
            job = 'conf/' + str(species.chemid) + '_r' \
                  + str(conf_idx).zfill(self.zf) + '_' \
                  + str(dih_idx).zfill(self.zf)

        kwargs = self.get_qc_arguments(job, species.mult, species.charge, 
                                       ts=species.wellorts, step=1, max_step=1, 
                                       hir=1)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        kwargs.pop('opt', None)
        kwargs.pop('chk', None)
        kwargs['method'] = 'am1'
        kwargs['basis'] = ''

        if self.qc == 'gauss':
            code = 'gaussian'
            Code = 'Gaussian'
        elif self.qc == 'qchem':
            code = 'qchem'
            Code = 'QChem'
        elif self.qc == 'nwchem':
            code = 'nwchem'
            Code = 'NWChem'
        elif self.qc == 'nn_pes':
            kwargs.pop('method', None)
            kwargs.pop('basis', None)
            code = 'nn_pes'
            Code = 'Nn_surr'
        else:
            raise ValueError(f'Unexpected value for qc parameter: {self.qc}')
        template_file = f'{kb_path}/tpl/ase_sella_ring_conf.tpl.py'
        # if self.use_sella:
        #     template_file = f'{kb_path}/tpl/ase_sella_ring_conf.tpl.py'
        # else:
        #     template_file = f'{kb_path}/tpl/ase_{self.qc}_ring_conf.tpl.py'
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(species.atom),
                                   geom=list([list(gi) for gi in geom]),
                                   fix=fix,
                                   change=change,
                                   ppn=self.ppn,
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   code=code,  # Sella
                                   Code=Code,  # Sella
                                   sella_kwargs=self.par['sella_kwargs']  # Sella
                                   )

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job)
        return 0

    def qc_conf(self, species, geom, index, semi_emp=0):
        '''
        Creates a geometry optimization input for the conformational search and runs it.
        qc: 'gauss' or 'nwchem' or 'qchem'
        wellorts: 0 for wells and 1 for saddle points
        index: >=0 for sampling, each job will get numbered with index
        '''
        add = ''
        if semi_emp:
            add = 'semi_emp_'
        if species.wellorts:
            job = 'conf/' + species.name + '_' + add + str(index).zfill(self.zf)
        else:
            job = 'conf/' + str(species.chemid) + '_' + add \
                  + str(index).zfill(self.zf)

        if species.wellorts:
            kwargs = self.get_qc_arguments(job, species.mult, species.charge, 
                                           ts=1, step=1, max_step=1)
        else:
            kwargs = self.get_qc_arguments(job, species.mult, species.charge)
            if self.qc == 'gauss' and not self.use_sella:
                if self.par['opt'].casefold() == 'Tight'.casefold(): 
                    kwargs['opt'] = 'CalcFC, Tight'
                else:
                    kwargs['opt'] = 'CalcFC'

        if self.qc == 'gauss':
            code = 'gaussian'
            Code = 'Gaussian'
            kwargs.pop('chk', None)
        elif self.qc == 'qchem':
            code = 'qchem'
            Code = 'QChem'
        elif self.qc == 'nwchem':
            code = 'nwchem'
            Code = 'NWChem'
        elif self.qc == 'nn_pes':
            code = 'nn_pes'
            Code = 'Nn_surr'
        else:
            raise ValueError(f'Unexpected value for qc parameter: {self.qc}')
        
        if semi_emp:
            kwargs['method'] = self.par['semi_emp_method']
            kwargs['basis'] = ''
        if species.natom < 3:
            kwargs.pop('Symm', None)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        if self.use_sella:
            kwargs.pop('opt', None)
            kwargs.pop('freq', None)
            template_file = f'{kb_path}/tpl/ase_sella_opt_well.tpl.py'
        else:        
            template_file = f'{kb_path}/tpl/ase_{self.qc}_opt_well.tpl.py'
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(species.atom),
                                   geom=list([list(gi) for gi in geom]),
                                   ppn=self.ppn,  # QChem and NWChem
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   code=code,    # Sella
                                   Code=Code,    # Sella
                                   order=species.wellorts,      # Sella
                                   sella_kwargs=self.par['sella_kwargs']  # Sella
                                   )    
        
        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job)

        return 0

    def qc_aie(self, species, geom, ext):
        '''
        Sets up AIE calculations. Currently results are not digested by KinBot.
        Only for well, not for saddles!
        qc: 'gauss' or 'nwchem' or 'qchem'
        index: the index of the conformer as in conf search
        '''
        job0 = f'aie/{str(species.chemid)}_AIE0_{ext}'  # neutral
        job1 = f'aie/{str(species.chemid)}_AIE1_{ext}'  # cation
        kwargs0 = self.get_qc_arguments(job0, species.mult, species.charge, aie=1)
        if species.mult == 1: m1 = 2
        elif species.mult == 2: m1 = 1
        elif species.mult == 3: m1 = 2
        kwargs1 = self.get_qc_arguments(job1, m1, species.charge + 1, aie=1)
        template_file = f'{kb_path}/tpl/ase_{self.qc}_opt_well.tpl.py'
        template = open(template_file, 'r').read()
        t0 = template.format(label=job0,
                             kwargs=kwargs0,
                             atom=list(species.atom),
                             geom=list([list(gi) for gi in geom]),
                             ppn=self.ppn,
                             qc_command=self.qc_command,
                             working_dir=os.getcwd())
        t1 = template.format(label=job1,
                             kwargs=kwargs1,
                             atom=list(species.atom),
                             geom=list([list(gi) for gi in geom]),
                             ppn=self.ppn,
                             qc_command=self.qc_command,
                             working_dir=os.getcwd())
        
        with open(f'{job0}.py', 'w') as f:
            f.write(t0)
        with open(f'{job1}.py', 'w') as f:
            f.write(t1)
        self.submit_qc(job0)
        self.submit_qc(job1)

        return 0

    def qc_opt(self, species, geom, high_level=0, mp2=0, bls=0, ext=None, 
               fdir=None, do_vdW=False):
        '''
        Creates a geometry optimization input and runs it.
        '''
        if do_vdW:
            if high_level:
                job = f'{species.name}_high'
        elif ext is None:
            job = str(species.chemid) + '_well'
            if high_level:
                job = str(species.chemid) + '_well_high'
            if mp2:
                job = str(species.chemid) + '_well_mp2'
            if bls:
                job = str(species.chemid) + '_well_bls'

        else:
            job = str(species.chemid) + ext

        if fdir is not None:
            job = f'{fdir}/{job}'

        # TODO: Code exceptions into their own function/py script that opt can call.
        # TODO: Fix symmetry numbers for calcs as well if needed
        mult = exceptions.get_multiplicity(species.chemid, species.mult)

        kwargs = self.get_qc_arguments(job, mult, species.charge,
                                       high_level=high_level)
        if 'opt' not in kwargs:
            kwargs['opt'] = ''

        if self.qc == 'gauss':
            code = 'gaussian'
            Code = 'Gaussian'
            if self.par['opt'].casefold() == 'Tight'.casefold(): 
                kwargs['opt'] = 'CalcFC, Tight'
            else:
                kwargs['opt'] = 'CalcFC'
        elif self.qc == 'qchem':
            code = 'qchem'
            Code = 'QChem'
        elif self.qc == 'nwchem':
            code = 'nwchem'
            Code = 'NWChem'   
        elif self.qc == 'nn_pes':
            code = 'nn_pes'
            Code = 'Nn_surr'
        else:
            raise ValueError(f'Unexpected value for qc parameter: {self.qc}')
        
        if mp2:
            kwargs['method'] = self.scan_method
            kwargs['basis'] = self.scan_basis
        if high_level and self.qc == 'gauss' and self.opt:
            kwargs['opt'] = 'CalcFC, {}'.format(self.opt)
        if species.natom < 3:
            kwargs.pop('Symm', None)
        # the integral is set in the get_qc_arguments parts, bad design
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        if self.use_sella:
            kwargs.pop('opt', None)
            kwargs.pop('freq', None)
            template_file = f'{kb_path}/tpl/ase_sella_opt_well.tpl.py'
        else:
            template_file = f'{kb_path}/tpl/ase_{self.qc}_opt_well.tpl.py'
        
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(species.atom),
                                   geom=list([list(gi) for gi in geom]),
                                   ppn=self.ppn,  # QChem and NWChem
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   code=code,    # Sella
                                   Code=Code,    # Sella
                                   order=0,      # Sella
                                   sella_kwargs=self.par['sella_kwargs'] # Sella
                                   )

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job)
        return 0

    def qc_opt_ts(self, species, geom, high_level=0, ext=None, fdir=None):
        '''Creates a ts optimization input and runs it
        '''

        job = str(species.name)
        if ext is None:
            if high_level:
                job += '_high'
        else:
            job += ext

        if fdir is not None:
            job = f'{fdir}/{job}'

        kwargs = self.get_qc_arguments(job, species.mult, species.charge, ts=1, 
                                       step=1, max_step=1, high_level=1)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        if self.qc == 'gauss':
            code = 'gaussian'
            Code = 'Gaussian'
        elif self.qc == 'qchem':
            code = 'qchem'
            Code = 'QChem'
        elif self.qc == 'nwchem':
            code = 'nwchem'
            Code = 'NWChem'
        elif self.qc == 'nn_pes':
            code = 'nn_pes'
            Code = 'Nn_surr'
        else:
            raise ValueError(f'Unrecognized qc option: {self.qc}')
        
        if self.use_sella:
            kwargs.pop('addsec', None)
            kwargs.pop('opt', None)
            kwargs.pop('freq', None)
            template_file = f'{kb_path}/tpl/ase_sella_ts_end.tpl.py'
        else:
            template_file = f'{kb_path}/tpl/ase_{self.qc}_ts_end.tpl.py'
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(species.atom),
                                   geom=list([list(gi) for gi in geom]),
                                   ppn=self.ppn,
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   code=code, # Sella
                                   Code=Code, # Sella
                                   sella_kwargs=self.par['sella_kwargs'] # Sella
                                   )

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job)

        return 0

    def qc_vts_frag(self, frag):
        '''
        Runs VTS fragment optimization
        '''

        job = f'vrctst/{str(frag.chemid)}_vts'
        mult = exceptions.get_multiplicity(frag.chemid, frag.mult)
        kwargs = self.get_qc_arguments(job, mult, frag.charge, vts=1)
        # Add chk to fragments only
        kwargs['chk'] = f'{str(frag.chemid)}_vts'

        if self.qc != 'gauss':
            raise ValueError(f'Only implemeted for Gaussian. Instead I got: {self.qc}')

        if frag.natom < 3:
            kwargs.pop('Symm', None)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)

        template_file = f'{kb_path}/tpl/ase_{self.qc}_opt_well.tpl.py'
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(frag.atom),
                                   geom=list([list(gi) for gi in frag.geom]),
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd(),
                                   )  # TODO add sella keywords

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job)
        
        return job 

    def qc_vts(self, reac, geom, step, equiv, asymptote, step0_geom):
        '''
        Creates a geometry optimization along a scan and runs it.
        reac: full reaction object
        '''

        if not asymptote:
            job = f'vrctst/{reac.instance_name}_vts_pt{str(step).zfill(2)}'
        else:
            job = f'vrctst/{reac.instance_name}_vts_pt_asymptote'
        mult = exceptions.get_multiplicity(reac.species.chemid, reac.species.mult)
        kwargs = self.get_qc_arguments(job, mult, reac.species.charge, vts=1)

        if self.qc == 'gauss' and not self.par['vrc_tst_scan_sella']:
            kwargs['addsec'] = f'{reac.scan_coo[0]+1} {reac.scan_coo[1]+1} F\n'
        
        if reac.species.natom < 3:
            kwargs.pop('Symm', None)
        if self.par['calc_kwargs']:
            kwargs = self.merge_kwargs(kwargs)
        
        if self.par['vrc_tst_scan_sella']:  # TODO sella globally
            kwargs.pop('opt', None)
            if self.qc == 'gauss':
                code = 'gaussian'
                Code = 'Gaussian'
            else:
                raise ValueError(f'Currently only Gaussian is supported: {self.qc}')
            template_file = f'{kb_path}/tpl/ase_sella_vts.tpl.py'
            template = open(template_file, 'r').read()
            template = template.format(label=job,
                                       scan_coo=reac.scan_coo,
                                       bonds=reac.irc_prod.bondlist,
                                       kwargs=kwargs,
                                       atom=list(reac.species.atom),
                                       init_geom=list([list(gi) for gi in geom]),
                                       qc_command=self.qc_command,
                                       working_dir=os.getcwd(),
                                       vts_ang_dev=self.par['vts_ang_dev'],
                                       scan_deviation=self.par['vrc_tst_scan_deviation'],
                                       frag_maps=list([list(rm) for rm in reac.maps]),
                                       froz_A_geom=list([list(gi) for gi in reac.products[0].geom]),
                                       froz_A_atom=reac.products[0].atom,
                                       froz_B_geom=list([list(gi) for gi in reac.products[1].geom]),
                                       froz_B_atom=reac.products[1].atom,
                                       code=code,
                                       Code=Code,
                                       sella_kwargs=self.par['sella_kwargs'],  # Sella
                                       equiv=equiv,
                                       asymptote=asymptote,
                                       frag_bonds_0=list([list(b) for b in reac.products[0].bond01]),
                                       frag_bonds_1=list([list(b) for b in reac.products[1].bond01]),
                                       step0_geom=step0_geom
                                       )

        else:
            template_file = f'{kb_path}/tpl/ase_{self.qc}_vts.tpl.py'
            template = open(template_file, 'r').read()
            template = template.format(label=job,
                                       scan_coo=reac.scan_coo,
                                       bonds=reac.irc_prod.bondlist,
                                       kwargs=kwargs,
                                       atom=list(reac.species.atom),
                                       init_geom=list([list(gi) for gi in geom]),
                                       qc_command=self.qc_command,
                                       working_dir=os.getcwd(),
                                       scan_deviation=self.par['vrc_tst_scan_deviation'],
                                       frag_maps=list([list(rm) for rm in reac.maps]),
                                       froz_A_geom=list([list(gi) for gi in reac.products[0].geom]),
                                       froz_A_atom=reac.products[0].atom,
                                       froz_B_geom=list([list(gi) for gi in reac.products[1].geom]),
                                       froz_B_atom=reac.products[1].atom,
                                       )

        with open(f'{job}.py', 'w') as f:
            f.write(template)

        self.submit_qc(job)
        return job 

    def submit_qc(self, job, singlejob=1, jobtype=None):
        '''Submit a job to the queue, unless the job:
            * has finished with Normal termination
            * has finished with Error termination
            * is currently running
        However, if the optional parameter singlejob is set to zero, then
        the job is run only if it has finished earlier with normal termination.
        This is for continuations, when the continuing jobs overwrite each other.
        If the number of jobs in the queue is larger than the user-set limit,
        KinBot will park here until resources are freed up.
        '''
        # if the logfile already exists, copy it with another name

        check = self.check_qc(job)
        if singlejob == 1:
            if check != 0:
                return 0
        else:
            if check == 'running':
                return 0

        if self.queue_job_limit > 0:
            self.limit_jobs()

        try:
            if jobtype == 'am1' and self.par['q_temp_am1']:
                template_head_file = self.par['q_temp_am1']
            elif job.endswith('mp2') and self.par['q_temp_mp2']:
                template_head_file = self.par['q_temp_mp2']
            elif job.endswith('high') and self.par['q_temp_hi']:
                template_head_file = self.par['q_temp_hi']
            elif self.par['queue_template'] == '':
                template_head_file = f'{kb_path}/tpl/{self.queuing}.tpl'
            else:
                template_head_file = self.par['queue_template']
        except OSError:
            logger.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logger.error('Or no file is found at {}.'.format(self.par['queue_template']))
            logger.error('Exiting')
            sys.exit()

        if self.queuing == 'local':
            err_msg = f'Job {job} is missing in the database or the output ' \
                      'file is not present. Unable to run calculations when ' \
                      'queuing is \'local\'. To ignore calculations that have ' \
                      'not finished set "error_missing_local" to False.'
            if self.par['error_missing_local']:
                logger.error(err_msg)
                raise FileNotFoundError(err_msg)
            else:
                logger.warning(err_msg)
                return -1

        template_file = f'{kb_path}/tpl/{self.queuing}_python.tpl'
        python_file = f'{job}.py'
        job_template = open(template_head_file, 'r').read() + open(template_file, 'r').read()

        if self.queuing == 'pbs':
            job_template = job_template.format(name=job, ppn=self.ppn, queue_name=self.queue_name,
                                               errdir='perm', python_file=python_file, arguments='')
        elif self.queuing == 'slurm':
            job_template = job_template.format(name=job, ppn=self.ppn, queue_name=self.queue_name, errdir='perm',
                                               slurm_feature=self.slurm_feature, python_file=python_file, arguments='')
        else:
            logger.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logger.error('Exiting')
            sys.exit()

        qu_file = '{}{}'.format(job, constants.qext[self.queuing])
        with open(qu_file, 'w') as f_out_qu:
            f_out_qu.write(job_template)

        command = [constants.qsubmit[self.queuing], job + constants.qext[self.queuing]]
        process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = process.communicate()
        out = out.decode()
        err = err.decode()
        try:
            if self.queuing == 'pbs':
                pid = out.split('\n')[0].split('.')[0]
            elif self.queuing == 'slurm':
                pid = out.split('\n')[0].split()[3]
        except:
            msg = 'Something went wrong when submitting a job'
            msg += 'This is the standard output:\n' + out
            msg += '\nThis is the standard error:\n' + err
            logger.error(msg)
            sys.exit()
        self.job_ids[job] = pid

        now = datetime.now()
        logger.debug(f'SUBMITTED {job} on {now.ctime()}')
        return 1  # important to keep it 1, this is the natural counter of jobs submitted

    def get_qc_geom(self,
                    job,
                    natom,
                    wait=0,
                    allow_error=0,
                    previous=0,
                    reorder: bool = False):
        '''
        Get the geometry from the ase database file.
        Returns it, with the following conditions about the status of the job.
        If wait = 0, return an (1, empty array) when the job is still running (instead of the geometry).
        If wait = 1, wait for the job to finish.
        If wait = 2, return the last geometry while the job is still running.
            This option is to monitor the formation of bimolecular products.
        If allow_error = 0, do not read the final geometry if the status is not 'normal'
        if allow_error = 1, read the geometry even though there is an error in the output file
            This option is to read the final IRC geometry when it did not converge
        If previous = 0, read the last geometry, this is the normal behaviour
        if previous = 1, read the geometry before the last one, this is needed in scan types so
            that the max energy point is taken, not the one after that
        '''
        geom = np.zeros((natom, 3))
        atoms = np.full(natom, 'H')
        check = self.check_qc(job)
        if check == 'error' and not allow_error:
            if reorder:
                return -1, geom, atoms
            else:
                return -1, geom
        status = 0
        while 1:
            check = self.check_qc(job)
            if check == 'running':
                if wait == 1:
                    time.sleep(1)
                elif wait == 2:
                    status = 2
                    break
                else:
                    if reorder:
                        return 1, geom, atoms
                    else:
                        return 1, geom
            else:
                break
        if check != 'normal':
            if not allow_error:
                if wait != 2:
                    if reorder:
                        return -1, geom, atoms
                    else:
                        return -1, geom

        # open the database
        rows = self.db.select(name=job)

        found_entry = 0
        # take the last entry
        for row in rows:
            if found_entry:
                prev_geom = geom  # saves the previous 
            mol = row.toatoms()
            geom = mol.positions
            atoms = mol.symbols
            found_entry = 1

        if found_entry and previous == 0:
            if reorder:
                return status, geom, atoms
            else:
                return status, geom
        elif found_entry and previous == 1:
            return status, prev_geom
        else:
            if reorder:
                return -1, geom, atoms
            else:
                return -1, geom

    def get_qc_freq(self, job, natom, wait=0, allow_error=0):
        '''
        Get the frequencies from the ase database file
        If wait is set to 1, it will wait for the job to finish.
        '''

        check = self.check_qc(job)
        if check == 'error':
            return -1, [0]
        while 1:
            check = self.check_qc(job)
            if check == 'running':
                if wait == 1:
                    time.sleep(1)
                else:
                    return 1, []
            else:
                break

        if check != 'normal':
            return -1, [0]

        freq = []

        # open the database
        rows = self.db.select(name=job)

        # take the last entry
        for row in rows:
            if hasattr(row, 'data'):
                if not row.data.get('frequencies') is None:
                    freq = list(row.data.get('frequencies'))

        if len(freq) == 0 and natom > 1:
            return -1, freq

        return 0, freq

    def get_qc_energy(self, job, wait=0):
        '''
        Read the last energy from a job.
        For Gaussian currently works for DFT and HF only.
        For NWChem it works for optimization jobs, using the @ handle.
        If wait is set to 1, it will wait for the job to finish, otherwise
        just reads the last energy in the file.
        Returns the error code and the energy
        Error codes:
        -1: error
         0: normal and done
         1: running
        '''

        check = self.check_qc(job)
        if check == 'error':
            return -1, 0.
        while 1:
            check = self.check_qc(job)
            if check == 'running':
                if wait == 1:
                    time.sleep(1)
                else:
                    return 1, 0.
            else:
                break

        # Get last entry
        *_, last_row = self.db.select(name=job)
        if hasattr(last_row, 'data') and 'energy' in last_row.data:
            energy = last_row.data.get('energy')
        else:
            logger.warning(f'No energy found in the database for {job}. '
                           'This will lead to erroneous energies.')
            energy = 0

        # ase energies are always in ev, convert to hartree
        energy *= constants.EVtoHARTREE

        return 0, energy

    def get_qc_zpe(self, job, wait=1):
        '''
        Read the zero point energy.
        If wait is set to 1 (default), it will wait for the job to finish.
        '''

        check = self.check_qc(job)
        if check == 'error':
            return -1, 0.
        while 1:
            check = self.check_qc(job)
            if check == 'running':
                if wait == 1:
                    time.sleep(1)
                else:
                    return 0, 0.
            else:
                break
        zpe = 0.0  # set as default
        # open the database
        *_, last_row = self.db.select(name=job)
        # take the last entry
        
        if hasattr(last_row, 'data'):
            zpe = last_row.data.get('zpe')
        else:
            zpe = 0.0
            logger.warning('{} has no zpe in database. ZPE SET TO 0.0'.format(job))

        if zpe is None:
            zpe = 0.00

        return 0, zpe

    def read_qc_hess(self, job, natom):
        '''
        Read the hessian of a gaussian chk file
        '''

        check = self.check_qc(job)
        if check != 'normal':
            return []

        if self.use_sella:
            db = connect('kinbot.db')
            for row in db.select(name=job):
                last_row = row
            hess = last_row.data['hess']
        elif self.qc == 'gauss':
            hess = np.zeros((3 * natom, 3 * natom))
            fchk = str(job) + '.fchk'
            if not os.path.exists(fchk):
                # create the fchk file using formchk
                os.system('formchk ' + job + '.chk > /dev/null')

            with open(fchk) as f:
                lines = f.read().split('\n')
            nvals = 3 * natom * (3 * natom + 1) / 2

            for index, line in enumerate(reversed(lines)):
                if re.search('Cartesian Force Constants', line) != None:
                    hess_flat = []
                    n = 0
                    while len(hess_flat) < nvals:
                        hess_flat.extend([float(val) for val in lines[-index + n].split()])
                        n += 1
                    n = 0
                    for i in range(3 * natom):
                        for j in range(i + 1):
                            hess[i][j] = hess_flat[n]
                            hess[j][i] = hess_flat[n]
                            n += 1
                    break
            else:
                # Try to see if the Hessian is the database.
                db = connect('kinbot.db')
                for row in db.select(name=job):
                    last_row = row
                if 'hess' in last_row.data:
                    hess = last_row.data['hess']
                    logger.warning(f'Hessian matrix not found on {fchk}. '
                                   'Reading from kinbot.db database.')
                else:
                    raise FileNotFoundError(f'Hessian matrix not found on {fchk}.')
        elif self.qc == 'qchem':
            hess = []
            row = 0
            do_read = False
            if natom == 1:
                return np.zeros([3, 3])
            with open(job + '_freq.out') as f:
                for line in f:
                    if not do_read and line.startswith(' Mass-Weighted Hessian Matrix'):
                        do_read = True
                    elif do_read and 'Translations and Rotations' in line:
                        do_read = False
                    elif do_read:
                        if len(line) > 5:
                            if len(hess) < 3 * natom:
                                hess.append([])
                            hess[row].extend([float(val) for val in line.split()])
                            row += 1
                        else:
                            row = 0
                    else:
                        continue
            hess = np.array(hess)
        else:
            raise NotImplementedError()
        return hess

    def is_in_database(self, job):
        '''
        Checks if the current job is in the database:
        '''
        # open the database
        rows = self.db.select(name=job)

        mol = None

        # take the last entry
        for row in rows:
            mol = row.toatoms()

        if mol is None:
            return 0

        return 1

    def check_qc(self, job):
        '''
        Checks the status of the qc job.
        Possible returns:
        running - the job is either running or is in the queue
        status - this can be normal or error, read from the database.
                 Only happens if log file had a done stamp and it was in the db.
        0 - job is not in the db or log file is not there with a done stamp or both.
            ==> this one resets the step number to 0
        '''
        logger.debug('Checking job {}'.format(job))
        devnull = open(os.devnull, 'w')
        if self.queuing == 'pbs':
            command = 'qstat -f | grep ' + '"Job Id: ' + self.job_ids.get(job, '-1') + '"' + ' > /dev/null'
            if int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull)) == 0:
                logger.debug('Job is running')
                return 'running'
        elif self.queuing == 'slurm':
            # command = 'scontrol show job ' + self.job_ids.get(job,'-1') + ' | grep "JobId=' + self.job_ids.get(job,'-1') + '"' + ' > /dev/null'
            command = 'squeue'
            process = subprocess.Popen(command,
                                       shell=True,
                                       stdout=subprocess.PIPE,
                                       stdin=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            out, err = process.communicate()
            out = out.decode().split('\n')[1:]
            for line in out:
                if len(line) > 0:
                    pid = line.split()[0]
                    if pid == self.job_ids.get(job, '-1'):
                        logger.debug('Job is running')
                        return 'running'
        elif self.queuing == 'local':
            pass
        else:
            logger.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logger.error('Exiting')
            sys.exit()
        # if int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull)) == 0:
        #     return 'running'

        logger.debug('Checking for job {} in db'.format(job))
        if self.is_in_database(job):
            logger.debug('{} is in db'.format(job))
            for i in range(1):
                if self.qc == 'gauss':
                    log_file = job + '.log'
                elif self.qc == 'nwchem':
                    log_file = job + '.out'
                elif self.qc == 'qchem':
                    log_file = job + '.out'
                elif self.qc == 'nn_pes':
                    log_file_exists = False
                else:
                    raise ValueError('Unknown code')
                if self.qc != 'nn_pes':
                    log_file_exists = os.path.exists(log_file)
                if log_file_exists:
                    with open(log_file, 'r') as f:
                        try:
                            last_line = f.readlines()[-1]
                            if 'done' not in last_line:
                                if self.queuing == 'local' and not self.par['error_missing_local']:
                                    return 'error'
                                logger.debug(f'Log file {log_file} is present, '
                                             'but has no "done" stamp.')
                                return 0
                        except IndexError:
                            logger.debug(f'Log file {log_file} is present, but it is empty.')
                            if self.queuing == 'local' and not self.par['error_missing_local']:
                                return 'error'
                            return 0
                    logger.debug('Log file is present after {} iterations'.format(i))
                elif self.qc == 'nn_pes':
                    pass
                else:
                    if self.queuing == 'local' and not self.par['error_missing_local']:
                        return 'error'
                    logger.debug('Checking againg for log file')
                    time.sleep(1)
                    log_file_exists = os.path.exists(log_file)
                    if not log_file_exists:
                        logger.debug(f'Log file {log_file} is not present.')
                        return 0
            
                # by deleting a log file, you allow restarting a job
                # open the database
                rows = self.db.select(name=job)
                data = None
                # take the last entry
                for row in rows:
                    if hasattr(row, 'data'):
                        data = row.data
                if data is None:
                    logger.debug('Data is not in database...')
                    return 0
                else:
                    logger.debug('Returning status {}'.format(data['status']))
                    return data['status']

            logger.debug('log file {} does not exist'.format(log_file))
            return 0
        else:
            logger.debug('job {} is not in database'.format(job))
            if self.queuing == 'local' and not self.par['error_missing_local']:
                logger.debug('local qu error')
                return 'error'
            return 0

    def limit_jobs(self):
        '''
        Check how many jobs are in the queue from the user, and if larger than the limit,
        then wait for resources to free up.
        '''

        while 1:
            if self.queuing == 'slurm':
                command = ['squeue', '-h', '-u', '{}'.format(self.username)]
            elif self.queuing == 'pbs':
                command = ['qselect', '-u', '{}'.format(self.username)]
            elif self.queuing == 'local':
                command = command = ['echo', '']
            jobs = subprocess.check_output(command)

            if len(jobs.split(b'\n')) < self.queue_job_limit:
                return 0
            time.sleep(1)

    def add_dummy(self, spatom, geom, spbond):
        '''
        Add a dummy atom for each close to linear angle.
        Obsolete.
        '''
        atom = copy.deepcopy(spatom)
        dummy = geometry.is_linear(geom, spbond)
        if len(dummy) > 0:
            for d in dummy:
                atom = np.append(atom, ['X'])
                geom = np.concatenate((geom, [d]), axis=0)
        dummy = [d.tolist() for d in dummy]
        return atom, geom, dummy

    def merge_kwargs(self, kwargs):
        user_kwargs = self.par['calc_kwargs']
        for key, uvalue in user_kwargs.items():
            if key in kwargs:
                user_val_dict, val_dict = {}, {}
                for field in uvalue.split(','):
                    if len(field.split('=')) == 1:
                        user_val_dict[field] = None
                    elif len(field.split('=')) == 2:
                        user_val_dict[field.split('=')[0]] = field.split('=')[1]
                    else:
                        raise ValueError(f'Unexpected value: {uvalue}')
                for field in kwargs[key].split(','):
                    if len(field.split('=')) == 1:
                        val_dict[field] = None
                    elif len(field.split('=')) == 2:
                        val_dict[field.split('=')[0]] = field.split('=')[1]
                    else:
                        raise ValueError(f'Unexpected value: {uvalue}')
                val_dict.update(user_val_dict)
                kwargs[key] = ','.join([f'{k}={v}' if v is not None else k for k, v in val_dict.items() ])
                
            else:
                kwargs[key] = uvalue
        return kwargs
