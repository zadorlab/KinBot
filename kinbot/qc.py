import os
import sys
import subprocess
import logging
import numpy as np
import re
import time
import copy
import pkg_resources
from shutil import copyfile

from ase.db import connect
from kinbot import constants
from kinbot import geometry


class QuantumChemistry:
    """
    This class provides the link between KinBot and the qc code
    It choses the write options, submits the jobs and checks
    the jobs for success or failure
    """

    def __init__(self, par):
        self.par = par
        self.qc = par.par['qc']
        self.method = par.par['method']
        self.basis = par.par['basis']
        self.scan_method = par.par['scan_method']
        self.scan_basis = par.par['scan_basis']
        self.bls_method = par.par['barrierless_saddle_method']
        self.bls_basis = par.par['barrierless_saddle_basis']
        self.high_level_method = par.par['high_level_method']
        self.high_level_basis = par.par['high_level_basis']
        self.bls_high_level_method = par.par['barrierless_saddle_method_high']
        self.bls_high_level_basis = par.par['barrierless_saddle_basis_high']
        self.integral = par.par['integral']
        self.opt = par.par['opt']
        self.ppn = par.par['ppn']
        self.queuing = par.par['queuing']
        self.queue_name = par.par['queue_name']
        self.slurm_feature = par.par['slurm_feature']
        self.zf = par.par['zf']
        self.db = connect('kinbot.db')
        self.job_ids = {}
        self.irc_maxpoints = par.par['irc_maxpoints']
        self.irc_stepsize = par.par['irc_stepsize']
        self.qc_command = par.par['qc_command']
        if par.par['slurm_feature'] == '':
            self.slurm_feature = ''
        else:
            self.slurm_feature = '#SBATCH -C ' + par.par['slurm_feature']
        self.queue_job_limit = par.par['queue_job_limit']
        self.username = par.par['username']

    def get_qc_arguments(self, job, mult, charge, ts=0, step=0, max_step=0, irc=None, scan=0,
                         high_level=0, hir=0, start_from_geom=0):
        """
        Method to get the argument to pass to ase, which are then passed to the qc codes.
        Job: name of the job
        mult: multiplicity of the job
        charge: charge of the job
        ts: 1 for transition state searches, 0 for wells and products
        step: Step number in the transition state search
        max_step: total number of necessary steps to get to the final transition state structure
            This is different for every reaction family
        irc: direction of the irc, None if this is not an irc job
        scan: is this calculation part of a scan of a bond length to find a maximum energy
        """
        if self.qc == 'gauss':
            # arguments for Gaussian
            kwargs = {
                'method': self.method,
                'basis': self.basis,
                'nprocshared': self.ppn,
                'mem': '700MW',
                'chk': job,
                'label': job,
                'NoSymm': 'NoSymm',
                'multiplicity': mult,
                'charge': charge,
                'scf': 'xqc'
            }
            if self.par.par['guessmix'] == 1 or 'barrierless_saddle' in job:
                kwargs['guess'] = 'Mix,Always'
            if ts:
                # arguments for transition state searches
                kwargs['method'] = 'am1'
                kwargs['basis'] = ''

                if step == 0:
                    kwargs['opt'] = 'ModRedun,Tight,CalcFC,MaxCycle=999'
                elif step < max_step:
                    kwargs['opt'] = 'ModRedun,Tight,CalcFC,MaxCycle=999'
                    kwargs['guess'] = 'Read'
                    if self.par.par['guessmix'] == 1 or 'barrierless_saddle' in job:
                        kwargs['guess'] = 'Read,Mix'
                else:
                    kwargs['method'] = self.method
                    kwargs['basis'] = self.basis
                    if self.par.par['calcall_ts'] == 1:
                        kwargs['opt'] = 'NoFreeze,TS,CalcAll,NoEigentest,MaxCycle=999'
                        # not sending the frequency calculation for CalcAll
                    else:
                        kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,MaxCycle=999'
                        kwargs['freq'] = 'freq'
            else:
                kwargs['freq'] = 'freq'
            if scan or 'R_Addition_MultipleBond' in job:
                kwargs['method'] = self.scan_method 
                kwargs['basis'] = self.scan_basis
            if 'barrierless_saddle' in job or 'bls' in job:
                kwargs['method'] = self.bls_method
                kwargs['basis'] = self.bls_basis
            if irc is not None:
                # arguments for the irc calculations
                if start_from_geom == 0:
                    kwargs['geom'] = 'AllCheck,NoKeepConstants'
                    if self.par.par['guessmix'] == 1 or 'barrierless_saddle' in job:
                        kwargs['guess'] = 'Read,Mix'  # Always is illegal here
                    else:
                        kwargs['guess'] = 'Read'
                    kwargs['irc'] = 'RCFC,{},MaxPoints={},StepSize={}'.format(irc, self.irc_maxpoints, self.irc_stepsize)
                else:
                    kwargs['irc'] = 'RCFC,CalcFC,{},MaxPoints={},StepSize={}'.format(irc, self.irc_maxpoints, self.irc_stepsize)
                del kwargs['freq']
            if high_level:
                if 'barrierless_saddle' in job:
                    kwargs['method'] = self.bls_high_level_method
                    kwargs['basis'] = self.bls_high_level_basis
                else:
                    kwargs['method'] = self.high_level_method
                    kwargs['basis'] = self.high_level_basis
                if len(self.opt) > 0:
                    kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,MaxCycle=999,{}'.format(self.opt)  # to overwrite possible CalcAll
                else:
                    kwargs['opt'] = 'NoFreeze,TS,CalcFC,NoEigentest,MaxCycle=999'  # to overwrite possible CalcAll
                kwargs['freq'] = 'freq'
                if len(self.integral) > 0:
                    kwargs['integral'] = self.integral
            if hir:
                kwargs['opt'] = 'ModRedun,CalcFC'
                if (not ts) or (ts and (not self.par.par['calcall_ts'])):
                    del kwargs['freq']
                if ts:
                    kwargs['opt'] = 'ModRedun,CalcFC,TS,NoEigentest,MaxCycle=999'
            return kwargs

        if self.qc == 'nwchem':
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

    def qc_hir(self, species, geom, rot_index, ang_index, fix):
        """
        Creates a constrained geometry optimization input and runs it.
        wellorts: 0 for wells and 1 for saddle points
        rot_index: index of the rotor in the molecule
        ang_index: index for the current size of the angle
        fix: four atoms of the dihedral that is currently fixed
        """
        if species.wellorts:
            job = 'hir/' + species.name + '_hir_' + str(rot_index) + '_' + str(ang_index).zfill(2)
        else:
            job = 'hir/' + str(species.chemid) + '_hir_' + str(rot_index) + '_' + str(ang_index).zfill(2)

        kwargs = self.get_qc_arguments(job, species.mult, species.charge, ts=species.wellorts, step=1, max_step=1, high_level=1, hir=1)
        kwargs['fix'] = fix
        del kwargs['chk']

        atom = copy.deepcopy(species.atom)

        dummy = geometry.is_linear(geom, species.bond)
        if len(dummy) > 0:  # add a dummy atom for each close to linear angle
            for d in dummy:
                atom = np.append(atom, ['X'])
                geom = np.concatenate((geom, [d]), axis=0)
        dummy = [d.tolist() for d in dummy]

        template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_hir.tpl.py'.format(qc=self.qc))
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(atom),
                                   geom=list([list(gi) for gi in geom]),
                                   ppn=self.ppn,
                                   dummy=dummy,
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd())

        f_out = open('{}.py'.format(job), 'w')
        f_out.write(template)
        f_out.close()

        self.submit_qc(job)

        return 0

    def qc_ring_conf(self, species, geom, fix, change, conf_nr, scan_nr):
        """
        Creates a constrained geometry optimization input for the
        conformational search of cyclic structures and runs it.
        Make use of the ASE optimizer PCOBFGS

        qc: 'gauss' or 'nwchem'
        scan: list of dihedrals to be scanned and their values
        wellorts: 0 for wells and 1 for saddle points
        conf_nr: number of the conformer in the conformer search
        scan_nr: number of the scan for this conformer
        """
        if species.wellorts:
            job = 'conf/' + species.name + '_r' + str(conf_nr).zfill(self.zf) + '_' + str(scan_nr).zfill(self.zf)
        else:
            job = 'conf/' + str(species.chemid) + '_r' + str(conf_nr).zfill(self.zf) + '_' + str(scan_nr).zfill(self.zf)

        kwargs = self.get_qc_arguments(job, species.mult, species.charge, ts=species.wellorts, step=1, max_step=1, hir=1)

        del kwargs['opt']
        del kwargs['chk']
        kwargs['method'] = 'am1'
        kwargs['basis'] = ''

        atom = copy.deepcopy(species.atom)

        dummy = geometry.is_linear(geom, species.bond)
        if len(dummy) > 0:  # add a dummy atom for each close to linear angle
            for d in dummy:
                atom = np.append(atom, ['X'])
                geom = np.concatenate((geom, [d]), axis=0)
        dummy = [d.tolist() for d in dummy]

        template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_ring_conf.tpl.py'.format(qc=self.qc))
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(atom),
                                   geom=list([list(gi) for gi in geom]),
                                   fix=fix,
                                   change=change,
                                   ppn=self.ppn,
                                   dummy=dummy,
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd())

        f_out = open('{}.py'.format(job), 'w')
        f_out.write(template)
        f_out.close()

        self.submit_qc(job)
        return 0

    def qc_conf(self, species, geom, index=-1, ring=0, semi_emp=0):
        """
        Creates a geometry optimization input for the conformational search and runs it.
        qc: 'gauss' or 'nwchem'
        wellorts: 0 for wells and 1 for saddle points
        index: >=0 for sampling, each job will get numbered with index
        """
        if index == -1:
            job = 'conf/' + str(species.chemid) + '_well'
        else:
            add = ''
            if ring:
                add = 'r'
            if semi_emp:
                add = 'semi_emp_'
            if species.wellorts:
                job = 'conf/' + species.name + '_' + add + str(index).zfill(self.zf)
            else:
                job = 'conf/' + str(species.chemid) + '_' + add + str(index).zfill(self.zf)

        if species.wellorts:
            kwargs = self.get_qc_arguments(job, species.mult, species.charge, ts=1, step=1, max_step=1)
        else:
            kwargs = self.get_qc_arguments(job, species.mult, species.charge)
            if self.qc == 'gauss':
                kwargs['opt'] = 'CalcFC, Tight'

        del kwargs['chk']
        if semi_emp:
            kwargs['method'] = self.par.par['semi_emp_method']
            kwargs['basis'] = ''
        atom = copy.deepcopy(species.atom)

        dummy = geometry.is_linear(geom, species.bond)
        if len(dummy) > 0:  # add a dummy atom for each close to linear angle
            for d in dummy:
                atom = np.append(atom, ['X'])
                geom = np.concatenate((geom, [d]), axis=0)
        dummy = [d.tolist() for d in dummy]

        template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_opt_well.tpl.py'.format(qc=self.qc))
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(atom),
                                   geom=list([list(gi) for gi in geom]),
                                   ppn=self.ppn,
                                   dummy=dummy,
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd())
        f_out = open('{}.py'.format(job), 'w')
        f_out.write(template)
        f_out.close()

        self.submit_qc(job)

        return 0

    def qc_opt(self, species, geom, high_level=0, mp2=0, bls=0):
        """
        Creates a geometry optimization input and runs it.
        """
        job = str(species.chemid) + '_well'
        if high_level:
            job = str(species.chemid) + '_well_high'
        if mp2:
            job = str(species.chemid) + '_well_mp2'
        if bls:
            job = str(species.chemid) + '_well_ls'

        # TODO: Code exceptions into their own function/py script that opt can call.
        # TODO: Fix symmetry numbers for calcs as well if needed
        # O2
        if species.chemid == "320320000000000000001":
            mult = 3
        # CH2
        elif species.chemid == "140260020000000000001":
            mult = 3
        else:
            mult = species.mult

        kwargs = self.get_qc_arguments(job, mult, species.charge, high_level=high_level)

        if self.qc == 'gauss':
            kwargs['opt'] = 'CalcFC, Tight'
        if mp2:
            kwargs['method'] = self.scan_method
            kwargs['basis'] = self.scan_basis
        if high_level:
            if self.opt:
                kwargs['opt'] = 'CalcFC, {}'.format(self.opt)
        # the integral is set in the get_qc_arguments parts, bad design

        atom = copy.deepcopy(species.atom)
        dummy = geometry.is_linear(geom, species.bond)
        if len(dummy) > 0:  # add a dummy atom for each close to linear angle
            for d in dummy:
                atom = np.append(atom, ['X'])
                geom = np.concatenate((geom, [d]), axis=0)
        dummy = [d.tolist() for d in dummy]

        template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_opt_well.tpl.py'.format(qc=self.qc))
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(atom),
                                   geom=list([list(gi) for gi in geom]),
                                   ppn=self.ppn,
                                   dummy=dummy,
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd())

        f_out = open('{}.py'.format(job), 'w')
        f_out.write(template)
        f_out.close()

        self.submit_qc(job)
        return 0

    def qc_freq(self, species, geom, high_level=0):
        """
        Creates a frequency input and runs it.
        """

        job = str(species.chemid) + '_fr'
        if high_level:
            job = str(species.chemid) + '_fr_high'

        kwargs = self.get_qc_arguments(job, species.mult, species.charge, high_level=high_level)
        if self.qc == 'gauss':
            kwargs['freq'] = 'freq'
            kwargs['ioplist'] = ['7/33=1']
        elif self.qc == 'nwchem':
            kwargs['task'] = 'frequencies'

        atom = copy.deepcopy(species.atom)

        dummy = geometry.is_linear(geom, species.bond)
        if len(dummy) > 0:  # add a dummy atom for each close to linear angle
            for d in dummy:
                atom = np.append(atom, ['X'])
                geom = np.concatenate((geom, [d]), axis=0)
            # switch on the symmetry of gaussian
            if 'NoSymm' in kwargs:
                del kwargs['NoSymm']
        dummy = [d.tolist() for d in dummy]

        template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_freq_well.tpl.py'.format(qc=self.qc))
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(atom),
                                   geom=list([list(gi) for gi in geom]),
                                   ppn=self.ppn,
                                   dummy=dummy,
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd())

        f_out = open('{}.py'.format(job), 'w')
        f_out.write(template)
        f_out.close()

        self.submit_qc(job)

        return 0

    def qc_opt_ts(self, species, geom, high_level=0):
        """
        Creates a ts optimization input and runs it
        """

        job = str(species.name)
        if high_level:
            job += '_high'

        kwargs = self.get_qc_arguments(job, species.mult, species.charge, ts=1, step=1, max_step=1, high_level=1)

        template_file = pkg_resources.resource_filename('tpl', 'ase_{qc}_ts_end.tpl.py'.format(qc=self.qc))
        template = open(template_file, 'r').read()
        template = template.format(label=job,
                                   kwargs=kwargs,
                                   atom=list(species.atom),
                                   geom=list([list(gi) for gi in geom]),
                                   ppn=self.ppn,
                                   qc_command=self.qc_command,
                                   working_dir=os.getcwd())

        f_out = open('{}.py'.format(job), 'w')
        f_out.write(template)
        f_out.close()

        self.submit_qc(job)

        return 0

    def submit_qc(self, job, singlejob=1):
        """
        Submit a job to the queue, unless the job:
            * has finished with Normal termination
            * has finished with Error termination
            * is currently running
        However, if the optional parameter singlejob is set to zero, then
        the job is run only if it has finished earlier with normal termination.
        This is for continuations, when the continuing jobs overwrite each other.
        If the number of jobs in the queue is larger than the user-set limit,
        KinBot will park here until resources are freed up.
        """
        # if the logfile already exists, copy it with another name
        if self.queue_job_limit > 0:
            self.limit_jobs()

        check = self.check_qc(job)
        if singlejob == 1:
            if check != 0:
                return 0
        else:
            if check == 'running':
                return 0

        try:
            if self.par.par['queue_template'] == '':
                template_head_file = pkg_resources.resource_filename('tpl', self.queuing + '.tpl')
            else:
                template_head_file = self.par.par['queue_template']
        except OSError:
            logging.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logging.error('Or no file is found at {}.'.format(self.par.par['queue_template']))
            logging.error('Exiting')
            sys.exit()

        template_file = pkg_resources.resource_filename('tpl', self.queuing + '_python.tpl')
        python_file = '{}.py'.format(job)
        name = job.split('/')[-1]
        python_template = open(template_head_file, 'r').read() + open(template_file, 'r').read()

        if self.queuing == 'pbs':
            python_template = python_template.format(name=job, ppn=self.ppn, queue_name=self.queue_name,
                                                        dir='perm', python_file=python_file, arguments='')
        elif self.queuing == 'slurm':
            python_template = python_template.format(name=job, ppn=self.ppn, queue_name=self.queue_name, dir='perm',
                                                        slurm_feature=self.slurm_feature, python_file=python_file, arguments='')
        else:
            logging.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logging.error('Exiting')
            sys.exit()

        qu_file = '{}{}'.format(job, constants.qext[self.queuing])
        with open(qu_file, 'w') as f_out_qu:
            f_out_qu.write(python_template)

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
            logging.error(msg)
            sys.exit()
        self.job_ids[job] = pid

        return 1  # important to keep it 1, this is the natural counter of jobs submitted

    def get_qc_geom(self, job, natom, wait=0, allow_error=0, previous=0):
        """
        Get the geometry from the ase database file.
        Returns it, with the following conditions about the status of the job.
        If wait = 0, return an (1, empty array) when the job is still running (instead of the geometry).
        If wait = 1, wait for the job to finish.
        If wait = 2, return the last geometry while the job is still running.
            This option is to monitor the formation of bimolecular products.
        If allow_error = 0, do not read the final geometry if the status is not "normal"
        if allow_error = 1, read the geometry even though there is an error in the output file
            This option is to read the final IRC geometry when it did not converge
        If previous = 0, read the last geometry, this is the normal behavious
        if previous = 1, read the geometry before the last one, this is needed in scan types so
            that the max energy point is taken, not the one after that
        """
        geom = np.zeros((natom, 3))

        check = self.check_qc(job)
        if check == 'error' and not allow_error:
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
                    return 1, geom
            else:
                break
        if check != 'normal':
            if not allow_error:
                if wait != 2:
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
            found_entry = 1

        if found_entry and previous == 0:
            return status, geom
        elif found_entry and previous == 1:
            return status, prev_geom
        else:
            return -1, geom

    def get_second_to_last_geom(self, job, natom, wait=0, allow_error=0):
        """
        Get the geometry from the ase database file.
        Returns it, with the following conditions about the status of the job.
        If wait = 0, return an (1, empty array) when the job is still running (instead of the geometry).
        If wait = 1, wait for the job to finish.
        If wait = 2, return the last geometry while the job is still running.
            This option is to monitor the formation of bimolecular products.
        If allow_error = 0, do not read the final geometry if the status is not "normal"
        if allow_error = 1, read the geometry even though there is an error in the output file
            This option is to read the final IRC geometry when it did not converge
        """

        geom = np.zeros((natom, 3))

        check = self.check_qc(job)
        if check == 'error' and not allow_error:
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
                    return 1, geom
            else:
                break
        if check != 'normal':
            if not allow_error:
                if wait != 2:
                    return -1, geom

        # open the database
        rows = self.db.select(name=job)

        found_entry = 0
        geoms = []
        # take the last entry
        for row in rows:
            mol = row.toatoms()
            geoms.append(mol.positions)
            found_entry = 1

        if found_entry:
            return status, geoms[-2]
        else:
            return -1, geom

    def get_qc_freq(self, job, natom, wait=0, allow_error=0):
        """
        Get the frequencies from the ase database file
        If wait is set to 1, it will wait for the job to finish.
        """

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
        """
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
        """

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

        # open the database
        rows = self.db.select(name=job)
        energy = 0
        # take the last entry
        for row in rows:
            if hasattr(row, 'data'):
                if row.data.get('energy') is not None:
                    energy = row.data.get('energy')

        # ase energies are always in ev, convert to hartree
        energy *= constants.EVtoHARTREE

        return 0, energy

    def get_qc_zpe(self, job, wait=1):
        """
        Read the zero point energy.
        If wait is set to 1 (default), it will wait for the job to finish.
        """

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
        rows = self.db.select(name=job)
        # take the last entry
        for row in rows:
            if hasattr(row, 'data'):
                zpe = row.data.get('zpe')
            else:
                zpe = 0.0
                logging.warning("{} has no zpe in database. ZPE SET TO 0.0".format(job))

        return 0, zpe

    def read_qc_hess(self, job, natom):
        """
        Read the hessian of a gaussian chk file
        """

        check = self.check_qc(job)
        if check != 'normal':
            return []

        # initialize hessian matrix
        hess = np.zeros((3 * natom, 3 * natom))

        if self.qc == 'gauss':

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
        return hess

    def is_in_database(self, job):
        """
        Checks if the current job is in the database:
        """
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
        """
        Checks the status of the qc job.
        """
        logging.debug('Checking job {}'.format(job))
        devnull = open(os.devnull, 'w')
        if self.queuing == 'pbs':
            command = 'qstat -f | grep ' + '"Job Id: ' + self.job_ids.get(job, '-1') + '"' + ' > /dev/null'
            if int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull)) == 0:
                logging.debug('Job is running')
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
            out = out.decode()
            for line in out.split('\n'):
                if len(line) > 0:
                    while line.startswith(' '):
                        line = line[1:]
                    pid = line.split()[0]
                    if pid == self.job_ids.get(job, '-1'):
                        logging.debug('Job is running')
                        return 'running'
        else:
            logging.error('KinBot does not recognize queuing system {}.'.format(self.queuing))
            logging.error('Exiting')
            sys.exit()
        # if int(subprocess.call(command, shell=True, stdout=devnull, stderr=devnull)) == 0:
        #     return 'running'

        if self.is_in_database(job):
            for i in range(10):

                if self.qc == 'gauss':
                    log_file = job + '.log'
                elif self.qc == 'nwchem':
                    log_file = job + '.out'
                log_file_exists = os.path.exists(log_file)
                if log_file_exists:
                    logging.debug('Log file is present after {} iterations'.format(i))
                    # by deleting a log file, you allow restarting a job
                    # open the database
                    rows = self.db.select(name=job)
                    data = None
                    # take the last entry
                    for row in rows:
                        if hasattr(row, 'data'):
                            data = row.data
                    if data is None:
                        logging.debug('Data is not in database...')
                        return 0
                    else:
                        logging.debug('Returning status {}'.format(data['status']))
                        return data['status']
                else:
                    logging.debug('Checking againg for log file')
                    log_file_exists = os.path.exists(log_file)
                    time.sleep(1)

            logging.debug('log file {} does not exist'.format(log_file))
            return 0
        else:
            logging.debug('job {} is not in database'.format(job))
            return 0

    def limit_jobs(self):
        """
        Check how many jobs are in the queue from the user, and if larger than the limit,
        then wait for resources to free up.
        """

        while 1:
            if self.queuing == 'slurm':
                command = ['squeue', '-h', '-u', '{}'.format(self.username)]
            elif self.queuing == 'pbs':
                command = ['qselect', '-u', '{}'.format(self.username)]
            jobs = subprocess.check_output(command)

            if len(jobs.split(b'\n')) < self.queue_job_limit:
                return 0
            time.sleep(30)
