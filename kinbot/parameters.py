"""
File that contains all the default and user-defined parameters for the
reaction search.
The defaults are listed in this file and the user defined parameters
are defined in a json file that needs to be given as an argument to the
initializer.
"""
from __future__ import with_statement
import sys
import json
import logging


class Parameters:
    """
    This class initiates all parameters to their defaults and reads in the
    user-defined variables, which overwrite the defaults
    """
    def __init__(self, inpfile=None):
        """
        Initialize all the variable and read the file which is the user input
        file
        """
        # user input file
        self.input_file = inpfile

        self.par = {
            # User should not set this
            'input_file': self.input_file,
            # GENERAL INFO
            # title of the current calculations
            'title': '',
            # verbose log file
            'verbose': 0,

            # INPUT SPECIES INFOR
            # SMILES of the species
            'smiles': '',
            # geometry of the species
            'structure': [],
            # Charge of the species
            'charge': 0,
            # Multiplicity of the species
            'mult': 0,
            # Whether it is a dimer
            'dimer': 0,

            # WHICH STEPS TO TAKE
            # Do a reaction search
            'reaction_search': 1,
            # Which reaction families to include in the search
            'families': ['all'],
            # Which reaction families to skip in the search
            'skip_families': ['none'],
            # Which chemids to skip kinbot runs for during PES calculations
            'skip_chemids': ['none'],
            # break all single bonds to find the barriers
            # of potential homolytic scissions
            'homolytic_scissions': 0,
            # perform variational calculations for the homolytic scissions
            'variational': 0,
            # break specific bonds in the homolytic search
            # this is a dictionary written as:
            # {chemid1: [[atom1, atom2], [atom3, atom4], ...], [chemid2: [..]]}
            'barrierless_saddle': {},
            # starting distance for barrierless_saddle searches, A
            'barrierless_saddle_start': 2.0,
            # step size in A
            'barrierless_saddle_step': 0.2,
            'homolytic_bonds': [],
            # if requested with specific_reaction = 1
            # then only these bonds are broken and formed
            # atom index for break/form bonds starts at 0
            'specific_reaction': 0,
            'break_bonds': [],
            'form_bonds': [],
            # Threshold above which barriers are deemed unimportant
            'barrier_threshold': 100.,
            # Number of 0.1 Angstrom steps in bond scans
            'scan_step': 30,
            # Do a full PES scan instead of one well
            'pes': 0,
            # Maximum number of simultaneous kinbot runs in a pes search
            'simultaneous_kinbot': 5,
            # Perform high level optimization and freq calculation (L2)
            'high_level': 0,
            # Do a conformational search
            'conformer_search': 0,
            # The angular grid for dihedrals, angle = 360 / grid
            'conf_grid': 3,
            # Do a semi empirical conformational search and select the lowest conformers
            # for the L1 conformer search
            'semi_emp_conformer_search': 0,
            # Do a hindered rotor scan
            'rotor_scan': 0,
            # Number of points along the rotor scan
            'nrotation': 12,
            # Make figures of the HIR profiles
            'plot_hir_profiles': 0,
            # Do master equation calculations
            'me': 0,
            # Number of HIR restarts in case a lower energy point gets found
            'rotation_restart': 3,
            # Maximum number of diherals for which exhaustive
            # comformation searches are done
            'max_dihed': 5,
            # Number of random conformers in case no exhaustive search is done
            'random_conf': 500,
            # Maximum number of diherals for which exhaustive
            # comformation searches are done at semi empirical level
            'max_dihed_semi_emp': 5,
            # Number of random conformers in case no exhaustive search is done
            # at semi empirical level
            'random_conf_semi_emp': 500,
            # threshold of conformers at semi empirical level to take to the L1 level
            # in kcal/mol
            'semi_emp_confomer_threshold' : 5,
            # For the combinatorial search, minimum number of bonds to break
            # this value is decreased by 1 for radical reactions
            'min_bond_break': 2,
            # For the combinatorial search, maximum number of bonds to break
            # this value is decreased by 1 for radical reactions
            'max_bond_break': 3,
            # For the combinatorial search, include molecular pathways
            'comb_molec': 1,
            # For the combinatorial search, include reactions with radicals
            'comb_rad': 1,
            # For the combinatorial search, include reactions with lone electron pairs
            'comb_lone': 1,
            # For the combinatorial search, include reactions with pi electrons
            'comb_pi': 1,
            # For the combinatorial search, allow the breaking of valence
            'break_valence': 1,
            # Search for one specific reaction using combinatorial approach
            'one_reaction_comb': 0,
            # Search for one specific reaction using family approach
            'one_reaction_fam': 0,
            # For cyclic transition states, look at this ring size range
            'ringrange': [3, 9],

            # QUANTUM CHEMISTRY INFO
            # Which quantum chemistry code to use
            'qc': 'gauss',  # or nwchem
            # nwchem-specific parameter
            'methodclass': 'dft',  # or scf or mp2
            # Command for the quantum chemistry code
            'qc_command': 'g09',
            # Quantum chemistry method to use as L1
            'method': 'b3lyp',
            # Basis set to use
            'basis': '6-31G',
            # Quantum chemistry method to use as L1 for scans
            'scan_method': 'mp2',
            # Basis set to use for scans
            'scan_basis': '6-31G',
            # Method to scan bonds in barrierless_saddle family
            'barrierless_saddle_method': 'b3lyp',
            # Basis set to scan bonds in barrierless_saddle family
            'barrierless_saddle_basis': '6-31G',
            # Method to scan bonds in barrierless_saddle family
            'barrierless_saddle_method_high': 'b3lyp',
            # Basis set to scan bonds in barrierless_saddle family
            'barrierless_saddle_basis_high': '6-31G',
            # for Gaussian, request CalcAll for TS optimization
            'calcall_ts': 0,
            # Quantum chemistry method to use for high-level L2
            'high_level_method': 'M062X',
            # Basis set to use for high-level
            'high_level_basis': '6-311++G(d,p)',
            # method for semi empirical conformer search
            'semi_emp_method' : 'am1',
            # Integral grid for Gaussian, only for the high-level calculations
            'integral': '',
            # Optimization threshold
            'opt': '',
            # for Gaussian irc: IRC(MaxPoints=n)
            'irc_maxpoints': 30,
            # for Gaussian irc, IRC(StepSize=n)
            'irc_stepsize': 20,
            # for Gaussian, allow Guess=(Mix,Always)
            'guessmix': 0,
            # Turn off/on (0/1) molpro L3 calculations
            'L3_calc' : 0,
            # name of the single point code's name
            'single_point_qc': 'molpro',
            # Name of the template for the single-point calculation (L3)
            # If not specified, then the tpl/[single_point_qc].inp is used
            'single_point_template': '',
            # The keyword to be searched for in Molpro for the desired
            # energy. Compulsory if Molpro energies are used.
            'single_point_key': '',
            # L3 for barrierless template, CASPT2-like molpro is expected
            'barrierless_saddle_single_point_template': '',
            # L3 for barrierless template for product, CASPT2-like molpro is expected
            'barrierless_saddle_prod_single_point_template': '',
            # active orbitals
            'barrierless_saddle_norbital': 0,
            # number of electrons on them 
            'barrierless_saddle_nelectron': 0,
            # number of states
            'barrierless_saddle_nstate': 0,
            # single point key for barrierless
            'barrierless_saddle_single_point_key': '',
            # Command string to be used for single point energy calculation
            'single_point_command': '',
            # Hindered rotor max optimization steps
            'hir_maxcycle': None,
            # Non-rigid or rigid hir
            'rigid_hir': 0,

            # COMPUTATIONAL ENVIRONEMNT
            # Which queuing system to use
            'queuing': 'pbs',  # or slurm
            # Template for queue:
            'queue_template': '',
            # Name of the queue
            'queue_name': 'medium',
            # E.g. the type of node or anything that comes with -C in SLURM
            'slurm_feature': '',
            # Number of cores to run the L0-L2 qc jobs on
            'ppn': 1,
            # Number of cores to run the L3 qc jobs on
            'single_point_ppn': 4,
            # This many spaces can be used for numbering files, e.g., in ga
            'zf': 4,
            # delete intermediate files
            'delete_intermediate_files': 0,
            # Scratch directory
            'scratch': '',
            # User name
            'username': '',
            # Max. number of job from user in queue, if negative, ignored
            'queue_job_limit': -1,

            # MASTER EQUATION
            # Which ME code to use:
            'me_code': 'mess',  # or mesmer
            # collision parameters
            'collider': 'He',
            'epsilon': 0.0,
            'sigma': 0.0,
            # MESS specific keywords
            'mess_command': 'mess',
            'TemperatureList': [300. + 100. * i for i in range(18)],
            'PressureList': [7.6, 76., 760., 7600., 76000.],
            'EnergyStepOverTemperature': .2,
            'ExcessEnergyOverTemperature': 30,
            'ModelEnergyLimit': 400,
            'CalculationMethod': 'direct',  # or low-eigenvalue
            'ChemicalEigenvalueMax': 0.2,
            'EnergyRelaxationFactor': 200,
            'EnergyRelaxationPower': .85,
            'EnergyRelaxationExponentCutoff': 15,
            # MESMER specific keywords
            'mesmer_command': 'mesmer',

            # UQ KEYWORDS
            # UQ off/on (0/1)
            'uq': 0,
            # Number of mess input files generated for each structure
            'uq_n': 100,
            # Max number of mess calculations running at once
            'uq_max_runs': 5,
            # Uncertainty in stable intermediate  energy, +/- 0.5 kcal/mol
            'well_uq': 0.5,
            # Uncertainty in saddle point (TS) energy, +/- 1.0 kcal/mol
            'barrier_uq': 1.0,
            # Uncertainty in positive frequency values, mult/div by a maximum factor of 1.2.
            # factor of 1.2 corresponds to values ranging from 0.833 to 1.2 times the original frequency
            'freq_uq': 1.2,
            # Uncertainty in negative frequency values, mult/div by a maximum factor of 1.1.
            # factor of 1.2 corresponds to values ranging from 0.909 to 1.1 times the original frequency
            'imagfreq_uq': 1.1,

            # for development
            'test': 0,
        }

        if self.input_file is not None:
            self.read_user_input()

        if self.par['me'] == 1:
            if self.par['epsilon'] == 0. or self.par['sigma'] == 0.:
                err = 'If you want to run a ME, you need to provide sigma and epsilon for the complexes.'
                logging.error(err)
                sys.exit(-1)
            if self.par['rotor_scan'] == 0:
                err = 'If you want to run a ME, the rotor_scan needs to be turned on.'
                logging.error(err)
                sys.exit(-1)
                
        if self.par['families'] != 'all' and self.par['skip_families'] != ['none']:
            err = 'Only one of the "families" or "skip_families" parameters can be defined.'
            logging.error(err)
            sys.exit(-1)

        if self.par['pes'] and self.par['specific_reaction']:
            logging.error('Specific reaction cannot be searched in PES mode.')
            sys.exit(-1)

        if self.par['uq'] == 0:
            self.par['uq_n'] = 1

        self.par['well_uq'] = float(self.par['well_uq'])
        self.par['barrier_uq'] = float(self.par['barrier_uq'])
        self.par['freq_uq'] = float(self.par['freq_uq'])
        self.par['imagfreq_uq'] = float(self.par['imagfreq_uq'])

    def read_user_input(self):
        """
        Read the user input file and overwrite the default values
        """
        try:
            with open(self.input_file) as json_file:
                try:
                    user_data = json.load(json_file)
                except ValueError:
                    msg = 'There is an error in the json input file'
                    raise ValueError(msg)
        except IOError:
            msg = 'Input file {} does not exist'.format(self.input_file)
            raise IOError(msg)
        for key in user_data:
            if key in self.par:
                self.par[key] = user_data[key]
            else:
                err = 'KinBot does not recognize option {} with value {}'
                logging.error(err.format(key, user_data[key]))

    def print_parameters(self):
        """
        Make a string out of the parameters
        """
        s = ''
        for key in self.par:
            s += '{}\t{}\n'.format(key, self.par[key])
        return s
