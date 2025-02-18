"""
File that contains all the default and user-defined parameters for the
reaction search.
The defaults are listed in this file and the user defined parameters
are defined in a json file that needs to be given as an argument to the
initializer.
"""
import sys
import json
import logging
import numpy as np
from ase import units
from kinbot import kb_path
from kinbot import pp_tables
from kinbot import constants

logger = logging.getLogger('KinBot')


class Parameters:
    """
    This class initiates all parameters to their defaults and reads in the
    user-defined variables, which overwrite the defaults.
    """
    def __init__(self, inpfile=None, show_warnings=False):
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

            # INPUT SPECIES INFO
            # SMILES of the species
            'smiles': '',
            # geometry of the species
            'structure': [],
            # Charge of the species
            'charge': 0,
            # Multiplicity of the species
            'mult': 0,
            # Whether it is a bimolecular reaction
            'bimol': 0,
            # Cluster
            'cluster' : 0,
            # Atom numbers of the solute, only used in cluster mode
            'solute': [],
            # H-bond recognition as bond
            'hbond': 0,

            # WHICH STEPS TO TAKE
            # Do a reaction search
            'reaction_search': 1,
            # Which reaction families to include in the search
            'families': ['all'],
            # Which reaction families to include in the search bimolecular reactions
            'bimol_families': ['all'],
            # Which reaction families to skip in the search
            'skip_families': ['none'],
            # Which chemids to skip kinbot runs for during PES calculations
            'skip_chemids': ['none'],
            # Which chemids to keep kinbot runs for during PES calculations
            'keep_chemids': ['none'],
            # Skip specific reactions, usually makes sense once the search is done
            'skip_reactions': ['none'],
            # perform variational calculations for the homolytic scissions
            'variational': 0,
            # break specific bonds in the homolytic search
            # this is a dictionary written as:
            # {chemid1: [[atom1, atom2], [atom3, atom4], ...], [chemid2: [..]]}
            'barrierless_saddle': {},
            # starting distance for barrierless_saddle searches in A
            'barrierless_saddle_start': 2.0,
            # step size in A 
            'barrierless_saddle_step': 0.2,
            # for the hom_sci family, using the same format as in barrierless_saddle
            'homolytic_bonds': {},
            # if requested with specific_reaction = 1
            # then only these bonds are broken and formed
            # atom index for break/form bonds starts at 0
            'specific_reaction': 0,
            'break_bonds': [],
            'form_bonds': [],
            # Threshold above which barriers are deemed unimportant at L1
            'barrier_threshold': None,
            # Additional threshold if L2 screening is desired instead of L1
            'barrier_threshold_L2': None,
            # Threshold allowance for L1 relative to L2
            'barrier_threshold_add': 10.,
            # Additional barrier allowance for homolytic scissions
            'hom_sci_threshold_add': 5.,
            # Number of 0.1 Angstrom steps in bond scans
            'scan_step': 30,
            # Block certain reactions beyond this distace in A 
            'rigid_reaction_cutoff': 3.0,
            # Do a full PES scan instead of one well
            'pes': 0,
            # Maximum number of simultaneous kinbot runs in a pes search
            'simultaneous_kinbot': 5,
            # Perform high level optimization and freq calculation (L2)
            'high_level': 0,
            # Calculate AIE for each conformer - requires conformer search
            'calc_aie': 0,
            # Detect vdW wells deeper than threshold (kcal/mol)
            'vdW_detection': 0.5,
            

            # CONFORMATIONAL SEARCH
            # Do a conformational search
            'conformer_search': 0,
            # Threshold to differentiate two structures, kcal/mol
            'difference_threshold': 0.1,
            # The angular grid for dihedrals, angle = 360 / grid
            'conf_grid': 3,
            # Do a semi empirical conformational search and select the lowest conformers
            # for the L1 conformer search
            'semi_emp_conformer_search': 0,
            # Do a hindered rotor scan
            'rotor_scan': 0,
            # Energy threshold for free rotor in kcal/mol
            'free_rotor_thrs': 0.1,
            # Number of points along the rotor scan
            'nrotation': 12,
            # Make figures of the HIR profiles
            'plot_hir_profiles': 0,
            # Number of HIR restarts in case a lower energy point gets found
            'rotation_restart': 3,
            # Max. threshold in kcal/mol that is allowed for a HIR restart to be triggered
            # To avoid invalid rotors in restart
            'rotation_restart_threshold': 5.0,
            # Maximum number of diherals for which exhaustive
            # comformation searches are done
            'max_dihed': 5,
            # Number of random conformers in case no exhaustive search is done
            'random_conf': 500,
            # Dihedral angle to consider a section of a ring flat
            'flat_ring_dih_angle': 5.,
            # Maximum number of diherals for which exhaustive
            # comformation searches are done at semi empirical level
            'max_dihed_semi_emp': 5,
            # Number of random conformers in case no exhaustive search is done
            # at semi empirical level
            'random_conf_semi_emp': 500,
            # threshold of conformers at semi empirical level to take to the L1 level
            # in kcal/mol
            'semi_emp_confomer_threshold': 5,
            # multi conformer TST
            'multi_conf_tst': 0,
            # temperature in K
            'multi_conf_tst_temp': 300.0,
            # percent of Boltzmann to include
            'multi_conf_tst_boltz': 0.05, 
            # force print conformational info, might be slow 
            'print_conf': 0,

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

            # CALCULATION PARAMETERS
            # Which quantum chemistry code to use
            'qc': 'gauss',  # or nwchem or nn_pes
            # nwchem-specific parameter
            'methodclass': 'dft',  # or scf or mp2
            # Command for the quantum chemistry code
            'qc_command': 'g16',
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
            'calcall_ts': 1,
            # Quantum chemistry method to use for high-level L2
            'high_level_method': 'M062X',
            # Basis set to use for high-level
            'high_level_basis': '6-311++G(d,p)',
            # method for semi empirical conformer search
            'semi_emp_method': 'am1',
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
            'L3_calc': 0,
            # name of the single point code's name
            'single_point_qc': 'molpro',
            # Name of the template for the single-point calculation (L3)
            # If not specified, then the tpl/[single_point_qc].inp is used
            'single_point_template': '',
            # The keyword to be searched for in Molpro for the desired
            # energy. Compulsory if Molpro energies are used.
            'single_point_key': 'mytza',
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
            # Whether to use sella as optimizer or not
            'use_sella': False,
            # Sella hyperparameters
            'sella_kwargs': {},
            # calc_kwargs
            'calc_kwargs': {},
            # Threshold to accept negative frequencies for floppy structures, this is a positive number
            'imagfreq_threshold': 50.,
            # List of files containing the parameters for the NN model. 
            'nn_model': None,

            # VRC-TST PARAMETERS
            # Amount (Mb) of memory to use in rotdPy for each job during the sampling.
            'rotdPy_mem': 300,
            # Define the species and the reactions for which scans are requested
            # {chemid1: ["reaction_name1", "reaction_name2"], chemid2: [...]}
            'vrc_tst_scan': {},
            # for these, write rotdpy input, but don't do scan
            'vrc_tst_noscan': {},
            # using sella for scan
            'vrc_tst_scan_sella': 0,
            # Method to scan bonds in vrc_tst_scan
            'vrc_tst_scan_method': 'ub3lyp',
            # Basis set to scan bonds in vrc_tst_scan
            'vrc_tst_scan_basis': '6-31+G(d)',
            # Energy calculations
            'vrc_tst_sample_method': 'caspt2(2,2)',
            'vrc_tst_high_method': 'caspt2(2,2)',
            'vrc_tst_sample_basis': 'vdz',
            'vrc_tst_high_basis': 'avtz',
            # Parameters for the vrc_tst scan
            'vrc_tst_scan_points': list(np.arange(2.5, 20.0, 0.2)),
            'vrc_tst_scan_molpro_key': 'MYENERGY',
            # Must be provided
            'vrc_tst_scan_molpro_tpl': '',
            # Max. rmsd deviation allowed
            'vrc_tst_scan_deviation': 100.,
            # Max angle deviation allowed
            'vts_ang_dev': 10,
            # Explicit reaction center for a fragment, {'frament chemid': [atomids of centers]}
            'vrc_tst_scan_reac_cent': {},
            # Dictionary of distances in bohr at which pivot points are generated for each atom
            'pp_length': None,
            # List [start, stop] in angstrom of the pp_oriented procedure
            'pp_oriented': None,
            # List [start, stop] in angstrom of the pp_on_atom procedure
            'pp_on_atom': None,
            # Start value in angstrom of the pp_on_COM procedure
            'pp_on_COM': 10.0,
            # mode of pivot point placement. Accepted values are geometric and homo
            'pp_orient': 'homo',
            # list of distances used for the surfaces
            'rotdpy_dist': list(np.arange(3, 20.0, 0.2)),

            # COMPUTATIONAL ENVIRONEMNT
            # Which queuing system to use
            'queuing': 'pbs',  # or slurm
            # Template for queue:
            'queue_template': '',
            # Queue template for AM1 jobs
            'q_temp_am1': None,
            # Queue template for high jobs
            'q_temp_hi': None, 
            # Queue template for MP2 jobs
            'q_temp_mp2': None,
            # Queue template for MP2 jobs
            'q_temp_l3': None,
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
            # Whether to raise an error when 'queuing' is set to 'local' and the 
            # files and db entries are missing, otherwise just show a warning.
            'error_missing_local': True,
            # Whether to perform the initial cleanup of files.
            'do_clean': True,

            # MASTER EQUATION
            # Assemble the ME
            'me': 0,
            # Run the ME
            'run_me': 0,
            # Which ME code to use:
            'me_code': 'mess',  # or mesmer
            # collision parameters
            'collider': 'He',
            'epsilon': 0.0,
            'epsilon_unit': 'K',  # can be K or J/mol or cm-1
            'sigma': 0.0,
            # correct submerged barrier
            'correct_submerged': 0,
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
            # factor of 1.1 corresponds to values ranging from 0.909 to 1.1 times the original frequency
            'imagfreq_uq': 1.1,
            # LJ parameters
            'epsilon_uq': 1.2,
            'sigma_uq': 1.2,
            # Collisional parameters
            'enrelfact_uq': 1.2,
            'enrelpow_uq': 1.2,
            # PST symmetry number
            'pstsymm_uq': 2.0,
        }

        if self.input_file is not None:
            self.read_user_input()

        err = None
        if self.par['me'] == 1:
            if self.par['calc_aie'] and not self.par['conformer_search']:
                err = 'AIE calculation requires a conformer search.'
            if self.par['epsilon'] == 0. or self.par['sigma'] == 0.:
                err = 'If you want to run a ME, you need to provide sigma and epsilon for the complexes.'
            if self.par['rotor_scan'] == 0 and self.par['multi_conf_tst'] == 0:
                err = 'If you want to run a ME, the rotor_scan needs to be turned on.'
            # convert to cm-1 units
            if self.par['epsilon_unit'] == 'cm-1':  # units in MESS
                pass
            elif self.par['epsilon_unit'] == 'K':  # CHEMKIN units, default
                self.par['epsilon'] *= units.kB / units.invcm
            elif self.par['epsilon_unit'] == 'J/mol':  # units in RMG
                self.par['epsilon'] *= units.J / units.mol / units.invcm
            else:
                err = 'Unknown unit for epsilon, has to be K, J/mol or cm-1'
                
        if self.par['families'] != ['all'] and self.par['skip_families'] != ['none']:
            err = 'Only one of the "families" or "skip_families" parameters can be defined.'

        if self.par['pes'] and self.par['specific_reaction']:
            err = 'Specific reaction cannot be searched in PES mode.'

        if self.par['high_level'] == 1 and self.par['conformer_search'] == 0:
            err = 'Conformer search has to be done before L2.'

        if self.par['high_level'] == 0 and self.par['rotor_scan'] == 1:
            if show_warnings:
                logger.warning('L1 level of theory (here set to '
                               f'{self.par["method"].upper()}/{self.par["basis"].upper()}) is '
                               'designed to be used for exploratory purposes only. '
                               'Running the hindered rotor calculations at this '
                               'level of theory is highly discouraged.')
            if self.par['method'].lower() != self.par['high_level_method'].lower() \
                    or self.par['basis'].lower() != self.par['high_level_basis'].lower():
                err = 'When running the hindered rotors at L1, "high_level_method" ' \
                      'must be the same as "method" and "high_level_basis" ' \
                      'must be the same as "basis".'

        if self.par['uq'] == 0:
            self.par['uq_n'] = 1

        if self.par['bimol'] == 1 and self.par['method'] == 'b3lyp':
            logger.warning('B3LYP is not recommended as L1 for bimolecular reactions.')
            logger.warning('Choose for instance M06-2X or wB97XD.')
            print('B3LYP is not recommended as L1 for bimolecular reactions.')
            print('Choose for instance M06-2X.')

        if self.par['bimol'] and len(self.par['structure']) != 2:
            err = 'For bimolecular reactions two fragments need to be defined.'

        if self.par['multi_conf_tst'] and not self.par['conformer_search']:
            err = 'For multi conformer tst calculation conformer search needs to be activated.'

        if self.par['imagfreq_threshold'] < 0:
            err = 'The threshold for imaginary freqency has to be a positive number, interpreted as an imaginary value.'

        if not self.par['multi_conf_tst']:
            self.par['multi_conf_tst_temp'] = None
            self.par['multi_conf_tst_boltz'] = 0.05

        if self.par['multi_conf_tst']:
            self.par['rotor_scan'] = 0

        for par in ['q_temp_am1', 'q_temp_hi', 'q_temp_mp2', 'q_temp_l3']:
            if self.par[par] is None:
                self.par[par] = self.par['queue_template']

        self.par['well_uq'] = float(self.par['well_uq'])
        self.par['barrier_uq'] = float(self.par['barrier_uq'])
        self.par['freq_uq'] = float(self.par['freq_uq'])
        self.par['imagfreq_uq'] = float(self.par['imagfreq_uq'])

        # Check user input
        if self.par['pp_length'] is not None \
                and not isinstance(self.par['pp_length'], dict):
            err = 'User defined pp_length should be a dict. Using default values.'
            self.par['pp_length'] = pp_tables.pp_length_table()
        # Check keys
        elif self.par['pp_length'] is not None \
                and isinstance(self.par['pp_length'], dict):
            i = 0
            X_is_defined = False
            for key, value in self.par['pp_length'].items():
                if i == 0:
                    length = len(value)
                    i += 1
                if key == 'X':
                    X_is_defined = True
                if len(value) != length:
                    err = 'All lists in pp_length should have the same length.'
            if not X_is_defined:
                err = 'X (default) should be included in pp_length if the user defines a list for any atom.'
        elif self.par['pp_length'] is None:
            self.par['pp_length'] = pp_tables.pp_length_table()

        if self.par['pp_oriented'] is not None and\
           not isinstance(self.par['pp_oriented'], list):
            err = 'User defined pp_oriented should be a list. Using default values.'
            self.par['pp_oriented'] = [1.5, 6]
        elif self.par['pp_oriented'] is None:
            self.par['pp_oriented'] = [1.5, 6]
        if self.par['pp_on_atom'] is not None and\
           not isinstance(self.par['pp_on_atom'], list):
            err = 'User defined pp_on_atom should be a list. Using default values.'
            self.par['pp_on_atom'] = [5.0, 12.0]
        elif self.par['pp_on_atom'] is None:
            self.par['pp_on_atom'] = [5.0, 12.0]
        if self.par['pp_orient'].casefold() not in ['homo', 'geometric']:
            err = "pp_orient not well defined. 'homo' and 'geometric' are currently the only possible options"

        if self.par['barrier_threshold'] == 'none':
            self.par['barrier_threshold'] = None
        if self.par['barrier_threshold_L2'] == 'none':
            self.par['barrier_threshold_L2'] = None
        if not self.par['barrier_threshold'] and not self.par['barrier_threshold_L2']:
            err = 'One of barrier_threshold or barrier_threshold_L2 needs to be set.'
        elif self.par['barrier_threshold'] and self.par['barrier_threshold_L2']:
            logger.warning('L1 threshold is overwritten.')
        elif self.par['barrier_threshold_L2']:
            self.par['barrier_threshold'] = self.par['barrier_threshold_L2'] + self.par['barrier_threshold_add']

        try:
            if isinstance(self.par['vrc_tst_scan_points'][0], list):
                # throw an error if incorrect user input
                self.par['vrc_tst_scan_points'][0][0]
                tmp = []
                tmp_dist = []
                for sp in self.par['vrc_tst_scan_points']:
                    tmp.extend(
                        np.round(np.arange(sp[0], sp[1], sp[2]), 3).tolist())
                # avoids duplicates in the list
                for i in tmp:
                    if i not in tmp_dist:
                        tmp_dist.append(i)
                self.par['vrc_tst_scan_points'] = tmp_dist
        except (TypeError, IndexError):
            pass
        try:
            if isinstance(self.par['rotdpy_dist'][0], list):
                self.par['rotdpy_dist'][0][0]
                tmp = []
                tmp_dist = []
                for sp in self.par['rotdpy_dist']:
                    tmp.extend(
                        np.round(np.arange(sp[0], sp[1], sp[2]), 3).tolist())
                # avoids duplicates in the list
                for i in tmp:
                    if i not in tmp_dist:
                        tmp_dist.append(i)
                self.par['rotdpy_dist'] = tmp_dist
        except (TypeError, IndexError):
            pass

        if err is not None:
            logger.error(err)
            sys.exit(-1)

        return

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
                msg = f'KinBot does not recognize option {key} with value {user_data[key]}'
                raise IOError(msg)

        return

    def print_parameters(self):
        """
        Make a string out of the parameters
        """
        s = ''
        for key in self.par:
            s += '{}\t{}\n'.format(key, self.par[key])
        return s
