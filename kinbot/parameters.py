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
from kinbot import pp_tables
from kinbot.minimum_recovery import (RECOVERY_DEFAULTS, ANL_RECOVERY_DEFAULTS,
                                     validate_anl_recovery)
from kinbot.theory import (automatic_cbh_anl_requested, profiled_requested,
                           convergence_exponent, resolve_profiles,
                           validate_anl_size_policy)

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
        # Keep explicit user choices separate from defaults.  Automatic ANL
        # can then derive subsystem limits from the ordinary KinBot queue cap
        # without overwriting an expert override.
        self.user_keys = set()

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
            # Optional per-species reaction-family/motif selection.  The
            # special key ``initial`` addresses the input well without
            # depending on its generated connectivity identifier.
            'selected_reactions': {},
            # In PES mode, retain and optimize products but do not launch a
            # recursive reaction search from newly discovered product wells.
            'stop_product_reaction_search': 0,
            # Restrict certain families to just the more likely cases
            'family_restrict': 0, 
            # Which reaction families to include in the search bimolecular reactions
            'bimol_families': ['all'],
            # Which reaction families to skip in the search
            'skip_families': ['none'],
            # Which chemids to skip kinbot runs for during PES calculations
            'skip_chemids': ['none'],
            # Saved declared configuration for generated PES inputs (including MC mirror selections).
            'stereo_reference': None,
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
            # Expert option: allow L2 to run with conformer_search off. L2
            # normally requires a conformer search so that it refines the
            # lowest conformer; with this on, L2 refines whatever conformer it
            # is given and the energies are not guaranteed to be the lowest.
            # Useful when L1 was done elsewhere, e.g. reusing geometries from a
            # different code or level. A warning is logged when it is used.
            'allow_l2_without_conf': 0,
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
            # Fixed stereoisomers are the default; racemic means the global mirror pair.
            'optical_population': 'specified',
            # Optional ground/excited electronic levels keyed by molecular
            # SMILES and charge. These states control product multiplicities
            # before QC submission and the matching MESS partition functions.
            'electronic_state_models': [],
            # Optional native JKFIT, MP2FIT, and RI/CABS identifiers.
            # Explicit model choice for unresolved harmonic/HIR coverage only.
            'optical_factor_assumptions': {},
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
            # Independent, opt-in ASE calculator profiles for L1/L2.
            'profiled_theory': 0,
            'theory_preset': '',
            'l1': '',
            'l2': '',
            'l1_profile': {},
            'l2_profile': {},
            'composite_method': '',
            # Native convergence contracts used by every automatic ANL
            # geometry, frequency, and single-point calculation.  Values are
            # exact thresholds, rather than backend-specific exponents.
            'anl_scf_convergence': 1.e-8,
            'anl_cc_convergence': 1.e-7,
            'anl_geometry_fmax': 5.e-4,
            'anl_minimum_recovery': dict(ANL_RECOVERY_DEFAULTS),
            # Optional full augmentation of the F12 reference branch.
            # Conventional ANL reference bases are always fully augmented,
            # including hydrogen.
            'anl_f12_augmented': 0,
            # Require the intended L3 geometry before recipe single points.
            'anl_require_l3_geometry': 1,
            'anl_l3_max_retries': 1,
            # Automatic CBH planning chooses the highest structurally valid
            # rung (up to CBH-3), then falls down by rung and common ANL tier
            # as required.  An integer remains accepted for reproducible
            # legacy/diagnostic runs, but ordinary users do not select it.
            'cbh_rung': 'auto',
            # Automatic composite-cost ceilings.  The values are ordinary
            # user policy, not molecule-specific exceptions; each CBH state
            # records the metrics and selected tier.
            'anl_size_policy': {
                'full_anl0_f12': {
                    'max_atoms': 16, 'max_heavy_atoms': 4,
                    'max_heavy_bonds': 4},
                'l3_geometry': {
                    'max_atoms': 30, 'max_heavy_atoms': 8,
                    'max_heavy_bonds': 8},
                'l3_single_point': {
                    'max_atoms': 50, 'max_heavy_atoms': 12,
                    'max_heavy_bonds': 12},
            },
            'atct_version': '1.220',
            'atct_cache': '~/.cache/kinbot/atct',
            'l3_overrides': {},
            'l3_resource_overrides': {},
            # Persistent profiled-theory job routing. Separate continuation
            # phases can select another ledger without deleting the original
            # run's provenance.
            'profiled_route_file': '.kinbot_theory_jobs.json',
            # A continuation may reconstruct a completed L1 reaction network
            # while permitting new L2/HIR work. Missing L1 jobs then fail
            # locally instead of being submitted again.
            'profiled_l1_reuse_only': 0,
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
            # for Gaussian irc, IRC(MaxCycle=n)
            'irc_maxcycle': 999,
            # for Gaussian, allow Guess=(Mix,Always)
            'guessmix': 0,
            # orcablocks
            'orcablocks': '',
            # orcasimpleinput
            'orcasimpleinput': '',
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
            # Turn off bad rotors based on 0 angle energies
            'rotor_0_test': 1,
            # Whether to use sella as optimizer or not
            'use_sella': False,
            # Sella hyperparameters
            'sella_kwargs': {},
            # fmax for Sella
            'sella_fmax': 0.005,
            # steps for Sella
            'sella_steps': 250,
            # Native Gaussian Hessians may expose a saddle after small forces.
            # Retry minima only, along a downhill imaginary mode; 0 disables.
            'minimum_recovery_attempts': RECOVERY_DEFAULTS['minimum_recovery_attempts'],
            # Full Cartesian displacement norm in Angstrom, doubled on retry 2.
            'minimum_recovery_displacement': RECOVERY_DEFAULTS['minimum_recovery_displacement'],
            # Continue a Sella HIR point from its last finite geometry when
            # the constrained optimization exhausts its step budget.  The
            # original calculation remains immutable evidence; recovery uses
            # a distinct job name and a larger step budget.
            'hir_recovery_attempts': 2,
            'hir_recovery_step_multiplier': 2,
            # calc_kwargs
            'calc_kwargs': {},
            # Threshold to accept negative frequencies for floppy structures, this is a positive number
            'imagfreq_threshold': 50.,
            # List of files containing the parameters for the NN model.
            'nn_model': None,
            # FairChem MLIP model configuration (used when qc = 'fc')
            # Local checkpoint path or FairChem pretrained model name.
            'fc_model_path': '',
            # Task name for FAIRChemCalculator: omol, omat, oc20, oc22, odac, omc
            'fc_task_name': 'omol',
            # Device to run on: cpu, cuda, cuda:0, etc.
            'fc_device': 'cpu',

            # VRC-TST PARAMETERS
            # Molpro memory stack in MW per MPI rank for each rotdPy sample.
            # rotdPy converts this to a total Slurm allocation for all ranks.
            'rotdPy_mem': 300,
            # Execute generated rotdPy inputs as part of the KinBot workflow.
            # Input generation remains available when this is disabled.
            'rotdpy_run': 0,
            # rotdPy submits the electronic-structure samples itself.  These
            # settings are passed to its calculator without imposing a
            # site-specific launcher in KinBot.
            'rotdpy_processors': 1,
            'rotdpy_max_jobs': 2000,
            'rotdpy_max_retries': 2,
            'rotdpy_poll_interval': 10,
            # ``partition_max`` follows the selected Slurm partition. This
            # lets a deliberately short validation partition exercise the
            # same recorded timeout/retry behavior as the ANL dispatcher.
            'rotdpy_walltime': 'partition_max',
            # Optional for sites that require a whole node for VRC samples.
            # ANL L3 jobs remain exclusive independently of this setting.
            'rotdpy_exclusive': 0,
            # [start, initial interval, interval multiplier, point count]
            'rotdpy_temperature_grid': [10., 10., 1.05, 70],
            'rotdpy_energy_grid': [0., 10., 1.05, 190],
            'rotdpy_angular_grid': [0., 1., 1.1, 80],
            'rotdpy_flux_parameters': {
                'pot_smp_max': 1000,
                'pot_smp_min': 100,
                'tot_smp_max': 15000,
                'tot_smp_min': 100,
                'flux_rel_err': 5,
                'smp_len': 1,
            },
            # Standard local recrossing correction used by VRC-TST. ROTD_py
            # applies it directly to the sampled flux, so MESS must not also
            # divide its Rotd symmetry factor by this value.
            'rotdpy_dynamical_correction': 0.85,
            # ``production`` requires a multipoint trusted correction,
            # multiple dividing surfaces, and completed Monte Carlo
            # convergence records before a ROTD result can enter MESS.
            # ``interface`` preserves inexpensive execution smoke tests.
            'rotdpy_validation_profile': 'interface',
            # Define the species and the reactions for which scans are requested
            # {chemid1: ["reaction_name1", "reaction_name2"], chemid2: [...]}
            'vrc_tst_scan': {},
            # for these, write rotdpy input, but don't do scan
            'vrc_tst_noscan': {},
            # Maximum concurrent exclusive nodes for Molpro correction jobs
            'vrc_tst_max_nodes': 1,
            # Wall time for each Molpro correction job
            'vrc_tst_walltime': '24:00:00',
            # Minimum Molpro memory stack (MW) available to every MPI rank
            'vrc_tst_min_stack_mw': 3200,
            # Explicit whole-node allocation ceiling for VRC jobs; None lets
            # site discovery choose within the available node memory.
            'vrc_tst_memory_mb': None,
            # Relaxed VRC paths use a multireference surface by default.
            # Molpro's RS2 implementation supplies analytic CASPT2 gradients
            # to Sella while the reaction-coordinate distance is constrained.
            'vrc_tst_scan_sella': 1,
            # ``legacy`` retains historical levels. ``anl1-f12`` expands one
            # compact request into an augmented-basis production VRC model.
            # The active space must be chosen explicitly; expert overrides
            # are preserved.
            'vrc_tst_accuracy_preset': 'legacy',
            # explicit preserves the historical user-supplied CAS;
            # radical_pair is a guarded CAS(2,2) shortcut for two doublets;
            # full_valence builds a whole-molecule AVAS space.
            'vrc_tst_active_space': 'explicit',
            'vrc_tst_scan_method': 'caspt2(2,2)',
            'vrc_tst_scan_basis': 'aug-cc-pVTZ',
            # auto uses an analytic gradient when the selected Molpro method
            # supplies the requested energy derivative and a numerical
            # gradient otherwise.  Expert inputs may force either contract.
            'vrc_tst_scan_gradient': 'auto',
            'vrc_tst_scan_shift': 0.1,
            'vrc_tst_scan_ipea': 0.25,
            # Optional validated AVAS/state-average model shared by the
            # relaxed scan and multireference correction energies.  This is
            # required for production work where an energy-window CAS is not
            # enough to identify the intended electronic states.
            'vrc_tst_state_model': None,
            # Optional JSON produced by molpro_reconnaissance select.  A
            # relative path is resolved beside the main KinBot input.  The
            # file is required only when the VRC branch begins, allowing the
            # same input to be prepared before the reviewed WFU exists.
            'vrc_tst_state_model_file': '',
            # Reviewed physical state contributions, separate from SA orbitals.
            'vrc_tst_surface_model_file': '',
            'vrc_tst_surface_model': None,
            'vrc_tst_selected_roots': {},
            'vrc_tst_high_extrapolation_power': 3.7,
            'vrc_tst_state_average_roots': {},
            'vrc_tst_scf_convergence': 1.e-8,
            'vrc_tst_correlation_convergence': 1.e-7,
            # Maximum constrained force (eV/angstrom) for Molpro relaxed scans.
            'vrc_tst_scan_fmax': 5.e-4,
            # Auxiliary fragment orbitals are used only for HOMO-directed
            # pivot points.  They never supply a VRC scan geometry or energy.
            'vrc_tst_pivot_backend': 'gaussian',
            'vrc_tst_pivot_method': 'ub3lyp',
            'vrc_tst_pivot_basis': '6-31+G(d)',
            # Molpro pivot backend: native OPTG and CUBE controls. fmax is
            # eV/angstrom, grid step is Bohr; method uses UKS/ or RKS/ prefix.
            'vrc_tst_pivot_fmax': 5.e-4,
            'vrc_tst_pivot_grid_step': 0.2,
            'vrc_tst_pivot_maxit': 100,
            'vrc_tst_pivot_scf_maxit': 200,
            # Calculator keywords scoped only to VRC Gaussian calculations
            'vrc_tst_calc_kwargs': None,
            # Energy calculations
            'vrc_tst_sample_method': 'caspt2(2,2)',
            'vrc_tst_high_method': 'caspt2(2,2)',
            'vrc_tst_sample_basis': 'vdz',
            'vrc_tst_sample_extrapolation_power': 3.7,
            'vrc_tst_high_basis': 'avtz',
            # ``direct`` preserves the original single trusted calculation.
            # ``spin_splitting`` reconstructs the selected lower-spin energy
            # from a high-spin RHF-UCCSD(T)-F12b anchor and configured
            # multireference spin gaps (CASPT2-F12 in the production preset).
            'vrc_tst_correction_model': 'direct',
            # Optional overrides to the documented spin-splitting preset.
            'vrc_tst_spin_splitting': None,
            # Parameters for the vrc_tst scan
            'vrc_tst_scan_points': list(np.arange(2.5, 20.0, 0.2)),
            # Separated-fragment energy reference, reviewed against recon.
            'vrc_tst_asymptote_distance': 30.,
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
            # Wall time for ordinary QC jobs. ``partition_max`` discovers the
            # selected Slurm partition's maximum, avoiding an often much
            # shorter site default without adding a required user setting.
            'queue_walltime': 'partition_max',
            # E.g. the type of node or anything that comes with -C in SLURM
            'slurm_feature': '',
            # Optional site-level exclusions. KinBot also records a node here
            # for the current run when Slurm reports an I/O launch failure
            # (exit signal 53) before the batch script starts.
            'slurm_exclude_nodes': [],
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
            # Omit optional product complexes from final MESS networks only.
            'me_skip_vdW': 0,
            # Raise submerged barriers to connected bound-well ground energies.
            'correct_submerged': 0,
            # MESS specific keywords
            # Executable and optional arguments; quote paths that contain spaces.
            'mess_command': 'mess',
            'messpf_command': 'messpf',
            # Run and parse the standalone partition-function calculations
            # after a CBH/ANL MESS model has staged them.
            'run_messpf': 0,
            # Optional schema-1 RRHO/HR model library for CBH reference
            # species that are not part of the active reaction network.
            'cbh_thermochemistry_models': '',
            # Optional output used by automatic CBH reference calculations.
            # It exports the accepted L2 geometry, modes, symmetry, and HIR
            # model without changing an ordinary KinBot run.
            'thermochemistry_model_output': '',
            # Restartable automatic CBH/ANL graph, relative to the active
            # KinBot calculation directory unless an absolute path is given.
            'anl_automatic_directory': 'anl_automatic',
            # Consume a matching, hash-verified ROTD_py number-of-states file
            # for homolytic VRC channels when it is present.  ``required``
            # refuses a phase-space fallback. ``production`` additionally
            # requires a production-validated correction and Monte Carlo run;
            # ``off`` preserves the phase-space fallback.
            'rotdpy_mess_mode': 'auto',
            # -1 selects ROTD_py's highest (most corrected) Ne_<n>.out file.
            'rotdpy_mess_energy_index': -1,
            'rotdpy_directory': 'rotdPy',
            # Optional schema-1 file that attaches completed, hash-verified
            # ANL/CBH records to the reconstructed KinBot network before MESS.
            'anl_handoff_manifest': '',
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
            # Reference frequency in cm-1 at which freq_uq is applied as is.
            # Below it the perturbation is amplified, above it dampened.
            'freq_uq_ref': 100.,
            # Cap on the amplification exponent, which bounds how hard the very
            # lowest frequencies are perturbed. The default of 4 freezes the
            # amplification below freq_uq_ref / 4 = 25 cm-1.
            'freq_uq_max_exp': 4.,
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
        # Serialize absolute checkpoint paths before PES workers change cwd.
        from pathlib import Path
        from kinbot.vrc_multistate import normalize_inline_model_paths
        base_input = self.par.get('input_file') or self.input_file
        normalize_inline_model_paths(self.par, Path(base_input).resolve().parent if base_input else Path.cwd())
        if self.par['theory_preset'] == 'b3lyp-b2plyp-anl' and 'sella_fmax' not in self.user_keys:
            self.par['sella_fmax'] = 5.e-4

        from kinbot.vrc_correction import apply_vrc_accuracy_preset
        apply_vrc_accuracy_preset(self.par, self.user_keys)

        self.par['profiled_theory'] = int(profiled_requested(self.par))
        self.theory_profiles = resolve_profiles(self.par)
        if self.theory_profiles:
            # Legacy scan defaults must not silently change a selected L1/L2
            # profile (or combine MP2 with a DFT dispersion keyword). Keep
            # explicit legacy overrides and serialize the resolved values so
            # child PES/CBH inputs retain exactly the same choices.
            for method_key, basis_key, level in (
                    ('scan_method', 'scan_basis', 'l1'),
                    ('barrierless_saddle_method', 'barrierless_saddle_basis', 'l1'),
                    ('barrierless_saddle_method_high', 'barrierless_saddle_basis_high', 'l2')):
                profile = self.theory_profiles[level]
                if profile.calculator == 'gaussian':
                    if method_key not in self.user_keys:
                        self.par[method_key] = profile.method
                    if basis_key not in self.user_keys:
                        self.par[basis_key] = profile.basis

        # One queue limit controls the complete automatic calculation.  This
        # includes KinBot QC jobs, VRC correction jobs, ROTD_py samples, and
        # the composite DAG.  Existing nonautomatic inputs and explicit
        # subsystem overrides retain their historical behavior.
        queue_limit = self.par['queue_job_limit']
        if (automatic_cbh_anl_requested(self.par)
                and not isinstance(queue_limit, bool)
                and isinstance(queue_limit, int) and queue_limit > 0):
            if 'vrc_tst_max_nodes' not in self.user_keys:
                self.par['vrc_tst_max_nodes'] = queue_limit
            if 'rotdpy_max_jobs' not in self.user_keys:
                self.par['rotdpy_max_jobs'] = queue_limit

        err = None
        if self.par['me'] == 1:
            if self.par['calc_aie'] and not self.par['conformer_search']:
                err = 'AIE calculation requires a conformer search.'
            if self.par['epsilon'] == 0. or self.par['sigma'] == 0.:
                err = 'If you want to run a ME, you need to provide sigma and epsilon for the complexes.'
            if self.par['rotor_scan'] == 0 and self.par['multi_conf_tst'] == 0:
                err = 'If you want to run a ME, the rotor_scan needs to be turned on.'
        if self.par['me'] in (1, 2):
            # convert to cm-1 units
            if self.par['epsilon_unit'] == 'cm-1':  # units in MESS
                pass
            elif self.par['epsilon_unit'] == 'K':  # CHEMKIN units, default
                self.par['epsilon'] *= units.kB / units.invcm
            elif self.par['epsilon_unit'] == 'J/mol':  # units in RMG
                self.par['epsilon'] *= units.J / units.mol / units.invcm
            else:
                err = 'Unknown unit for epsilon, has to be K, J/mol or cm-1'
            if err is None:
                # PES writes resolved parameters to worker inputs. Retain the
                # resolved unit so a worker does not convert epsilon twice.
                self.par['epsilon_unit'] = 'cm-1'
                
        if self.par['families'] != ['all'] and self.par['skip_families'] != ['none']:
            err = 'Only one of the "families" or "skip_families" parameters can be defined.'

        homolytic = self.par['homolytic_bonds']
        if not isinstance(homolytic, dict):
            raise ValueError('homolytic_bonds must be an object.')
        for key, bonds in homolytic.items():
            if (not isinstance(key, str) or not key.strip() or not isinstance(bonds, list)
                    or any(not isinstance(bond, list) or len(bond) != 2
                           or any(isinstance(atom, bool) or not isinstance(atom, int)
                                  or atom < 0 for atom in bond)
                           or bond[0] == bond[1] for bond in bonds)):
                raise ValueError('homolytic_bonds entries require zero-based atom pairs.')
            if len({tuple(sorted(bond)) for bond in bonds}) != len(bonds):
                raise ValueError('homolytic_bonds contains a duplicate bond.')

        selected = self.par['selected_reactions']
        if not isinstance(selected, dict):
            raise ValueError('selected_reactions must be an object.')
        for species_key, families in selected.items():
            if (not isinstance(species_key, str) or not species_key.strip()
                    or not isinstance(families, dict) or not families):
                raise ValueError(
                    'selected_reactions entries must map a species name to '
                    'a nonempty family object.')
            for family, motifs in families.items():
                if (not isinstance(family, str) or not family.strip()
                        or not isinstance(motifs, list) or not motifs):
                    raise ValueError(
                        'selected_reactions family entries must contain a '
                        'nonempty list of atom motifs.')
                normalized = []
                for motif in motifs:
                    if (not isinstance(motif, list) or not motif
                            or any(isinstance(atom, bool)
                                   or not isinstance(atom, int) or atom < 0
                                   for atom in motif)):
                        raise ValueError(
                            'selected_reactions motifs must be nonempty '
                            'zero-based lists of nonnegative integers.')
                    value = tuple(motif)
                    if value in normalized:
                        raise ValueError(
                            'selected_reactions contains a duplicate motif.')
                    normalized.append(value)
        if (isinstance(self.par['stop_product_reaction_search'], bool)
                or self.par['stop_product_reaction_search'] not in (0, 1)):
            raise ValueError('stop_product_reaction_search must be 0 or 1.')
        cbh_rung = self.par['cbh_rung']
        if not (cbh_rung == 'auto'
                or (not isinstance(cbh_rung, bool)
                    and isinstance(cbh_rung, int) and cbh_rung in range(4))):
            raise ValueError('cbh_rung must be auto, 0, 1, 2, or 3.')
        if cbh_rung != 'auto' and cbh_rung and not self.par['composite_method']:
            raise ValueError('an explicit CBH rung above zero requires '
                             'composite_method.')
        self.par['anl_size_policy'] = validate_anl_size_policy(
            self.par['anl_size_policy'])
        convergence_exponent(
            self.par['anl_scf_convergence'], 'anl_scf_convergence')
        convergence_exponent(
            self.par['anl_cc_convergence'], 'anl_cc_convergence')
        for key in ('sella_fmax', 'anl_geometry_fmax', 'vrc_tst_scan_fmax'):
            value = self.par[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value) or value <= 0:
                raise ValueError(f'{key} must be a positive finite force threshold in eV/angstrom.')
        self.par['anl_minimum_recovery'] = validate_anl_recovery(
            self.par['anl_minimum_recovery'])
        if (isinstance(self.par['anl_f12_augmented'], bool)
                or self.par['anl_f12_augmented'] not in (0, 1)):
            raise ValueError('anl_f12_augmented must be 0 or 1.')
        if (isinstance(self.par['anl_require_l3_geometry'], bool)
                or self.par['anl_require_l3_geometry'] not in (0, 1)):
            raise ValueError('anl_require_l3_geometry must be 0 or 1.')
        if type(self.par['anl_l3_max_retries']) is not int or self.par['anl_l3_max_retries'] < 0:
            raise ValueError('anl_l3_max_retries must be a nonnegative integer.')
        if (not isinstance(self.par['atct_version'], str)
                or not self.par['atct_version'].strip()):
            raise ValueError('atct_version must be a nonempty string.')
        if (not isinstance(self.par['atct_cache'], str)
                or not self.par['atct_cache'].strip()):
            raise ValueError('atct_cache must be a nonempty path string.')

        if self.par['pes'] and self.par['specific_reaction']:
            err = 'Specific reaction cannot be searched in PES mode.'

        if self.par['high_level'] == 1 and self.par['conformer_search'] == 0:
            if not self.par['allow_l2_without_conf']:
                err = 'Conformer search has to be done before L2.'
            elif show_warnings:
                logger.warning('Running L2 without a conformer search, '
                               'because allow_l2_without_conf is on. L2 will '
                               'refine the conformer it is given, and the '
                               'energies are not guaranteed to belong to the '
                               'lowest conformer.')

        if (self.par['high_level'] == 0 and self.par['rotor_scan'] == 1
                and self.par['qc'].lower() != 'fc'
                and not self.par['profiled_theory']):
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
        
        if self.par['qc'] == 'fc' and not self.par['fc_model_path']:
            err = 'fc_model_path must be set when using qc = "fc".'

        profiled_sella = any(
            profile.optimizer == 'sella'
            or profile.calculator in ('fairchem', 'orca')
            for profile in self.theory_profiles.values())
        if (self.par['use_sella'] == 1 or self.par['qc'] in ('fc', 'orca')
                or profiled_sella):
            if 'internal' not in self.par['sella_kwargs']:
                self.par['sella_kwargs']['internal'] = True
                logger.warning('Automatically turning on internal coordiantes for Sella.')
            if 'gamma' not in self.par['sella_kwargs'] and self.par['calcall_ts']:
                self.par['sella_kwargs']['gamma'] = 0.0
                logger.warning('Automatically setting gamma=0.0 in Sella to mimic CalcAll.')

        if self.par['uq'] == 0:
            self.par['uq_n'] = 1

        if (self.par['bimol'] == 1 and self.par['method'] == 'b3lyp'
                and self.par['qc'].lower() != 'fc'):
            logger.warning('B3LYP is not recommended as L1 for bimolecular reactions.')
            logger.warning('Choose for instance M06-2X or wB97XD.')
            print('B3LYP is not recommended as L1 for bimolecular reactions.')
            print('Choose for instance M06-2X.')

        if self.par['bimol'] and len(self.par['structure']) != 2:
            err = 'For bimolecular reactions two fragments need to be defined.'

        if self.par['optical_population'] not in ('specified', 'racemic'):
            raise ValueError('optical_population must be specified or racemic')
        from kinbot.electronic_states import electronic_state_models
        electronic_state_models(self.par)
        if not isinstance(self.par['vrc_tst_surface_model_file'], str):
            raise ValueError('vrc_tst_surface_model_file must be a path string')
        if self.par['vrc_tst_surface_model'] is not None and self.par['vrc_tst_surface_model_file']:
            raise ValueError('Specify an inline surface model or a surface-model file, not both.')
        for key in ('vrc_tst_selected_roots', 'vrc_tst_state_average_roots'):
            from kinbot.vrc_correction import root_overrides
            root_overrides(self.par, key)
        if self.par['me_skip_vdW'] not in (0, 1):
            raise ValueError('me_skip_vdW must be 0 or 1')
        if self.par['run_messpf'] not in (0, 1):
            raise ValueError('run_messpf must be 0 or 1')
        if self.par['run_messpf'] and self.par['me'] not in (1, 2):
            raise ValueError('run_messpf requires me = 1 or me = 2')
        if not isinstance(self.par['cbh_thermochemistry_models'], str):
            raise ValueError('cbh_thermochemistry_models must be a path string')
        if self.par['rotdpy_mess_mode'] not in (
                'auto', 'required', 'production', 'off'):
            raise ValueError('rotdpy_mess_mode must be auto, required, '
                             'production, or off')
        from kinbot.rotdpy import validate_mess_handoff
        validate_mess_handoff(self.par)
        if self.par['rotdpy_validation_profile'] not in (
                'interface', 'production'):
            raise ValueError('rotdpy_validation_profile must be interface '
                             'or production')
        from kinbot.vrc_correction import (correction_model,
                                           spin_splitting_settings,
                                           validate_vrc_levels)
        from kinbot.molpro import validate_vrc_state_model
        state_model = self.par['vrc_tst_state_model']
        if not (isinstance(state_model, dict) and state_model.get('schema') == 2):
            validate_vrc_state_model(state_model)
        validate_vrc_levels(self.par)
        if correction_model(self.par) == 'spin_splitting':
            spin_splitting_settings(self.par)
        if self.par['vrc_tst_scan_gradient'] not in (
                'auto', 'analytic', 'numerical'):
            raise ValueError('vrc_tst_scan_gradient must be auto, analytic, '
                             'or numerical')
        for key in ('vrc_tst_scan_shift', 'vrc_tst_scan_ipea'):
            value = self.par[key]
            if (isinstance(value, bool) or not isinstance(value, (int, float))
                    or not np.isfinite(value) or value < 0.):
                raise ValueError(f'{key} must be a finite nonnegative number')
        convergence_exponent(
            self.par['vrc_tst_scf_convergence'],
            'vrc_tst_scf_convergence')
        from kinbot.vrc_orbitals import pivot_backend, molpro_pivot_settings
        if pivot_backend(self.par) == 'molpro':
            molpro_pivot_settings(self.par)
        convergence_exponent(
            self.par['vrc_tst_correlation_convergence'],
            'vrc_tst_correlation_convergence')
        energy_index = self.par['rotdpy_mess_energy_index']
        if (isinstance(energy_index, bool)
                or not isinstance(energy_index, int) or energy_index < -1):
            raise ValueError(
                'rotdpy_mess_energy_index must be -1 or a nonnegative integer')
        if (not isinstance(self.par['rotdpy_directory'], str)
                or not self.par['rotdpy_directory'].strip()):
            raise ValueError('rotdpy_directory must be a nonempty path string')
        if not isinstance(self.par['anl_handoff_manifest'], str):
            raise ValueError('anl_handoff_manifest must be a path string')
        if not isinstance(self.par['thermochemistry_model_output'], str):
            raise ValueError('thermochemistry_model_output must be a path string')
        if (not isinstance(self.par['anl_automatic_directory'], str)
                or not self.par['anl_automatic_directory'].strip()):
            raise ValueError('anl_automatic_directory must be a nonempty path string')
        assumptions = self.par['optical_factor_assumptions']
        if not isinstance(assumptions, dict) or any(
                not isinstance(name, str) or not name or not isinstance(value, dict)
                or set(value) != {'multiplier', 'reason'}
                or isinstance(value.get('multiplier'), bool) or value.get('multiplier') not in (1, 2)
                or not isinstance(value.get('reason'), str) or not value['reason'].strip()
                or '\n' in value['reason'] or '\r' in value['reason']
                for name, value in assumptions.items()):
            raise ValueError('optical_factor_assumptions requires exact state names with multiplier 1 or 2 and a one-line reason')

        if self.par['multi_conf_tst'] and not self.par['conformer_search']:
            err = 'For multi conformer tst calculation conformer search needs to be activated.'

        attempts = self.par['minimum_recovery_attempts']
        if isinstance(attempts, bool) or not isinstance(attempts, int) or attempts < 0:
            raise ValueError('minimum_recovery_attempts must be a nonnegative integer.')
        displacement = self.par['minimum_recovery_displacement']
        if (isinstance(displacement, bool) or not isinstance(displacement, (int, float))
                or not np.isfinite(displacement) or displacement <= 0):
            raise ValueError('minimum_recovery_displacement must be finite and positive.')

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
        from kinbot.vrc_correction import asymptote_distance
        asymptote_distance(self.par)
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
        self.user_keys = set(user_data)

        return

    def print_parameters(self):
        """
        Make a string out of the parameters
        """
        s = ''
        for key in self.par:
            s += '{}\t{}\n'.format(key, self.par[key])
        return s
