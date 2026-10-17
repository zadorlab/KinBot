# Production input parameters and provenance

Baseline: `origin/master` as recorded at commit `5c6525a7` (5 October 2026). “Branch” includes earlier work on `composite`, not only the latest fix. Existing backend algorithms are distinguished from KinBot orchestration. Values are copied from the actual inputs; the complete JSON remains authoritative.

| Parameter | CH3S / CH3SS value | Hook origin and purpose |
| --- | --- | --- |
| `title` | `"ch3s_anl1_f12"` / `"ch3ss_anl1_f12"` | Master key. Run label; no change to chemistry. |
| `smiles` | `"[CH3][S]"` / `"C[S][S]"` | Master key. Initial molecular connectivity; explicit hydrogens are generated in the recorded atom order. |
| `charge` | `0` | Master key. Neutral total charge. |
| `mult` | `2` | Master key. Doublet parent; fragment multiplicities are assigned separately. |
| `barrier_threshold_L2` | `400.0` | Master key. ZPE-corrected L2 barrier ceiling in kcal/mol; a conservative screen for the requested 5000 K maximum. |
| `barrier_threshold_add` | `10.0` | Master key. L1 screening allowance above the L2 ceiling: resolved L1 limit is 410 kcal/mol. |
| `TemperatureList` | `300,400,...,5000` | Master MESS/PF temperature list in kelvin. |
| `rotdpy_temperature_grid` | `[300,100,1,48]` | Existing RotdPy grid hook: start, initial interval, interval multiplier, count. Matches the MESS temperature list. |
| `rotdpy_energy_grid` | `[0,10,1.05,190]` | Existing grid hook in kelvin; retains the usual energy grid, ending near 4018 kcal/mol. |
| `ModelEnergyLimit` | `800` | Master MESS partition-function/model energy ceiling, kcal/mol; separate from the reaction-search cutoff. |
| `ExcessEnergyOverTemperature` | `30` | Master MESS control; relaxation matrix extends 30 RT above the highest model barrier. |
| `conformer_search` | `1` | Master key. Search conformations before final thermochemistry. |
| `high_level` | `1` | Master key. Refine accepted structures at the L2 profile. |
| `rotor_scan` | `1` | Master key. Build hindered-rotor potentials where chemically applicable. |
| `theory_preset` | `"b3lyp-b2plyp-anl"` | Branch hook. L1 Gaussian B3LYP-D3(BJ)/cc-pVDZ exploration; L2 Gaussian B2PLYP-D3(BJ)/cc-pVTZ. Both use Sella optimization and native harmonic Hessians. Open-shell forms follow the supplied multiplicity. |
| `composite_method` | `"ANL1-F12"` | Branch hook. Automatic CBH/ATcT plus ANL1-F12 component workflow for accepted species and stationary points; not a CASPT2 kinetic sampling method. |
| `anl_f12_augmented` | `1` | Branch hook. Use fully augmented QZ-F12/5Z-F12 reference bases for ANL1-F12 thermochemistry. This is distinct from the VRC aDZ/aTZ pair. |
| `vrc_tst_accuracy_preset` | `"anl1-f12"` | Branch hook. CASPT2-only VRC defaults; explicit values below override the preset. It does not select a chemical active space automatically. |
| `vrc_tst_sample_basis` | `"avdz-f12"` | Master key. aDZ-F12 orbital basis for orientation sampling. |
| `vrc_tst_active_space` | nested object/list: detailed below / nested object/list: detailed below | Branch hook. Selected electron/orbital counts; per-bond surface overrides bind the actual checkpoint. See nested controls below. |
| `rotdpy_run` | `1` | Branch hook. Submit and complete RotdPy sampling rather than only writing input. |
| `rotdpy_validation_profile` | `"production"` | Branch hook. Require production correction/surface/convergence records before accepting flux. |
| `rotdpy_mess_mode` | `"production"` | Branch hook. Use validated Rotd state counts in MESS; do not silently fall back to PST for a failed VRC calculation. |
| `run_messpf` | `1` | Branch hook. Execute MESSPF and parse standalone species partition functions for thermochemistry. |
| `me` | `1` | Master key. Build the MESS network. PES workers use me=2 for intermediate files; the PES controller assembles the final network. |
| `run_me` | `1` | Master key. Execute MESS after accepted ANL/thermochemistry/VRC handoff. |
| `collider` | `"He"` | Master key. Helium bath gas. |
| `epsilon` | `300.0` / `350.0` | Master key. Lennard-Jones well depth of the reacting system, in epsilon_unit. These example transport estimates are provisional, not fitted kinetic evidence. |
| `epsilon_unit` | `"K"` | Master key. Input energy unit; resolved worker JSON records cm-1 after conversion to avoid converting twice. |
| `sigma` | `4.5` / `5.0` | Master key. Lennard-Jones collision diameter in angstrom. Provisional transport estimate; pressure-dependent rates depend on this choice. |
| `queuing` | `"slurm"` | Master key. Slurm submission backend. |
| `queue_name` | `"month-long-cpu"` | Master key. Blodgett partition; the runtime discovers its actual time limit. |
| `ppn` | `8` | Master key. Eight requested processor slots for ordinary QC jobs. |
| `queue_job_limit` | `10` | Master key. At most ten simultaneous allocations within a KinBot worker calculation tree. |
| `electronic_state_models` | nested object/list: detailed below / nested object/list: detailed below | Branch hook. General species-level electronic data: S, CH3S and SH measured levels for MESS/PF; S's separate explicit ANL centroid correction; spatial terms for state planning. Not an automatic VRC SOC potential. |
| `homolytic_bonds` | nested object/list: detailed below / nested object/list: detailed below | Master key. Only the initial-well C-S bond, or C-S and S-S bonds, are requested for hom_sci. Atom indices here are ZERO based. The initial alias was added on this branch. |
| `vrc_tst_scan_shift` | `0.1` | Branch hook. Real CASPT2 level shift in hartree for relaxed scans: 0.1. |
| `vrc_tst_scan_ipea` | `0.25` | Branch hook. CASPT2 IPEA shift in hartree for relaxed scans: 0.25. |
| `vrc_tst_spin_splitting` | nested object/list: detailed below | Branch hook. Ground high-spin CC anchor plus target-minus-reference CASPT2 gaps, extrapolated aDZ/aTZ with exponent 3.7. See nested controls. |
| `vrc_tst_scan_method` | `"caspt2-f12(5,4)" / "caspt2-f12(7,5)"` | Master key. Explicit selected method for constrained relaxed geometries. |
| `vrc_tst_scan_basis` | `"avtz-f12"` | Master key. aTZ-F12 basis for constrained relaxed scans. |
| `vrc_tst_scan_gradient` | `"numerical"` | Branch hook. Numerical derivatives of the declared target CASPT2-F12 root. The composite spin-gap-corrected energy is not differentiated. |
| `vrc_tst_sample_method` | `"caspt2-f12(5,4)" / "caspt2-f12(7,5)"` | Master key. Orientation-sampling method with the explicit selected CAS. |
| `vrc_tst_high_method` | `"caspt2-f12(5,4)" / "caspt2-f12(7,5)"` | Master key. Trusted MR method; the spin_splitting model uses its explicit gap levels and CC anchor. |
| `vrc_tst_high_basis` | `["avdz-f12","avtz-f12"]` | Master key. Two-basis MR estimate, aDZ-F12 and aTZ-F12. |
| `vrc_tst_high_extrapolation_power` | `3.7` | Branch hook. 3.7 in X^(-p) extrapolation. A user-requested two-basis estimate, not a demonstrated error bar for F12 spin gaps. |
| `vrc_tst_correction_model` | `"spin_splitting"` | Branch hook. Combine the spin-gap-corrected relaxed cut with the frozen-fragment sampling cut, including the configured-level geometry-relaxation contribution. |
| `single_point_ppn` | `8` | Master key. Eight processor slots for ordinary single-point work. |
| `vrc_tst_max_nodes` | `10` | Branch hook. VRC allocation ceiling; admission is shared with concurrent automatic CBH/ANL work. |
| `vrc_tst_min_stack_mw` | `1000` | Branch hook. Minimum Molpro stack per rank, in million 8-byte words; the dispatcher checks total memory including overhead. |
| `rotdpy_max_jobs` | `10` | Branch hook. Maximum concurrent Rotd jobs for this PES; ten, not an additional ten concurrent with its unfinished earlier branches. |
| `rotdPy_mem` | `1000` | Master key. Molpro sampling stack in million words per rank. |
| `rotdpy_processors` | `8` | Branch hook. Eight MPI processors per sampling allocation. |
| `slurm_exclude_nodes` | `["node1","node9","node10","node13","node14","node16","node21","node23","node24","node30","node33"]` | Branch hook. The eleven user-identified unusable nodes. The list is explicit site data, never hardcoded in general chemistry code. |
| `rotdpy_exclusive` | `1` | Branch hook. Request exclusive Slurm nodes for Rotd jobs. |
| `vrc_tst_memory_mb` | `104000` | Branch hook. 104000 MB total requested memory per VRC allocation, covering all rank stacks and overhead. |
| `simultaneous_kinbot` | `1` | Master key. One PES worker at a time, so ten jobs per worker remains ten per PES rather than multiplying across wells. |
| `rotdpy_flux_parameters` | nested object/list: detailed below | Branch hook. Monte Carlo minimum/maximum counts and 1 percent relative statistical target. A sample ceiling does not establish convergence. |
| `vrc_tst_surface_model` | `nested inline object` | Inline channels, physical surfaces, per-bond CAS counts and basis-specific orbital indices/checkpoints, states, weights and solver controls. No second configuration file. |
| `sella_fmax` | `5e-4` | Final L1/L2 maximum force norm 0.0005 eV/angstrom, followed by native Hessian validation. Gaussian Opt keywords do not control Sella. |
| `anl_geometry_fmax` | `5e-4` | Final automatic L3 geometry maximum force norm 0.0005 eV/angstrom; configurable for every molecule. |
| `anl_scf_convergence` | `1e-8` | SCF 1e-8 control across DFT L1 and L2: Gaussian density Conver=8; Molpro/CFOUR/MRCC backend-mapped threshold. |
| `anl_cc_convergence` | `1e-7` | Coupled-cluster convergence 1e-7, mapped to each ANL backend. |
| `vrc_tst_scf_convergence` | `1e-8` | Molpro VRC SCF convergence 1e-8. |
| `vrc_tst_scan_fmax` | `5e-4` | Maximum constrained force norm for Molpro relaxed scans, eV/angstrom. |
| `vrc_tst_correlation_convergence` | `1e-7` | Molpro VRC CC energy convergence 1e-7; local CASPT2 controls are specified in each basis model. |
| `minimum_recovery_attempts` | `2` | Branch hook. Bounded Gaussian L1/L2 minimum descent and reoptimization, with fresh native Hessian each time. |
| `minimum_recovery_displacement` | `0.02` | Branch hook. Full Cartesian displacement in angstrom; second attempt uses 0.04. Both signs must be evaluated. |
| `imagfreq_threshold` | `50.0` | Existing minimum/TS frequency policy, in cm^-1. No threshold relaxation was used to accept the failed minima. |
| `anl_minimum_recovery` | `attempts=2, displacement_angstrom=0.02, difference_step_angstrom=0.01, imaginary_threshold_cm_inverse=1.0` | Branch hook. Check unconstrained ANL minima before dependent energies; finite-difference force Hessian with rigid motions projected out, bounded downhill displacement/reoptimization. Adds 6N force evaluations per curvature check. Later native Molpro harmonic acceptance remains mandatory. `attempts=0` disables this extra check/recovery, not native harmonic validation. |

## Nested controls

**Active space.** `mode=selected` uses the explicit `active_electrons` and `active_orbitals`; `source` is provenance. The separate channel overrides select (5,4), (7,5), (7,5). A reviewed Molpro record supplies the ordered subspace; the integers alone do not identify orbitals.

**Spin-splitting.** `reference.method=uccsd(t)-f12b` and `reference.basis=avtz-f12` specify the restricted-reference high-spin CC anchor. Each `gap_levels` entry supplies a stable `label`, a CASPT2-F12 method with explicit CAS counts, the two `basis` values, `extrapolation_power=3.7`, `shift=0.1`, and `ipea=0.25`. `weighting.model=arithmetic_mean` combines alternative method estimates; with one CASPT2 estimate it is the identity. It does not average doublet and quartet energies or pretend that zero method spread means zero uncertainty.

**Monte Carlo.** `pot_smp_min=100`, `pot_smp_max=1000`, `tot_smp_min=100`, `tot_smp_max=15000`, `flux_rel_err=1` (percent), `smp_len=1` are forwarded to RotdPy. Inspect whether convergence was reached before the cap. Numerical convergence and total kinetic uncertainty are different quantities.

**Electronic species data.** `smiles`, optional `smiles_aliases`, `charge`, and `ground_multiplicity` match the species. `levels` contains measured `energy_cm_inverse`, physical `degeneracy`, spin `multiplicity`, and a descriptive `label`. The S atom has (0,5), (396.055,3), (573.640,1). The declared `spin_free_reference=term_centroid` gives -195.756111 cm^-1 to the atomic ground-J energy in ANL, separate from thermal populations in MESS/PF. `term_orbital_degeneracy=3` for S and 2 for CH3S/SH supplies separated-fragment spatial-root counting; it is not multiplied into resolved-level degeneracies. CH3S has (0,2), (255,2), and SH has (0,2), (376.835,2). Both molecules also explicitly set `spin_free_reference=term_centroid`: SH adds -188.4175 cm^-1 and CH3S adds -127.5 cm^-1 once in ANL, retaining scalar anharmonic ZPE. CH3S is the approved effective JT/SO approximation, not an exact coupled spin-vibronic treatment. `spin_orbit_correction_hartree` is the general explicit alternative to a level-derived correction; the two options are mutually exclusive. Current input corrections are used for cached references during assembly. `source` records the experimental citations. S2 retains its ground triplet declaration. The [spectroscopy review](../../../../docs/SULFUR_ELECTRONIC_LEVELS_2026-10-17.md) documents the RRHO approximation, JT limitations and inventory. No distance-dependent SOC potential is introduced.

**Inline `vrc_tst_surface_model` (branch hook to existing RotdPy/MESS).** `schema=1` is the format, `reviewed` and `review_note` record the explicit physical model. `channels[].bond` is ONE based; `[1,2]` is C-S and `[2,3]` is S-S. `combination=single` and a sole `bottlenecks[].name=single` select ordinary VRC. `coordinate=reactive` uses reactive-atom pivots; `distances_angstrom` lists 3.0 through 20.0 in 0.2 steps. `parameters` contains only the per-bond active-space and complete inline `vrc_tst_state_model` overrides. Each `states` entry gives doublet `multiplicity=2`, selected `root`, `weight=1`, and a common `energy_offset_cm_inverse=0` relative to ground separated fragments, explained by `energy_reference_note`. CH3S includes the two spatial components correlating with its ground 2E state; CH3SS uses its ground doublet per bond. This is a spin-free adiabatic model assumption. Extra SA roots and quartet reference roots are not extra reactive flux. Spin degeneracy is applied once by the assembler. There is no inner/outer stage, rapid equilibration assumption or automatically inserted intermediate well.

**Inline `vrc_tst_state_model` (branch hook to native Molpro features).** `schema=2`, `state_average_policy=reviewed`, `active_space`, and `molecule.symbols/charge` bind counts and atom order. `basis_models` selects a separate reviewed record for each basis. `orbital_selection.mode=molpro_record`, `path`, `sha256`, `record=2140.2`, and `orbital_type=natural` identify the exact seed bytes. `active_indices` explicitly lists the one-based natural-orbital indices in that record: `[11,12,13,14]`, `[18,19,20,21,22]`, or `[18,19,20,21,22]` in distinct channel records; `align_active_orbitals=false` preserves the native-confirmed protocol without adding DIAB. No seed permutation is repeated. `canonical_orbitals=true` writes 2141.2 for CASPT2 while retaining the natural record. `state_average` contains `multiplicity`, `roots`, and normalized within-spin `weights`: 3+3, 1+1, and 6+6 for the three channels; equal spin-block weights are implicit. `selected_roots` gives a fixed root unless overridden by an explicit surface or energy policy.

`energy_selection={"4":"lowest"}` evaluates the full declared correlated quartet manifold with the same shared reference and chooses the minimum for the ground-state CC correction anchor in all three channels. CH3S needs this at the reviewed center point; one-root CH3SS C–S adds no work, while S–S retains protection if the reference order changes during the relaxed cut. It is a reusable energy-only policy; it does not claim state-character tracking or permit a mismatched gradient. Reconnaissance always reports individual root energies, even when a production model contains this selector.

`casscf_options`: `solver=so_sci`, `exact_diagonalization=true`, `quasi_newton=true`, `energy_convergence=1e-8`, `gradient_convergence=1e-5`, `max_iterations=999`. `correlation_options`: `reference_policy=state_average` sends the full same-spin REFSTATE and equal WEIGHT values; `energy_convergence=1e-7`, `coefficient_convergence=1e-5` gives native THRVAR=1e-10, `max_iterations=2000`, `max_reference_iterations=2000`, `max_reference_vectors=60`. These are not the obsolete THRVAR=1e-20 settings. Molpro numerical failure invalidates the whole point.

## Relevant implicit settings

**Native orbital diagnostics (general branch hook).** Each inline basis state model may override `orbital_diagnostics`: `orbital_coefficient_threshold=0.05`, `ci_coefficient_threshold=0.05`, `virtual_orbitals=2`. Defaults enable documented native printing for CAS/MRCI without changing convergence or electronic settings. Each native output gets an adjacent `.orbital_audit.json`, extracting actual occupations, AO coefficients and states with source hashes. See [the print/report contract](../../../../docs/MOLPRO_NATIVE_ORBITAL_DIAGNOSTICS.md).

The new reusable DFT preset is explicit in both example inputs; the old UMA
presets and implicit ANL preset selection remain available. No new electronic
structure method is implemented. Gaussian supplies B3LYP/B2PLYP and D3(BJ);
KinBot selects and consistently forwards those settings. `SCF=(XQC,Conver=10)` and
`Integral=UltraFine` are retained. An unrestricted open-shell calculation is
selected by multiplicity, without forcing unrestricted singlet calculations.

For a Gaussian level profile, omitted `scan_method/scan_basis` and
`barrierless_saddle_method/basis` now inherit the appropriate level instead
of legacy MP2/6-31G or B3LYP/6-31G defaults. Explicit user overrides survive.
The `_well_mp2` filename remains a legacy label for the scan-reference job;
in these inputs its actual method is the L1 B3LYP-D3(BJ)/cc-pVDZ profile.
Preliminary constrained TS/ring paths honor explicit profiles, and IRC/path
writers retain their profile's dispersion/options. Legacy non-profiled AM1
searches are unchanged. The optional semiempirical conformer search is off.

`minimum_recovery_attempts=2` and `minimum_recovery_displacement=0.02` Å
are new general controls for Gaussian/Sella minima with native Hessians.
A rejected imaginary mode is tried in both signs with separate checkpoint
copies, followed by reoptimization and a new Hessian; attempt two uses 0.04 Å.
`imagfreq_threshold=50` cm⁻¹ remains the existing acceptance cutoff. Recovery
cannot accept a rejected spectrum, change method/spin, or descend a TS.
`minimum_recovery_attempts=0` disables it. See the
[bounded retry policy](../../../../docs/MINIMUM_RECOVERY_2026-10-16.md).


The constrained correction scan defaults to 2.5–19.9 angstrom in 0.2 steps plus a 30 angstrom asymptote. Numerical F12 gradients make this more expensive than the fixed-geometry confirmation. The minimum over the scanned dividing surfaces must be inspected for boundary minima; a finite grid cannot guarantee the global bottleneck.

RotdPy uses the existing 0.85 local recrossing factor once. Its temperature, energy and angular grids default to [10,10,1.05,70], [0,10,1.05,190], [0,1,1.1,80] in RotdPy grid notation. Slurm retry budgets and all native inputs/results retain execution history. This is not a new trajectory-based dynamical calculation.

MESS explicitly uses 300–5000 K in 100 K increments and ModelEnergyLimit 800 kcal/mol. Unchanged defaults: pressures 7.6,76,760,7600,76000 Torr; direct solver; EnergyStepOverTemperature 0.2; ChemicalEigenvalueMax 0.2. The explicit ExcessEnergyOverTemperature 30 retains its usual value. Exponential-down energy relaxation uses Factor 200 cm^-1, Power 0.85, ExponentCutoff 15. These transport/collision assumptions are distinct from the 400 kcal/mol L2 reaction-search cutoff and need scrutiny for quantitative pressure-dependent rates. The general temperature, energy-grid and barrier defaults are unchanged; these are per-example JSON choices.

All reaction families are searched on each explored well, but the explicit homolytic-bond map enables only the three reviewed initial-well VRC channels. Newly found isomers do not automatically inherit a sulfur seed with matching element counts. Geometry/TS/IRC checks determine the accepted network. The candidate inventory is not a prediction that every candidate is a real pathway.

The full resolved configuration, including defaults, is saved by the preparation audit. This document describes the explicitly supplied parameters and the consequential defaults; it does not claim that every numerical backend has run end to end on these molecules.

## Single-input selection and numerical criteria

All editable production configuration is in `ch3s.json` or `ch3ss.json`.
`vrc_tst_surface_model.channels[].bottlenecks[].parameters.vrc_tst_state_model.basis_models[BASIS].orbital_selection.active_indices`
is the exact location of the index list. The associated WFU is binary data,
not another configuration; its path, SHA-256, record and orbital type are in
the same JSON. Atom order and charge are checked. A raw MO number without
that source context does not specify an orbital. Existing external-file inputs
remain readable for compatibility; the production handoff and recon installer
write one JSON and do not require them.

For CAS(5,4) on a 25-electron molecule the general formula is
`closed=(25-5)/2=10`, `occ=10+4=14`. For CAS(7,5) with 41 electrons it
gives `closed=17`, `occ=22`. These are calculated from the input, not
hardcoded sulfur branches. Noncontiguous user-selected indices are packed
with documented Molpro `ROTATE,i.1,j.1,0` exchanges. Already packed records
require no exchange; continued checkpoints explicitly use the packed indices
and preserve the alignment choice. See the [general interface](../../../../docs/ORBITAL_REVIEW_WORKFLOW.md).

Final L1/L2 and automatic L3 use `fmax=5e-4` eV/angstrom (about
9.72e-6 Eh/bohr), a tight force criterion. Sella does not use Gaussian Berny's
complete force/displacement test. Native Hessian validation remains required;
preliminary constrained reaction-search steps keep their exploration criteria.
Molpro VRC relaxation uses the explicit `vrc_tst_scan_fmax=5e-4` force threshold.
SCF Conver=8 controls Gaussian density convergence, not a 1e-8 energy-error bound.
Molpro CC uses local energy 1e-7 and squared-amplitude 1e-10; CASPT2 uses
THRDEN=1e-7 and THRVAR=1e-10. CASSCF orbital convergence is unchanged.
See the [restart numerical contract](../../../../docs/PRODUCTION_CONVERGENCE_RESTART_2026-10-17.md).
