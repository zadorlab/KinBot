# Production input parameters and provenance

Baseline: `origin/master` as recorded at commit `5c6525a7` (5 October 2026). “Branch” includes earlier work on `composite`, not only the latest fix. Existing backend algorithms are distinguished from KinBot orchestration. Values are copied from the actual inputs; the complete JSON remains authoritative.

| Parameter | CH3S / CH3SS value | Hook origin and purpose |
| --- | --- | --- |
| `title` | `"ch3s_anl1_f12"` / `"ch3ss_anl1_f12"` | Master key. Run label; no change to chemistry. |
| `smiles` | `"[CH3][S]"` / `"C[S][S]"` | Master key. Initial molecular connectivity; explicit hydrogens are generated in the recorded atom order. |
| `charge` | `0` | Master key. Neutral total charge. |
| `mult` | `2` | Master key. Doublet parent; fragment multiplicities are assigned separately. |
| `barrier_threshold` | `1000.0` | Master key. Search cutoff in kcal/mol. 1000 keeps high candidate barriers available; it does not certify their stationary points or enlarge the MESS energy grid. |
| `conformer_search` | `1` | Master key. Search conformations before final thermochemistry. |
| `high_level` | `1` | Master key. Refine accepted structures at the L2 profile. |
| `rotor_scan` | `1` | Master key. Build hindered-rotor potentials where chemically applicable. |
| `theory_preset` | `"uma-b2plyp-anl"` | Branch hook. L1 UMA-s-1p2/OMol exploration; L2 Gaussian B2PLYP-D3(BJ)/cc-pVTZ with Sella optimization and native harmonic Hessians. |
| `composite_method` | `"ANL1-F12"` | Branch hook. Automatic CBH/ATcT plus ANL1-F12 component workflow for accepted species and stationary points; not a CASPT2 kinetic sampling method. |
| `anl_f12_augmented` | `1` | Branch hook. Use fully augmented QZ-F12/5Z-F12 reference bases for ANL1-F12 thermochemistry. This is distinct from the VRC aDZ/aTZ pair. |
| `vrc_tst_accuracy_preset` | `"anl1-f12"` | Branch hook. CASPT2-only VRC defaults; explicit values below override the preset. It does not select a chemical active space automatically. |
| `vrc_tst_sample_basis` | `"avdz-f12"` | Master key. aDZ-F12 orbital basis for orientation sampling. |
| `vrc_tst_active_space` | nested object/list: detailed below / nested object/list: detailed below | Branch hook. Selected electron/orbital counts; per-bond surface overrides bind the actual checkpoint. See nested controls below. |
| `vrc_tst_state_model_file` | `"../caspt2_production/ch3s_cs_state_model.json"` / `""` | Branch hook. Reviewed basis-specific seed/state model; a blank CH3SS global entry is filled by its two separate bond overrides. |
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
| `electronic_state_models` | nested object/list: detailed below / nested object/list: detailed below | Branch hook. General species-level electronic data: sulfur fine structure for MESS/PF and its explicit ANL centroid correction; CH3S spatial degeneracy for state planning. Not an automatic VRC SOC potential. |
| `homolytic_bonds` | nested object/list: detailed below / nested object/list: detailed below | Master key. Only the initial-well C-S bond, or C-S and S-S bonds, are requested for hom_sci. Atom indices here are ZERO based. The initial alias was added on this branch. |
| `vrc_tst_surface_model_file` | `"ch3s_surface_model.json"` / `"ch3ss_surface_model.json"` | Branch hook. Ordinary single-bottleneck channels, reactive doublet components and per-bond orbital records. It does not enable inner/outer routing. |
| `vrc_tst_scan_shift` | `0.1` | Branch hook. Real CASPT2 level shift in hartree for relaxed scans: 0.1. |
| `vrc_tst_scan_ipea` | `0.25` | Branch hook. CASPT2 IPEA shift in hartree for relaxed scans: 0.25. |
| `vrc_tst_spin_splitting` | nested object/list: detailed below | Branch hook. Ground high-spin CC anchor plus target-minus-reference CASPT2 gaps, extrapolated aDZ/aTZ with exponent 3.7. See nested controls. |
| `vrc_tst_scan_method` | `"caspt2-f12(2,2)"` | Master key. Method template for constrained relaxed geometries. CAS(2,2) here is replaced by the selected per-bond counts before rendering; it is not the executed sulfur active space. |
| `vrc_tst_scan_basis` | `"avtz-f12"` | Master key. aTZ-F12 basis for constrained relaxed scans. |
| `vrc_tst_scan_gradient` | `"numerical"` | Branch hook. Numerical derivatives of the declared target CASPT2-F12 root. The composite spin-gap-corrected energy is not differentiated. |
| `vrc_tst_sample_method` | `"caspt2-f12(2,2)"` | Master key. Orientation-sampling method template, replaced by CAS(5,4) or CAS(7,5). |
| `vrc_tst_high_method` | `"caspt2-f12(2,2)"` | Master key. Trusted MR method template; the spin_splitting model uses its explicit gap levels and CC anchor. |
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
| `rotdpy_flux_parameters` | nested object/list: detailed below | Branch hook. Monte Carlo minimum/maximum counts and 5 percent relative statistical target. A sample ceiling does not establish convergence. |

## Nested controls

**Active space.** `mode=selected` uses the explicit `active_electrons` and `active_orbitals`; `source` is provenance. The separate channel overrides select (5,4), (7,5), (7,5). A reviewed Molpro record supplies the ordered subspace; the integers alone do not identify orbitals.

**Spin-splitting.** `reference.method=uccsd(t)-f12b` and `reference.basis=avtz-f12` specify the restricted-reference high-spin CC anchor. Each `gap_levels` entry supplies a stable `label`, a CASPT2-F12 method template, the two `basis` values, `extrapolation_power=3.7`, `shift=0.1`, and `ipea=0.25`. `weighting.model=arithmetic_mean` combines alternative method estimates; with one CASPT2 estimate it is the identity. It does not average doublet and quartet energies or pretend that zero method spread means zero uncertainty.

**Monte Carlo.** `pot_smp_min=100`, `pot_smp_max=1000`, `tot_smp_min=100`, `tot_smp_max=15000`, `flux_rel_err=5` (percent), `smp_len=1` are forwarded to RotdPy. Inspect whether convergence was reached before the cap. Numerical convergence and total kinetic uncertainty are different quantities.

**Electronic species data.** `smiles`, optional `aliases`, `charge`, and `ground_multiplicity` match the species. `levels` contains measured `energy_cm_inverse`, physical `degeneracy`, spin `multiplicity`, and a descriptive `label`. The S atom has (0,5), (396.055,3), (573.640,1). The declared `spin_free_reference=term_centroid` gives -195.756111 cm^-1 to the atomic ground-J energy in ANL, separate from thermal populations in MESS/PF. `term_orbital_degeneracy=3` for S and 2 for CH3S supplies separated-fragment spatial-root counting; it is not multiplied into already resolved J degeneracies. `source` records the data citation. S2 has a ground triplet declaration. No molecular spin-vibronic level spectrum or distance-dependent SOC correction is invented.

**Surface files (branch hook to existing RotdPy/MESS).** `schema=1` is the format, `reviewed` and `review_note` record the explicit physical model. `channels[].bond` is ONE based; `[1,2]` is C-S and `[2,3]` is S-S. `combination=single` and a sole `bottlenecks[].name=single` select ordinary VRC. `coordinate=reactive` uses reactive-atom pivots; `distances_angstrom` lists 3.0 through 20.0 in 0.2 steps. `parameters` contains only the per-bond active-space and state-file overrides. Each `states` entry gives doublet `multiplicity=2`, selected `root`, `weight=1`, and a common `energy_offset_cm_inverse=0` relative to ground separated fragments, explained by `energy_reference_note`. CH3S includes the two spatial components correlating with its ground 2E state; CH3SS uses its ground doublet per bond. This is a spin-free adiabatic model assumption. Extra SA roots and quartet reference roots are not extra reactive flux. Spin degeneracy is applied once by the assembler. There is no inner/outer stage, rapid equilibration assumption or automatically inserted intermediate well.

**State files (branch hook to native Molpro features).** `schema=2`, `state_average_policy=reviewed`, `active_space`, and `molecule.symbols/charge` bind counts and atom order. `basis_models` selects a separate reviewed record for each basis. `orbital_selection.mode=molpro_record`, `path`, `sha256`, `record=2140.2`, and `orbital_type=natural` identify the exact seed bytes; `align_active_orbitals=false` preserves the native-confirmed protocol without adding DIAB. No seed permutation is repeated. `canonical_orbitals=true` writes 2141.2 for CASPT2 while retaining the natural record. `state_average` contains `multiplicity`, `roots`, and normalized within-spin `weights`: 3+3, 1+1, and 6+6 for the three channels; equal spin-block weights are implicit. `selected_roots` gives a fixed root unless overridden by an explicit surface or energy policy.

`energy_selection={"4":"lowest"}` evaluates the full declared correlated quartet manifold with the same shared reference and chooses the minimum for the ground-state CC correction anchor in all three channels. CH3S needs this at the reviewed center point; one-root CH3SS C–S adds no work, while S–S retains protection if the reference order changes during the relaxed cut. It is a reusable energy-only policy; it does not claim state-character tracking or permit a mismatched gradient. Reconnaissance always reports individual root energies, even when a production model contains this selector.

`casscf_options`: `solver=so_sci`, `exact_diagonalization=true`, `quasi_newton=true`, `energy_convergence=1e-8`, `gradient_convergence=1e-5`, `max_iterations=999`. `correlation_options`: `reference_policy=state_average` sends the full same-spin REFSTATE and equal WEIGHT values; `energy_convergence=1e-8`, `coefficient_convergence=1e-5` gives native THRVAR=1e-10, `max_iterations=2000`, `max_reference_iterations=2000`, `max_reference_vectors=60`. These are not the obsolete THRVAR=1e-20 settings. Molpro numerical failure invalidates the whole point.

## Relevant implicit settings

The constrained correction scan defaults to 2.5–19.9 angstrom in 0.2 steps plus a 30 angstrom asymptote. Numerical F12 gradients make this more expensive than the fixed-geometry confirmation. The minimum over the scanned dividing surfaces must be inspected for boundary minima; a finite grid cannot guarantee the global bottleneck.

RotdPy uses the existing 0.85 local recrossing factor once. Its temperature, energy and angular grids default to [10,10,1.05,70], [0,10,1.05,190], [0,1,1.1,80] in RotdPy grid notation. Slurm retry budgets and all native inputs/results retain execution history. This is not a new trajectory-based dynamical calculation.

MESS defaults: 300–2000 K in 100 K increments; pressures 7.6,76,760,7600,76000 Torr; direct solver; ModelEnergyLimit 400 kcal/mol; EnergyStepOverTemperature 0.2; ExcessEnergyOverTemperature 30; ChemicalEigenvalueMax 0.2. Exponential-down energy relaxation uses Factor 200 cm^-1, Power 0.85, ExponentCutoff 15. These transport/collision assumptions are distinct from the 1000 kcal/mol reaction-search cutoff and need scrutiny for quantitative pressure-dependent rates.

All reaction families are searched on each explored well, but the explicit homolytic-bond map enables only the three reviewed initial-well VRC channels. Newly found isomers do not automatically inherit a sulfur seed with matching element counts. Geometry/TS/IRC checks determine the accepted network. The candidate inventory is not a prediction that every candidate is a real pathway.

The full resolved configuration, including defaults, is saved by the preparation audit. This document describes the explicitly supplied parameters and the consequential defaults; it does not claim that every numerical backend has run end to end on these molecules.
