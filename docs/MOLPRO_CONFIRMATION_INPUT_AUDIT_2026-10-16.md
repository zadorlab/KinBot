# Molpro confirmation input audit and restart

The first confirmation exposed three input defects. Local Python tests did not
validate native Molpro execution, and the original readiness claim was too
strong. This record separates native evidence, documented syntax, repairs and
remaining native validation. The user handles SSH. No active-space conclusion
is changed by these input failures.

## Native evidence supplied by the user

Run: `/home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016`.
The paths below remain in place after the scoped restart. These are excerpts
from the user's native-output reports, not newly computed results.

| Native output relative to that run | Reported failure | Interpretation and repair |
| --- | --- | --- |
| `ch3s_cs_adz_caspt2/adz__far/adz__far.out` (52908060); corresponding `ch3ss_cs_adz_caspt2` (52908063) and `ch3ss_ss_adz_caspt2` (52908066) | `TRYING TO OVERWRITE DUMP RECORD 2140.2 WITH DIFFERENT BASIS DIMENSION` | Projection read the seed, then the destination collided with the source. Stage the foreign-basis copy on FILE,3 and write a fresh FILE,2. |
| `ch3s_cs_atz_caspt2/atz__center/atz__center.out` (52908069); corresponding `ch3ss_cs_atz_caspt2` (52908071) | `RHF mode requires a single determinant wavefunction!` in `rhf_initialize` | SO-SCI reported MS2=3 with one singly occupied orbital, despite local WF requesting a doublet. The exact internal override is not established. Reset global spin and explicitly supply the matching local WF/OCC/CLOSED determinant instead of relying on restart defaults. |
| `ch3s_cs_atz_mrci/atz__far/atz__far.out` (52908062) | `not sure if it works with nstatr/=nstati` in `cif12solve.F90` | A lone root-2 reference gave incompatible internal/reference counts. Use all declared same-spin references and explicit matching NSTATI. Reject unsupported lone excited-root MRCI-F12 requests during generation. |

The monitor now retains nonempty Molpro question-mark diagnostics, including
these actual causes, before MPI teardown messages. It still distinguishes a
live log warning from a terminal scheduler/execution failure.

## Official manual review

Checked against the online Molpro manual on October 16, 2026 UTC. The scope is
every feature emitted by this 39-point confirmation and its native executor,
not every feature elsewhere in KinBot. Some manual sections lag the installed
version; those discrepancies are explicitly identified below.

| Feature used | Audit outcome and source |
| --- | --- |
| `MEMORY`; permanent `FILE,2` / `FILE,3`; `NEW`; short basenames | Memory precedes files. The projected source is a separate staged copy; FILE,2,NEW cannot erase it. Original reviewed seeds are hash-checked and untouched. [File handling](https://www.molpro.net/manual/doku.php?id=file_handling). |
| Cartesian XYZ, `ANGSTROM`, `ORIENT,NOORIENT`, `SYMMETRY,NOSYM` | Fixed coordinates/orientation; electronic C1. C3v/Cs describes the reference geometries, not non-Abelian electronic irreps. [Geometry](https://www.molpro.net/manual/doku.php?id=molecular_geometry), [symmetry](https://www.molpro.net/manual/doku.php?id=general_program_structure#symmetry). |
| `RHF,SO-SCI`, local `ENERGY/MAXIT`, `WF/OCC/CLOSED` | Require `OCC+CLOSED=NELEC` and `OCC-CLOSED=2S`. Explicit occupations replace inherited defaults. [SCF program](https://www.molpro.net/manual/doku.php?id=the_scf_program). |
| `SET,SPIN` and local WF precedence | Reset restart spin before each RHF; subsequent local MULTI/CI WF cards specify their own spin. [Wavefunction definition](https://www.molpro.net/manual/doku.php?id=general_program_structure#defining_the_wavefunction_wf). |
| `START,record,TYPE=NATURAL,OVL` | Projection is a guess. Foreign-basis source 2140.3 and destination 2140.2 are distinct. [Orbital directives](https://www.molpro.net/manual/doku.php?id=general_program_structure#orbital_directives), [different-basis restriction](https://www.molpro.net/manual/doku.php?id=the_scf_program#starting_with_previous_orbitals). |
| MULTI `OCC/CLOSED`, WF/STATE/WEIGHT blocks | Explicit target CAS; weights follow their WF/STATE block. No CONFIG/SELECT restriction excludes mixed-spin determinant averaging. [MULTI](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi). |
| MULTI `SO-SCI`, `EXACT_DIAG=1`, `QN=1`, `ACCURACY` | Documented optimizer controls. E=1e-8 and gradient=1e-5. [MULTI options](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#options_of_the_mcscf_program). |
| MULTI `MAXIT,999` | Manual's legacy MAXITER paragraph says 40; this installed executable explicitly reported a 999 cap and has accepted 999. Do not claim portability of that cap to all releases. |
| `ROTATE,...,0`, `DIAB,2140.2` | Zero-angle exchange follows START once at the initial partition. DIAB permits equal source/save records at subsequent same-basis geometries. It does not prove state identity. [MULTI orbital handling](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#restoring_and_saving_the_orbitals_and_ci_vectors). |
| `CANONICAL,2141.2`; `NATORB,2140.2` | Separate pseudo-canonical/natural sets; explicit TYPE when selected. Default density is state averaged. [MULTI orbital output](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#saving_the_final_orbitals). |
| `RS2-F12,SHIFT=0.3,IPEA=0.25`; `STATE,1,k` | Single-state corrections on SA orbitals; not XMS. Pseudo-canonical orbitals are required for the non-invariant IPEA implementation. ENERGY(1) selects the sole computed state, including the level-shift energy correction. [RS perturbation theory](https://www.molpro.net/manual/doku.php?id=multireference_rayleigh_schroedinger_perturbation_theory). |
| `MRCI-F12,SCALF12=1,CABS_SINGLES=1` | Retain documented defaults. Prior CABS orthogonalization failure remains a native numerical risk, not repaired by this input audit. [MRCI-F12](https://www.molpro.net/manual/doku.php?id=the_mrci_program#explicitly_correlated_mrci-f12). |
| MRCI `STATE,1,k`; `REFSTATE,N,1,...,N`; `OPTION,NSTATI=N` | Full same-spin contraction-reference manifold, one correlated root per call. N is 3/1/6 by channel. Native restriction overrides the manual's broad compatibility statement. [MRCI state and reference controls](https://www.molpro.net/manual/doku.php?id=the_mrci_program#defining_reference_state_numbers). |
| `OPTION,CLUSTER=0`; `ENERGD(1)` | Unchanged relaxed-Davidson selector; not fixed/rotated corrections or the variational energy. Capture after the entire command, including CABS. [Cluster energies](https://www.molpro.net/manual/doku.php?id=the_mrci_program#cluster_corrections_for_multi-state_mrci). |
| `THRDEN=1e-8`; `THRVAR=1e-10`; `MAXIT/MAXITI=2000`; `MAXVI=60` | THRVAR is a squared coefficient-change criterion. Iteration and vector limits are distinct documented controls. Native `***` reflects integer print width, not a missing setting. [MRCI thresholds and options](https://www.molpro.net/manual/doku.php?id=the_mrci_program). |
| `basis=avdz-f12/avtz-f12`; automatic auxiliary bases | Both aliases resolved in the installed library before the reported failures. Preserve native basis-library/auxiliary selections; alias acceptance alone does not validate a wavefunction. [Basis input](https://www.molpro.net/manual/doku.php?id=basis_input), [F12 auxiliary selection](https://www.molpro.net/manual/doku.php?id=explicitly_correlated_methods). |
| `GTHRESH`, `GPRINT`, assignments and `SHOW` | Explicit local convergence overrides; printed markers retain each calculated root. Check native counts/convergence and complete marker sequence, not process exit alone. [General structure](https://www.molpro.net/manual/doku.php?id=general_program_structure). |
| `PUT,XML` and its `ORBITAL` subcommand | Keep high-angular-momentum orbital exports and explicitly select the natural set. [PUT](https://www.molpro.net/manual/doku.php?id=molecular_geometry#writing_files_for_postprocessing_put), [XML](https://www.molpro.net/manual/doku.php?id=post-processing_of_output_and_databases). |
| Launch `-d`, `-W`, `-t 1`, `-n 8`, `-o` | WFU repository is the job directory; scratch is separate. 1000 MW/rank gives nominal 64 GB across eight ranks; 104000 MB allocation leaves overhead below the reported 126000 MB capacity. GA/actual peak memory still needs native observation. [Launcher](https://www.molpro.net/manual/doku.php?id=running_molpro), [parallel memory](https://www.molpro.net/manual/doku.php?id=running_molpro_on_parallel_computers). |

The MRCI reference-policy change is explicit catalogue data (`state_average`),
not an automatic numerical retry. It changes the contraction space, so old
selected-reference MRCI values must not silently substitute for the new series.
It does not change CASSCF weights, CASPT2 reference policy, IPEA/shift, CABS
treatment or the Davidson selector. No new geometry optimization, XMS, SOC,
perturbative-CABS fallback or altered Hamiltonian has been added.

## Restart and validation limits

`restart_on_blodgett.sh` prepares a fresh `_r2` run, checks seed/input hashes,
then uses directory-scoped cleanup to stop the old confirmation controller
and allocations. All old native files, validated points and failed evidence
remain at their original paths. The replacement records the old manifests and
counts in `previous_confirmation.json`. All 39 points are recomputed under
the corrected, consistent protocol; old results are retained as comparison
evidence. An unavailable scheduler prevents cancellation/submission. Repeating
the script cannot start a second controller for an already-started replacement.

Regression coverage executes all 39 input/staging transformations with a fake
launcher, verifying source/destination separation, untouched source hashes,
RHF determinant counts, matching MRCI reference/internal counts, unchanged
energy selectors, and final-checkpoint discovery. Separate tests cover scoped
stop/preservation, scheduler/seed failures and the actual reported diagnostic
sentences. These are software checks, **not native Molpro runs**. No local
Molpro executable is available; the corrected native results remain pending.

The final repair regression passed **1,629 tests and 391 subtests** using
`.venv/bin/python -m pytest -q`. Two unrelated warnings concerned deprecated
Python `imp` and sandbox access to a Sella/JAX compilation cache. Shell syntax
and `git diff --check` also passed. The original handoff's smaller test count
is historical, not evidence that its native input defects were covered.

## Restart catalogue comparison repair

The subsequent Blodgett restart prepared all 39 points, then rejected its own
catalogue before cancellation or submission. `build()` writes normalized JSON,
but the restart compared that snapshot's hash to the differently formatted
source file. The original regression used the already-normalized copy and
missed the actual CLI path. The corrected regression first reproduced this
failure with three source formats. The repair compares parsed catalogue
contents; `validated_record()` still verifies the saved snapshot's exact hash.
Changed scientific settings and modified snapshots remain rejected before
stopping old work. Repeating the restart reuses the prepared `_r2` directory.
The targeted confirmation and cleanup suite passed **27 tests** after this
repair; no native Molpro input, energy convention or scientific setting changed.
