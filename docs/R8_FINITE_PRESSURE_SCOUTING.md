# R8 scouting with no explicit well cutoff

Use `WellCutoff=-1` as the scouting reference. The chemical model is unchanged: corrected B2PLYP L2 energies, the retained native geometry/frequency data, RRHO/PST surfaces, the documented symmetry/mirror factors, and the native triplet-S2 restoration. Positive cutoff 10 is a sensitivity control, not a prerequisite for using these results. Other native grid and collision-kernel limits still apply.

The 193 nm calculation initializes the parent well at 148.142152 kcal/mol (the native grid point is 148.1441874), with a 300 K, 760 torr helium bath. This is a statistical ground-state-PES assumption that puts that energy into the molecule. It is not a conversion of photon energy to a bath temperature or a prediction of an experimental quantum yield.

## Full-model 193 nm branching

| Parent | Products | Fraction (%) |
|---|---|---:|
| CH3S | CH2S + H | 99.586047 |
| CH3S | HSC + H2 | 0.209833 |
| CH3S | CH3 + S | 0.161861 |
| CH3SS | CH3 + triplet S2 | 95.610633 |
| CH3SS | CH2S + SH | 2.087050 |
| CH3SS | CH2SS + H | 1.058212 |
| CH3SS | HC(=S)SH + H | 0.531891 |

CH3S `pr_3` is **HSC**, with H bonded to S in the retained native geometry (`[H][S]#[C]`); it is not silently reassigned to HCS. All remaining product fractions and exact native identifiers are in `scouting_summary.json`.

## Cost and pruning sensitivity

Counts are MESS kinetic wells (including vdW complexes), product sets, and represented connections—not counts of unique QC calculations or approved VRC models. The selection combines the prior thermal HPL and fixed-energy scenarios, rather than selecting from the photochemical products alone. In particular, CH3S w2 remains important thermally despite little 193 nm flux.

| System / scope | Wells | Product sets | Connections | 193 nm distribution shift (percentage points) |
|---|---:|---:|---:|---:|
| CH3S full | 3 | 6 | 10 | 0 |
| CH3S 1pct candidate | 2 | 1 | 3 | 0.413953 |
| CH3S 0.2pct candidate | 2 | 3 | 5 | 0.042260 |
| CH3SS full | 17 | 18 | 91 | 0 |
| CH3SS 1pct candidate | 7 | 3 | 10 | 1.323978 |
| CH3SS 0.2pct candidate | 11 | 8 | 28 | 0.185680 |

The shift is total variation: half the sum of the absolute changes in all product fractions. It bounds the change for any grouping of products in this specific calculation. It is not a bound on electronic-structure/PST errors or a guarantee that every small branch has a small relative error. Removed sinks contribute zero only in this explicitly pruned-model distribution comparison; native `***` rate entries are never filled with zero.

The lean CH3SS candidate changes CH2S + SH from **2.087050% to 2.857446%**, a 36.9% relative shift, while CH3 + S2 changes from 95.610633% to 96.164215%. It is useful for identifying the dominant chemistry, but is less suitable for resolving the few-percent competing branches. The broader candidate preserves more of those competing exits. Neither scope has been installed into R9 or used to approve a CAS model.

The broader CH3SS candidate gives CH2S + SH = **2.210517%** (versus 2.087050%), CH2SS + H = **0.972797%** (versus 1.058212%), and HC(=S)SH + H = **0.575245%** (versus 0.531891%). Its maximum relative change among these three branches is about 8.2%. For planning accurate competition among the few-percent branches, prefer this broader scope; use the lean scope only if the stated losses are acceptable. It removes 63 of 91 CH3SS connections (69.2%) and six of 17 kinetic wells, while the CH3S broader scope removes five of ten connections and one of three wells.

The earlier collisionless sensitivity panel gave worst combined-removal distribution shifts of 0.944/3.380 percentage points (CH3S/CH3SS, lean candidate) and 0.0876/0.5349 percentage points (broader candidate). Those scenarios include approximate 0.5/2 channel-flux variations; they are not rigorous chemical uncertainty bounds.

## What was verified

- Native MESS hot-energy branching agrees with independent propagation of the exported native, unreduced transition rates to maximum absolute product differences of `4.71e-7` (CH3S) and `7.05e-7` (CH3SS). The native output is retained, including tiny signed entries. This conclusion does not rely on an alternative propagator alone.
- Native MESS also completed both pruned candidates for each PES. Their 193 nm product fractions agree with independent propagation of each corresponding native generator to a maximum absolute difference below `6.7e-7`. The reduction recommendation therefore has native MESS support for the full and reduced models.
- At 193 nm, halving the energy-grain fraction from 0.5 to 0.25 changes the full-model product distribution by `1.13e-8` / `4.00e-8` in total variation. Both use the same native initial energy bin. The separate 140 kcal/mol comparison does not have that identical-bin property.
- Time-step refinement gives changes below `3.1e-16` in the terminal 193 nm product distribution; probability-conservation errors are below `4.7e-14`. No final renormalization or clipping was used in the independent propagation. Terminal-yield convergence does not assert identical transient populations at every time.
- All six full/candidate fine-grid generator exports have finite nonnegative transition and sink rates, complete state coverage, independently consistent collision diagonals and detailed-balance discrepancies below `1.2e-13`. This is an additional structural check, not chemical validation.
- The supplied initial distribution dissociates essentially completely on the tested `1e-7 s` interval. The bath collision model is helium; air or nitrogen has not been substituted.
- For the broader scope, the summed incoming-traversal estimates for removed wells are below 0.00390% (CH3S) and 0.0887% (CH3SS) on the tested `1e-7 s` interval. These are conservative first-visit upper-bound estimates from the native-generator dynamics; they are separate from the combined-removal product-distribution comparison.
- The optional integrated channel counts describe expected traversals. Repeated circulation can make them exceed one. They are not silently presented as first-visit probabilities.

## Thermal finite-pressure status

Both QD well-reduction controls completed all requested temperatures (300–2000 K) and pressures (0.01, 0.1, 1, 10, 100, 760 torr). Conservative-QD direct CH3S controls also completed all ten temperatures at 760 torr. These are numerical controls, not approval of every table.

The direct CH3S correction removes the gross 300 K loss/outgoing-rate inconsistency: the slow w1→w2 rate is about `4.1e-13 s^-1`, rather than a numerical loss floor near `1e-6 s^-1`. At 400 K the direct loss is `1.981250959e-7 s^-1`; the QD well-reduction outgoing sum agrees, but its `8.66e-7 s^-1` diagonal still does not. Tiny negative couplings remain in some tables and are kept as numerical evidence.

At high temperature, the selected chemical subspace disappears for some conditions (CH3S direct: 1800/2000 K at 760 torr; CH3SS well-reduction: 1500–2000 K on this grid). Separate original-well phenomenological rates then have no representation in that reduction. Their `***` entries remain undefined; this is not a reason to invent rates or delete the corresponding PES chemistry.

A separate ordinary-double direct-method control gave grossly nonconserving 193 nm branching (`2.45e15` and `7.02e18` in a product fraction). This is not a tiny sign error. That control is rejected; the reported 193 nm fractions use the native QD well-reduction result verified against the unreduced generator. The direct method remains a preferred benchmark, but the manual also recognizes its low-temperature numerical limitations. No method is selected merely because it prints positive coefficients.

The larger CH3SS high-precision direct tests have not supplied a converged replacement for the full thermal grid; bounded DD controls reached their time limits. Association rows and the lowest-temperature loss summaries have not all been validated. Thus the deliverable supports the stated 193 nm/1-atm scouting conclusion and the separately validated HPL/microcanonical scope, **not certification of the entire finite-pressure rate matrix**. Resolving that wider thermal-rate task requires completed precision/convergence checks and justified kinetic population definitions; cutoff choice or positivity alone cannot substitute for those checks.

## Evidence and native-source accounting

- Raw corrected chemical inputs remain `local_runs/mess_review_20261102/deliverables/{ch3s,ch3ss}/mess_0000.inp`; the old delivery manifest is unchanged and has been reverified.
- New evidence is under `local_runs/mess_review_20261102/finite_pressure`: `scouting_summary.json`, `native_hot_comparison.json`, `native_candidate_comparison.json`, `scouting_generator_audit.json`, `unreduced_hot_convergence.json`, `193nm_bin_alignment.json`, `thermal_initial_well_comparison.json`, `experiments.json`, `recommended_workscope.json`, and their native calculation directories.
- [Exact MESS source changes, mathematical explanations, per-file line counts and limitations](../validation/native_mess/README.md). Seven patch layers reconstruct the saved source hashes from the original archive; compiler/ABI changes are distinct from numerical-method experiments.
- [Short symmetry/stereochemistry additions and attribution](SYMMETRY_STEREO_ADDITIONS.md). The original QC export and the current chemical models are preserved.

No extra Blodgett chemistry or controller refresh is required for this postprocessing result. Updating KinBot for the published audit helpers uses the usual `git pull --ff-only origin composite`; the experimental native MESS patches are not silently installed by that update.
