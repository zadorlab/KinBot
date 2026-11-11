# R8 provisional pruning, 2026-11-10

The corrected local scouting models support a lean R9 planning scope. The full model retains the unresolved entrance complex provisionally. The recommended lean models omit that complex.

## Corrections and verification

- Coupled Gaussian energy and gradient requests and bounded Cartesian minimum cycling recovery are published in `fba04b7e6`, with the complete Sella/Gaussian entry-point regression in `0c836befe`. The captured cycling pattern is detected at step 82 of the 453-frame trajectory. A native Gaussian optimization with the new procedure remains unverified.
- The same commit separates the colliding legacy fragment ID `771210770330000000002` into the observed H-C-S-S and S-C-S-H identities. Their H2 product thresholds are 82.28 and 43.71 kcal/mol. Every energy candidate is checked against its complete identity before cross-worker selection.
- Both captured PESs were rebuilt through the normal postprocessing entry point with QC submission disabled. Their original calculation artifacts remain the source evidence. Composition and endpoint checks passed for all 19 CH3S and 133 CH3SS models in the final inputs.
- The full-run electronic-state declarations were applied locally, including the CH3S, SH, S and CH level ladders and centroid shifts. The earlier native triplet-S2 result replaces the singlet-S2 model. The CH3 + S2 threshold is 43.59 kcal/mol. Declared corrections are applied once.
- Two reverse searches share one saddle with a proper-map RMSD of 0.000524 angstrom, mass-weighted imaginary-mode overlap 0.99999065, and energy difference 0.0000056 kcal/mol. The saddle contributes once. The H2 + H-C-S-S reverse TS remains a distinct channel after the fragment-ID correction.
- A native fixed-width table formatting case was also found. Full-width negative numbers can touch the preceding value or `***`. The rate auditor now reads those fields from their native header boundaries and preserves their values and undefined entries.

The independent native energy-grained checks have nonnegative finite elementary rates, collision-diagonal and detailed-balance errors below 1e-13, and propagation mass errors below 1.3e-14. Halving the propagation time step changes terminal product probabilities by at most 1.2e-15.

## Full-model results at 193 nm

The parent receives 148.142152 kcal/mol in the statistical ground-state model, with a 300 K, 760 torr helium bath. Native MESS samples 148.144 kcal/mol on its grid. These fractions describe the supplied energized-parent model.

| PES | Product set | Native MESS fraction |
|---|---|---:|
| CH3S | CH2S + H | 99.5860% |
| CH3S | CH3 + S | 0.1619% |
| CH3S | HSC + H2 | 0.2099% |
| CH3SS | CH3 + triplet S2 | 95.5624% |
| CH3SS | HCS + H2S | 1.6236% |
| CH3SS | CH2SS + H | 1.0597% |
| CH3SS | CH2S + SH | 0.9481% |
| CH3SS | HC(=S)SH + H | 0.3375% |

The verified reverse abstraction connects the CH2S + SH encounter well to the HCS + H2S encounter well. In the independent full-model propagation, adding the reverse routes changes HCS + H2S from 0.0951% to 1.6279% and CH2S + SH from 2.1396% to 0.9495%. This route is retained in the lean scope.

Native hot branching and unreduced propagation agree within 8.9e-8 in product probability for CH3S and 1.46e-4 for CH3SS. The latter is 0.0146 percentage point. The native QD well-reduction result and unreduced dynamics therefore support the same pruning scale, while retaining a measurable reduction discrepancy.

## Recommended lean scope

| PES | Full wells / products / connections | Lean wells / products / connections | Native 193 nm distribution change |
|---|---|---|---:|
| CH3S | 3 / 6 / 10 | 2 / 2 / 4 | 0.2522 percentage point |
| CH3SS | 18 / 20 / 95 | 5 / 4 / 9 | 0.8884 percentage point |

The change is total variation between the full and explicitly pruned product distributions. At 140–155 kcal/mol it spans 0.172–0.318 percentage point for CH3S and 0.590–1.232 percentage points for CH3SS. These comparisons quantify the chosen model's pruning error. Electronic-structure, PST and experimental-energy-distribution uncertainties remain separate.

CH3S retains CH3S, CH2SH, CH2S + H and the user-requested CH3 + S channel. Its four connections are `nobar_1`, `nobar_4`, `rxn_1` and `rxn_4` in this pass's full input. The small HSC + H2 channel is a one-connection extension through `rxn_3`.

CH3SS retains the following five populations.

| Local label | Population | Native identity or source |
|---|---|---|
| w_1 | CH3SS | 791471490990000000002 |
| w_3 | CH2SSH | 791671450790040000002 |
| w_6 | HSCH2S | 791472110370000000002 |
| w_9 | CH2S + SH encounter complex | 791671450790040000002_intra_H_migration_suprafacial_1_4_IRC_F_prod |
| w_12 | HCS + H2S encounter complex | 791871560500020000002-sb6e34ded5baaf03b_12_shift_S_R_3_2_4_IRC_F_prod |

Its nine retained connections are listed below. These labels are specific to the new full input.

| Connection | Endpoints | Role |
|---|---|---|
| nobar_1 | w_1 → CH3 + S2 | Parent C-S dissociation |
| nobar_2 | w_1 → CH2SS + H | Parent H loss |
| rxn_2 | w_1 ↔ w_3 | Isomerization |
| rxn_9 | w_3 ↔ w_6 | Isomerization |
| rxn_37 | w_3 ↔ w_9 | Inner route into CH2S + SH complex |
| nobar_14 | w_6 → CH2S + SH | Existing represented exit |
| nobar_39 | w_9 → CH2S + SH | Outer capture/dissociation |
| rxn_reverse_93359f49 | w_9 ↔ w_12 | Native verified abstraction saddle |
| nobar_41 | w_12 → HCS + H2S | Outer capture/dissociation |

The lean input preserves the distinction between the two encounter complexes. The extra high-energy reversible side wells were removed as an explicit planning choice and tested in native MESS. The lean CH3SS native fractions are 96.0773% CH3 + S2, 1.8429% HCS + H2S, 1.1023% CH2S + SH and 0.9775% CH2SS + H. Minor-channel relative errors are larger than the overall distribution error.

A single-channel extension adds the scout's `nobar_16` exit from HSCH2S to HC(=S)SH + H. It adds one product set and one connection, with no additional well. Native 193 nm total variation falls to 0.6949 percentage point. That branch becomes 0.1935%, compared with 0.3375% in the full model, so this extension retains the channel for further refinement while leaving substantial branching uncertainty. The retained scout representation remains the starting approximation for a later tight-TS/VRC treatment.

## Provisional complex and thermal implications

The unresolved entrance uses its accepted L1 geometry and twelve positive L1 frequencies, with its 5.4924 kcal/mol L1 binding depth anchored to the L2 fragments and their declared spin–orbit shift. Its resulting energy is approximately 83.3622 kcal/mol above CH3SS. This is the explicit provisional substitution.

The upper bound from summed incoming traversals for visiting this complex at 193 nm is 0.00582%. The expected forward traversals of its new abstraction saddle are 0.00178 per hundred initial molecules. Changing the binding depth by plus or minus 2 kcal/mol changes any terminal product probability by less than 3.4e-8 across 140–155 kcal/mol. This is an illustrative sensitivity test. The lean scope excludes this complex and its dependent route.

Native HPL calculations cover all ten temperatures from 300 to 2000 K. Conditional on dissociation, CH2S + H accounts for at least 99.9480% of CH3S products. CH3 + S2 accounts for at least 99.8636% of CH3SS products. CH2SH remains useful in the CH3S thermal pathway even with a small hot-energy branch. Thermal HPL supports a smaller network than the energized-parent competition alone requires.

Native finite-pressure tables also cover all 60 requested temperature/pressure conditions for each PES. The first CH3SS calculation reached its 1200-second bound, and a completed local continuation supplied the remaining high-temperature conditions. Original outputs and execution receipts are retained. The tables contain 43 CH3S and 305 CH3SS negative entries, together with undefined reductions and known low-temperature loss inconsistencies. Their individual phenomenological coefficients remain uncertified. Pruning here uses the verified hot-energy observables and the separately checked HPL limits.

## Saved handoff

The local deliverables are under `local_runs/r8_reverse_integration_20261110/deliverables`.

- `ch3s/full.inp`, `ch3ss/full.inp` contain the corrected full models.
- `ch3s/lean.inp`, `ch3ss/lean.inp` contain the recommended lean models.
- `ch3ss/h_loss_extension.inp` contains the one-channel extension.
- `scouting_summary.json`, `thermal_audit.json`, `chemical_audit.json`, `model_provenance.json` files and the native run subdirectories retain the quantitative checks.
- `manifest.json` hashes the delivered artifacts. The original exports and unsuccessful diagnostic attempts remain in their original evidence directories.

The current checks cover identity, composition, endpoint references, native saddle duplication, declared electronic corrections, numerical propagation and explicit pruning sensitivity. The Gaussian convergence of the provisional complex and certification of every reduced finite-pressure coefficient remain outside the verified results. No additional MESS source changes were needed in this pass. The existing experimental native build and its earlier source accounting are documented in `validation/native_mess/README.md`.
