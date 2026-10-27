# Closure tests for the remaining known orbital requests

The [34-channel ledger](../r7_detached_review_20261027/selection_ledger.md)
preserves all installed and previously prepared choices. This panel does not
rerun them. Its ten screens address eleven unresolved requests, since the
two new HCS/H2S reaction labels have the same verified product endpoint.

All source files come from the completed detached export, bound by native
WFU/output hashes in each recipe. [Source evidence](source_evidence.json)
records the source manifests and archive digest. The [preceding native audit](../r7_detached_review_20261027/README.md)
explains which hypotheses failed and why. No recipe is a CAS approval.

| Test | Jobs, both bases combined | Question to close |
| --- | ---: | --- |
| Two HCSH+SH parents | 8 | Do their own CAS(5,4) checkpoints retain HCSH pi/pi-star and the complete SH pi pair, as in the already reviewed analogous exit? |
| New stereoisomer C–S homolysis | 8 | Does the own-source S2 donor restore the bonding frontier in CAS(5,4), compared with CAS(3,3)? |
| New stereoisomer S–S homolysis | 8 | Is atomic S 3p plus the molecular radical, CAS(5,4), sufficient; does occupied bound-S participation require CAS(7,5)? |
| New CH2/HSS outer exit | 8 | Does the earlier outer-exit CAS(3,3) logic hold for this exact reverse product endpoint, with its own donor control? |
| New HCSSH+H5 outer exit | 4 | Confirm the own H5 SOMO CAS(1,1), guided by the completed H4 and earlier H exits. No redundant one-hole donor controls. |
| CH/HSSH failed donor control | 12 | Fixed CASCI separately opens source occupied 18 and 19; optimized CAS(7,5) contains both alternatives and the same frontier. |
| HCS/H2S failed donor control | 12 | Same controlled tests distinguish the H2S lone pair from substituted occupied HCS character. |
| HCSH/SH virtual reorganization | 4 | Fixed CASCI adds source occupied 18 without allowing the existing active frontier to rotate out. |
| C(SH2)/SH branch correlation | 26 | Follow the same 4D/2Q CAS(5,4) manifold outward and back from independent near/far seeds; identify how the singlet-fragment D3/D4 asymptote connects. |
| Total | **90** | Three jobs maximum in the existing detached pool. |

Fixed CASCI controls hold the orbital partition and geometry fixed; natural
orbitals can still rotate within that partition. Explicit closed-orbital
swaps select the intended source direction before opening it. Molpro documents
both the [CAS-CI option and zero-angle orbital interchange](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi).
Those controls are not publishable production seeds. Optimized union spaces
are tests of the competing orbital choices, not a rule to select the largest
CAS or lowest total energy. Near/far inactive ordinals alone do not identify
the same chemical direction; comparisons must use the exported orbital character.

The branch panel reuses the already demonstrated SH retention. Its new
question is state correlation, including initialization dependence, not
whether the old two-doublet setup can be made to work. It does not assume
that the lowest doublet at every distance has the requested singlet-fragment
asymptote. If the native results contradict that assignment, report the
contradiction rather than silently relabeling the products or roots.

## Evaluation and handoff

Require native normal completion, the explicit state manifold, immutable
inputs/WFU/XML, occupations and CI character at every anchor. Check complete
fragment degeneracies and the association bonding directions. Compare active
spans only using an appropriate AO metric; same-geometry reconstruction is not
a cross-geometry wavefunction overlap. Compare relative potentials with each
model's own asymptote and distinguish CI-space changes from orbital relaxation.
Any final selection must name its own two-basis seeds, orbital roles, reactive
roots, supporting roots and retained limitations.

The staging tests simulate all 90 jobs, predecessor adoption, the shared cap,
repeated invocation and unchanged production/source files. They establish
orchestration behavior only. Native completion and scientific conclusions
require the exported results from Blodgett.

Run `bash examples/anl/orbital_review_batches/r7_closure_20261027/run_blodgett.sh`.
When complete, run the sibling `collect_blodgett.sh`. Neither script resumes
production or writes `_r7` calculations. They do not create `_r8`.
