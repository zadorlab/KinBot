# Completed compact-resolution panel: 25 October 2026

The supplied monitor at 15:42:50 PDT reports **60/60 completed calculations**,
15 candidates with four jobs each, and no recorded startup/manual/allocation
retries. The campaign remains `complete_review_required`: no selection was
published by completion. The supplied text is preserved in
[resolution_character_evidence.txt](resolution_character_evidence.txt), with
[provenance and checksum](resolution_evidence_provenance.json).

The subsequent native capture `compact_resolution_native.20261025T224531Z.txt`
was supplied inline in the conversation. It contains the four near calculations
discussed below; the numerical excerpts here are a transcription, not a claim
that its full native text or binary orbitals have been exported locally.

## Native initialization and two unresolved cases

All four native executions report `native_xml_scf_seed_same_basis_v8`, correct
SCF checkpoint hashes, the requested active indices, and normal convergence.
They use `START,2100.3,TYPE=CANONICAL` without forced OVL. No recurrence of the
earlier VIRTORB error or sulfur-core substitution is shown. These results
therefore require scientific subspace review, not an automatic numerical retry.

**CH2 + HSS, parent 791671450790040000002, C1–S2 homolysis.**
CAS(3,3) uses source SCF orbitals 20,21,22 from 30 Å, one quartet and one doublet,
with weights 0.5 each. The far active orbitals contain the two carbon frontier
directions and HSS radical, each with occupation 1. At 2 Å, active 20 is largely
S2 p, 21 is carbon p, and 22 is S2/S3-rich with little printed carbon character.
This does not demonstrate retention of the intended scission span.

| Near result | avdz-f12 | avtz-f12 |
| --- | ---: | ---: |
| Active occupations, 20/21/22 | 1.49685 / 0.99991 / 0.50323 | 1.49675 / 0.99991 / 0.50334 |
| Q1 energy / Eh | -834.491777860787 | -834.505223988596 |
| D1 energy / Eh | -834.647636437801 | -834.662966176376 |
| Q1−D1 / kcal mol^-1 | 97.803 | 98.985 |
| D1 coefficient of `2a0` | 0.99830413 | 0.99825479 |
| D1 coefficient of `0a2` | -0.05388571 | -0.05481726 |
| Native elapsed seconds | 22.11 | 191.82 |

The quartet is exactly `aaa` in this three-orbital space. Averaging its 1/1/1
occupation with the nearly 2/1/0 doublet explains the near 1.5/1/0.5 occupations.
Those fractional state-averaged occupations do not by themselves demonstrate
strong correlation within D1 or validate orbital character. A large near spin
gap likewise does not alone establish an error. Character and state evolution
between endpoints remain unresolved.

**C(SH2) + SH, parent 791871710370000000002,
R_Addition_CSm_R_5_1_2 product-complex exit.** The CAS(7,5) trial correctly starts
from far SCF 18,19,20,21,23. Molpro prints the requested 22↔23 exchange as a
90-degree rotation. It uses 17 closed orbitals, seven active electrons in five
orbitals, and two doublets with weights 0.5/0.5.

At 30 Å, natural active 20/21 carry the SH p pair with occupations 1.5/1.5;
active 18/19 are occupied molecular valence orbitals and 22 is carbon-rich
acceptor character. A carbon p frontier is inactive 17. At the L2 complex,
inactive 17 instead has SH p character, active 20 is carbon-p-rich, and 21/22
contain mixed carbon/SH character. Thus the added orbital did not remove the
near/far exchange. Agreement between DZ and TZ establishes reproducibility of
this solution, not its suitability for the intended outer potential.

| Near result | avdz-f12 | avtz-f12 |
| --- | ---: | ---: |
| Active occupations, 18–22 | 1.99452, 1.96810, 1.48756, 1.30512, 0.24470 | 1.99445, 1.96752, 1.48732, 1.30632, 0.24439 |
| D1 energy / Eh | -834.476172327048 | -834.490671654038 |
| D2 energy / Eh | -834.440662433313 | -834.456282538867 |
| D2−D1 / kcal mol^-1 | 22.283 | 21.579 |
| D1 coefficient of `22a20` | 0.87317861 | 0.87387584 |
| D2 coefficient of `222a0` | 0.93077318 | 0.93034406 |
| Native elapsed seconds | 15.89 | 155.43 |

These are fixed-geometry CASSCF energies, not correlated production barriers
or rate estimates. CI strings refer to the displayed natural orbital order.
For CH2/HSS `RECONENERGY=energy(1)` is Q1 because the quartet block comes first;
for C(SH2)/SH the final `MULTI/... energy` line is D2. Neither scalar is a generic
reactive D1 selector. Use the explicit multiplicity/root energies.

## Inventory of all twelve pending requests

Atom order is C1 S2 S3 H4 H5 H6. SCF indices below always refer to the particular
request's own basis-specific record 2100.2; optimized natural indices refer to
2140.2. Similar product names do not authorize sharing ordinal indices between
requests. The first two rows retain the previously reviewed 56-job evidence.
All other rows use the new 60-job report, with the earlier compact comparisons
where stated. None is newly published by this audit.

| Parent and reaction suffix | Compact evidence and indices | Assessment |
| --- | --- | --- |
| 791871710370000000002 / `hom_sci_1_2` | CH + H2SS, CAS(5,4), SCF 17,18,21,22 → natural 19–22; two equal doublets. Sulfur donor plus CH frontier span. | Endpoint character supported in both bases; no new trial needed. |
| Same / `hom_sci_2_3` | CHSH2 + S, CAS(5,4), SCF 18,20,21,22 → natural 19–22; three D and three Q roots. All three atomic-S p directions plus fragment radical. | Endpoint character supported; equal weights within each block and equal total D/Q block weights. |
| Same / `R_Addition_CSm_R_5_1_2` | C(SH2) + SH outer exit, CAS(7,5), SCF 18,19,20,21,23 → natural 18–22; two equal D roots. | Held: near/far character exchange persists; continue this same space before another size change. |
| Same / `r12_insertion_R_2_1_5` | HCS + H2S outer exit. CAS(3,3), SCF 19,21,22 → natural 20–22 retains HCS pi/radical/pi-star. HCS-only CAS(5,4), SCF 18,19,21,22 → natural 19–22 adds occupied in-plane HCS valence. One D root. | Both consistently retain HCS character in both bases. CAS(3,3) remains the smaller supported candidate. The larger space changes near-minus-far CASSCF energy by about -0.22/-0.21 kcal/mol; quantitative CASPT2 convergence is not established. |
| 791671450790040000002 / `hom_sci_1_2` | CH2 + HSS, CAS(3,3), SCF/natural 20–22; one D and one Q, equal total weights. | Held: near character is unresolved, as detailed above. |
| Same / `r12_insertion_R_2_1_5` | CH2 + HSS outer exit, CAS(3,3), SCF/natural 20–22; one D and one Q. | Carbon-frontier/HSS-radical span retained at both endpoints in both bases. |
| Same / `r12_insertion_R_2_1_6` | CH2 + HSS outer exit, same counts and indices from its own checkpoint. | Same endpoint support, with its distinct atom map retained. |
| Same / `hom_sci_2_3` | CH2S + SH, CAS(5,4), pi start 18,19,21,22 versus donor start 19,20,21,22 → natural 19–22; two equal D roots. | Starts converge to the same solutions. Far pi/pi-star plus SH pair and near valence bond mixing support the pi-start compact candidate. Detailed state continuity remains a production scan check. |
| Same / `r12_insertion_R_2_3_4` | CH2S + SH outer exit, same two hypotheses from this exit's own seeds. | Both converge to the same pi/pi-star plus SH-pair solutions. Near/far character agrees in both bases. |
| Same / `r12_insertion_R_5_1_6` | HCSSH + H5 outer exit; CAS(3,2), SCF 19,21 versus 20,21 → natural 20/21; one D. | Both give the same CASSCF energy and 2/1 occupations, while retaining different occupied donors. The data do not establish that either donor is necessary. |
| Same / `r12_insertion_R_6_1_2` | HCSSH + H6 outer exit; same alternatives from this endpoint's own seeds. | Same donor ambiguity; distinct L2 internal geometry retained. |
| Same / `r12_insertion_R_6_1_5` | HCSSH + H6 outer exit; same alternatives from this endpoint's own seeds. | Same donor ambiguity; distinct L2 internal geometry retained. |

The HCSSH/H data support a single-reference description at the tested endpoints.
An H-only CAS(1,1) is the formal radical minimum, but it was not run in this panel;
the donor comparisons do not validate its CASPT2 potential. Neither a unique
donor selection nor its necessity is claimed here. No extra HCSSH calculations
are submitted by the two-case continuation plan.

Optimization roots are orbital/state support, not extra reactive flux. These
tables do not assign additional quartet flux or multiply a ground-parent rate
by the number of averaged roots. The outer exits remain product-side separation
problems, with their inner TSs represented separately.

## Next bounded diagnostic and publication status

The [continuation plan](../r7_compact_continuation_20261025/README.md) has
**14 CASSCF jobs for only the two unresolved subspaces**, starting from each
basis's completed 30 Å natural checkpoint. It preserves CAS sizes, state
weights and internal fragment geometry. It follows intermediate distances to
the same near targets to compare with the independent SCF-start results.
It is a finite diagnostic, not a relaxed production scan or automatic retry.

The current data do not justify saying all twelve selections are finalized.
Publication remains a single later operation after resolving the two character
questions and documenting the final HCSSH/H compact model choice. The already
published five models and the earlier CH3S selections are not changed here.
No fresh production monitor was included in this capture, so their current
scan/reconciliation progress and live cap23 are not independently verified by
this report. The continuation inherits the existing shared cap and exclusions;
it does not reset them or cancel allocations.
