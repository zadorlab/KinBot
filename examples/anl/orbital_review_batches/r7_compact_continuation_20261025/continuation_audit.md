# Completed inward-continuation review, 25 October 2026

The supplied monitor at **16:20:12 PDT** reports **14/14 completed jobs**:
eight CH2/HSS calculations and six C(SH2)/SH calculations. Startup, manual
startup and allocation retry counts are all zero. The controller reports
`complete_review_required`; no production model was accepted. These were
CASSCF-only calculations, so the monitor's zero correlated-root count is
expected. This report does not establish completion of production VRC scans,
CASPT2 corrections, RotdPy or either PES.

The exact supplied text is saved in
[continuation_character_evidence.txt](continuation_character_evidence.txt).
[continuation_evidence_provenance.json](continuation_evidence_provenance.json)
contains its checksum and a transcription of all 14 energy/occupation records.
The full native archive has now been transferred and verified locally: all
4,056 members match the manifest. The [completed native audit](../r7_compact_union_20261025/README.md)
records the archive checksum, complete coefficient comparisons, latest driver
status, the coordinate-discovery fix and the bounded next probes. The earlier
text below describes the evidence available before that transfer.

## CH2 + HSS: CAS(3,3), one quartet and one doublet

Parent `791671450790040000002`, reaction `hom_sci_1_2`, C1–S2 separation.
The model and 0.5/0.5 state weights are unchanged. Natural orbitals 20–22
from 30 Å were propagated through 6, 3.5, 2.5 and 2 Å in each basis.

| Separation / Å | DZ D1 / Eh | TZ D1 / Eh | DZ occupations 20/21/22 |
| --- | ---: | ---: | --- |
| 6 | -834.58998041 | -834.60451300 | 1.00229 / 1.00000 / 0.99771 |
| 3.5 | -834.59138179 | -834.60589833 | 1.09898 / 0.99997 / 0.90105 |
| 2.5 | -834.61856782 | -834.63343039 | 1.36908 / 0.99978 / 0.63114 |
| 2 | -834.64763644 | -834.66296617 | 1.49685 / 0.99991 / 0.50323 |

The final energies agree with the independent-start near calculations to
about 1e-8 Eh, and the printed character agrees. A direct 30 Å-to-near SCF
initialization jump therefore does not explain away the earlier result.

At 6, 3.5 and 2.5 Å, active 20 and 22 have appreciable carbon-frontier/HSS
radical mixing; active 21 retains the other carbon p direction. At 2 Å,
20 is mainly S2 p, 21 remains carbon p, and 22 is S2/S3-rich with little carbon
among the printed leading coefficients. Inactive 17 has carbon/sulfur character.
The change in the leading active character is bracketed by 2.5 and 2 Å.
These six-coefficient summaries do not quantify subspace overlap or identify
whether the change is continuous orbital relaxation or exchange with inactive
orbitals. They do not justify declaring a discontinuity, an incorrect physical
state, or successful preservation of the intended scission span.

The previous full near CI excerpt showed a nearly `2a0` doublet and an `aaa`
quartet. Consequently, state-averaged occupations near 1.5/1/0.5 must not be
interpreted as a strongly correlated doublet on their own. The inward result
reproduces the numerical solution; the state/subspace interpretation remains
to be resolved from the complete native evidence.

## C(SH2) + SH outer exit: CAS(7,5), two equal doublets

Parent `791871710370000000002`, reaction `R_Addition_CSm_R_5_1_2`.
This is the product-complex outer separation, with its own product-side L2
geometry. The final point is that actual complex, not an imposed 2 Å geometry.
The active indices are 18–22 in the natural checkpoints. Counts and two
doublet weights of 0.5 each are unchanged.

| Point | DZ D1 / Eh | TZ D1 / Eh | DZ occupations 18–22 |
| --- | ---: | ---: | --- |
| 6 Å | -834.47772946 | -834.49228309 | 1.99958, 1.93741, 1.50000, 1.50000, 0.06300 |
| 3.5 Å | -834.47587419 | -834.49044709 | 1.99956, 1.93772, 1.50000, 1.49999, 0.06273 |
| L2 complex | -834.44481459 | -834.46011449 | 1.99810, 1.93716, 1.50129, 1.49831, 0.06514 |

The inward near solution retains both SH p orbitals as active 20/21. Active
19 is carbon-s-rich, 22 is a carbon-p-rich acceptor, and 18 is occupied
S2/carbon valence character. A carbon-p-rich orbital remains inactive 17.
The TZ near occupations are 1.99814, 1.93659, 1.50124, 1.49835, 0.06567.

The independent-start near solution instead put an SH p orbital in inactive
17 and carbon p character in active 20. It also had substantially lower
energies. The two routes therefore give different converged solutions at the
same requested geometry, active-space size and state average.

| Inward minus independent near energy / kcal mol^-1 | DZ | TZ |
| --- | ---: | ---: |
| Root D1 within each solution | +19.677 | +19.175 |
| Root D2 within each solution | +22.179 | +22.803 |
| Equal-weight state average | +20.928 | +20.989 |

These differences use 627.5094740631 kcal mol^-1 per Eh. Root labels denote
energy ordering within each calculation, not demonstrated correspondence of
the electronic wavefunctions between solutions. Relative to the completed
30 Å D1 values, the inward complex D1 lies +20.398/+19.944 kcal mol^-1
(DZ/TZ), versus +0.720/+0.769 for the independent-start solution. These are
fixed-geometry CASSCF differences, not correlated binding energies or barriers.

This is reproducible initialization dependence in both bases. Retaining the
SH pair does not establish that the higher solution supplies the desired
adiabatic ground potential. Choosing the lower energy independently at every
point would also require a state/subspace continuity check. Neither branch is
newly accepted here.

The previous source SCF set 18,19,20,21,23 already included the carbon frontier
and SH p hypotheses. The optimized natural order is different. Repeating those
SCF indices is not a new hypothesis, and natural indices must not be treated
as chemically equivalent to the original SCF indices. No larger CAS or new
state average is inferred solely from this comparison.

## What is complete and what remains

The planned continuation diagnostic is complete; it need not be rerun. The
near initialization question is answered differently for the two channels.
No numerical convergence failure or retry condition is reported in this panel.

Molpro's `DIAB` aligns orbitals within the active space; it does not constrain
which orbitals enter that space. Therefore its presence is not an acceptance
criterion for orbital/state continuity. See the
[MULTI documentation](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#diabatic_orbitals).

That read-only review has now been completed; see the linked native audit.
The originally requested next step was a review of existing full native files, using
the [existing evidence export](README.md): initial/final orbital sets,
multiplicity-labeled CI vectors and energies, and predecessor checkpoint
provenance. That review should precede any further choice of space, state
weights or constrained optimization. Nothing in this report stages or launches
another calculation, changes production grids or publishes a model.

The [twelve-request inventory](../r7_compact_resolution_20261025/resolution_audit.md)
remains applicable: the two cases above are unresolved, and the HCSSH/H donor
choice still requires an explicit conclusion. The other endpoint assessments,
five previously published models and CH3S selections are unchanged. All-at-once
publication of the remaining selections has not occurred.

A subsequent production monitor at 16:31:04 PDT is included in the linked
native audit. It confirms two previously published channels scanning and
identifies the product-graph mutation blocking the other parent.
The review monitor displays its cap of three and shared node exclusions;
that alone does not verify every production allocation or controller.
