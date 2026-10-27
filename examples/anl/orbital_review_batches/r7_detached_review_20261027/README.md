# Disposition of the completed detached panels

The new native evidence changes the decisions below. The old statement that the
SH direction is still lost is superseded. Calculation completion, the physical
identity of the optimized space, and a publishable VRC surface are separate
questions; the latter is not established for every channel by this panel.

The verified archive is
`examples/anl/tools/detached-orbitals-20261027T040223Z/orbital-evidence.tar.gz`,
SHA-256 `e4c5bbe1dc5346b95db594f41d9d59319173bbc3fc8a0790892adaba2bcfe2a8`.
Its 1,742 files have no collection issues or changed-source flags. All 68 jobs
completed: **44 optimized CASSCF calculations and 24 SCF character calculations**.
The latter cover six requests and cannot constitute completed CAS comparisons.
The 44 CASSCF calculations cover six other requests, with the two HCS/H2S
reaction labels sharing their verified product endpoint.

[audit.py](audit.py) independently checks every job's output, XML and WFU against
its campaign manifest, validates native CASSCF dimensions/states/energies,
retains CI and printed orbital character, and calculates 28 same-geometry MO
comparisons. [audit.json](audit.json) records the results. The maximum MO
reconstruction error is 7.49e-15; maximum orthogonality error is 1.75e-12.
No cross-geometry MO overlap or squared-AO-population estimate is used.

## Decisions and remaining questions

| Question tested | Concrete result | Disposition |
| --- | --- | --- |
| Does six-state averaging retain the SH pair in C(SH2)+SH? | Yes, in both CAS(5,4) and CAS(9,6), at both endpoints and both bases. | **Close the original SH-orbital-loss issue.** Do not repeat the exhausted two-doublet trials or enlarge CAS merely to restore SH. |
| Can the requested singlet-C(SH2)+SH asymptote be assigned to D1/D2 after that optimization? | No. At 30 Å, D1/D2 have three-open-shell character and coincide with Q1/Q2; the singlet-fragment pair is D3/D4. | **Reject fixed D1/D2 as the requested asymptote.** Endpoint data identify the problem but do not establish the connecting near-to-far branches. Minimum production CAS/root selection remains held. |
| Is an H2S donor needed for the new HCS+H2S exit? | The nominal donor CAS(5,4) loses the H2S donor near the complex and substitutes occupied HCS character. | **Reject this as a successful donor-control test.** Its large binding-energy change does not answer the donor-necessity question. The CAS(3,3) HCS frontier remains a candidate, not a newly approved surface. |
| Is the extra occupied direction needed for CH+HSSH C–S cleavage? | CAS(3,3) frontier agrees closely with the larger space, but the extra orbital in nominal sulfur-donor CAS(5,4) becomes a different occupied combination near the bond. | **Do not infer donor necessity from the larger space's energy.** The intended donor comparison did not remain controlled. |
| Is CAS(5,4) or CAS(7,5) needed for HCSH+SH S–S cleavage? | Both contain the separated SH pair, but the molecular pair/virtual part reorganizes appreciably in the larger far calculation. | **No size choice from the energy difference alone.** This is not a comparison of the same frontier plus one unchanged spectator. |
| Is an extra occupied donor needed for HCSSH+H4 outer separation? | Explicit CAS(1,1) and both CAS(3,2) controls have identical reported energies and the same H SOMO at all four anchors. | **Prepare CAS(1,1), D1**, with the actual near natural-21 seeds. This selects the formal radical minimum; it is not a pair-correlation or CASPT2 convergence claim. |
| Have the two new HCSH+SH parents and four `s8fc` requests completed CAS selection? | Their 24 new jobs are SCF only. They establish their own orbital assignments, including the correct reverse product endpoint for the insertion exit. | **Reconnaissance complete; CAS comparisons not run.** Calling these completed orbital reviews would overstate the panel's scope. |

## SH retention and state ordering

For the original fixed-spectrum natural SH pair 19/20, the two singular values
of projection into the new optimized active space are:

| Candidate / anchor | DZ | TZ |
| --- | --- | --- |
| CAS(5,4), near | 0.998688, 0.994779 | 0.998604, 0.993853 |
| CAS(5,4), 30 Å | 1.000000, 1.000000 | 1.000000, 1.000000 |
| CAS(9,6), near | 0.999452, 0.998962 | 0.999449, 0.998893 |
| CAS(9,6), 30 Å | 1.000000, 1.000000 | 1.000000, 1.000000 |

The new natural 19/20 have SH p character; 21/22 carry the carbon frontier.
Thus the earlier inactive SH substitution is absent. The source and optimized
carbon directions also relax, so the result does not assert equality of the
entire old and new active spans.

At 30 Å in CAS(9,6), the dominant D3/D4 configurations are `22a220` and
`222a20`: the unpaired electron is in the SH pair and the carbon pair is closed.
D1/D2 instead have `22b2aa`, `222baa` and their spin-coupled partners. The
corresponding quartet configurations are `22a2aa` and `222aaa`. Together with
D1=Q1 and D2=Q2 to printed precision, these identify the competing
triplet-fragment branch. CAS(5,4) gives the corresponding four-orbital strings.
The near states mix these configurations. Endpoint root numbers alone cannot
assign a continuous requested singlet-fragment surface.

This closes a specific old failure and establishes a different, narrower
remaining task: correlate the identified states along separation before
publishing a selected-root model. Neither a larger active space nor another
two-endpoint calculation by itself establishes that correlation. These gaps
are state-averaged CASSCF diagnostics, not independently optimized fragment
singlet–triplet energies.

## Why the donor energy differences are not decisions

For HCS+H2S near the complex, natural orbital 19 of CAS(3,3) is the inactive
H2S lone pair. Its squared projection into the nominal donor CAS(5,4) active
space is only **6.64e-7 (DZ), 7.23e-7 (TZ)**; projection into the latter's
closed space is **0.99999826/0.99999821**. Added active orbital 19 instead has
**0.82123/0.81959** squared overlap with the smaller calculation's occupied
HCS orbital 18. At 30 Å, the H2S lone pair is active. Consequently the apparent
near-minus-far lowering of **8.714/8.536 kcal/mol** compares different added
orbital identities. It must not be described as a measured H2S-donor binding
contribution or as evidence that CAS(5,4) is sufficient.

For CH+HSSH, the three smaller active directions project into the larger
active space with minimum singular value **0.999803/0.999788** near, and unity
far. But the larger calculation's extra near orbital is predominantly a
combination of the smaller calculation's occupied 12–17, rather than its
S-rich occupied 19. The **1.472/1.525 kcal/mol** relative lowering cannot be
used as a clean sulfur-donor increment.

For HCSH+SH, the minimum singular value for the smaller active space inside
the larger one is **0.996804/0.996842** near but **0.765231/0.756653** far.
The two separated SH directions remain identical; the substantial change is
on the other fragment. The relative change is **+1.960/+1.867 kcal/mol**, and
does not provide an isolated extra-sigma-orbital contribution. Choosing a CAS
solely by the lower total energy would not satisfy the requested minimal
chemistry criterion.

These are failed selection hypotheses, not failed executions. A further test
must keep the chemically named direction present and distinguish its effect
from occupied/virtual substitution. Re-running the same converged decks is
not a resolution procedure. No replacement panel is launched by this audit.

## HCSSH+H4 minimum

All three candidates give these D1 energies (Eh):

| Basis | L2 exit geometry | 30 Å |
| --- | ---: | ---: |
| avdz-f12 | -834.52984915 | -834.53125888 |
| avtz-f12 | -834.54547995 | -834.54691530 |

The explicit CAS(1,1) CI vector is `a = 1`; both CAS(3,2) vectors are
`2a = 1` at printed precision. The H4 SOMO's squared overlap between candidates
is at least **0.9999999999993**. It retains H4 1s character near and far.
There is no reason in these results to keep either fully occupied donor as
active for this single-reference description. The
[prepared selection](selections/new_hcssh_h_exit.json) retains one doublet and
the actual reviewed natural checkpoint for each basis. It is not installed
in `_r7` or `_r8`.

The CAS(3,2) controls are one-hole spaces and give a single determinant after
orbital optimization. Their agreement must not be advertised as an independent
test of additional pair correlation. Dynamic correlation, intermediate scan
points and production state/geometry checks remain part of the production
calculation, as for the earlier prepared compact models.

## Reproduction and scope

After safely extracting the two verified archives:

```bash
PYTHONPATH="$PWD" .venv/bin/python -B \
  examples/anl/orbital_review_batches/r7_detached_review_20261027/audit.py \
  /path/to/detached-export /path/to/root-spectrum-export \
  --output /tmp/detached-review-audit.json
```

The analysis verifies native completed calculations, not simulated chemistry.
It changes no run, scheduler state, or production model. `_r7` remains retired.
The general `_r8` seeding and outage-recovery work is still unfinished; this
scientific audit does not complete those features or authorize an `_r8` launch.
