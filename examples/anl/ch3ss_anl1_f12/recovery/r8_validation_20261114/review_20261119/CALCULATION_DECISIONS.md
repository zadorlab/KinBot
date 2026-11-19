# R8 calculation necessity decisions — 19 November 2026

**Use the B2PLYP evidence to make B2PLYP judgments now.** A failed diagnostic
CC calculation is not, by itself, a reason to spend another allocation.
This ledger supersedes the earlier recommendation to finish TS2/TS3 CC
refinement and all independent native paths before making those judgments.

**Selected additional work: five initial allocations.** The exact retained-channel
review now admits four B2PLYP relaxed scans and the single missing reverse path
for `ch3ss_rxn_1_1`. No failed CC diagnostic is restarted. The user's current
Blodgett queue/process check returned no output; the owners still reconcile
live allocations if the site changes before launch.

## Exact admitted calculations

All new scans use B2PLYP-D3(BJ)/cc-pVTZ, the preserved neutral-doublet states,
UltraFine integration and SCF XQC/Conver=8. Each constrained optimization has a
150-cycle ceiling. The dispatcher retains input hashes, job identities and
native outcomes. The existing shared 20-job policy and node exclusions apply.

| Job | Missing evidence and reuse check | Bounded calculation | Decision affected |
|---|---|---|---|
| `w3_cs_out/scan` | Paper-seeded TS3 and native `12_2` connect W2 to encounters, not W3. The old W3 C–S(H) scan is L1 and drops abruptly after its maximum | Start at accepted W3; extend C–S(H) from 1.830805 to 4.230805 Å in 24 increments | Resolve the W3 β-scission inner branch and identify a possible saddle seed |
| `complex_a_cs_in/scan` | Valid B2PLYP encounter A exists, but its W3 approach is untraced | Contract its C···S(H) distance from 3.873747 to 1.85 Å in 21 increments | Identify whether A reaches W3, W2, or another branch |
| `complex_b_cs_in/scan` | Valid B2PLYP encounter B exists; the paper-seeded path instead reaches W2 | Contract C···S(H) from 2.808133 to 1.85 Å in 10 increments | Resolve the independent B-side approach and hysteresis/branch switching |
| `w4_ch_out/scan` | Direct CH3SS H loss is retained; its old L1 scan/guessed saddles do not settle the B2PLYP bottleneck | Extend C–H from 1.087994 to 5.087994 Å in 40 increments | Test direct CH3SS → CH2SS + H without borrowing TS5 from another parent |
| `ch3ss_rxn_1_1/saddle_IRC_R` | Retained source `791471490990000000002_intra_H_migration_2_4`; B2PLYP saddle is 70.730 kcal/mol above W4. F already converged in 95 steps; R has no allocation/result | Existing Sella owner, original frozen settings, this case only; stop after path/endpoint optimization | Resolve the other side of this retained competing H-migration route |

The two encounter-start scans are different initial branches, not repeated
forward paths from W3. The old L1 W3 trace falls from −836.2571703 to
−836.2802785 Eh between 2.930805 and 3.030805 Å. Following both encounter
approaches makes the potential branch change observable. Source coordinates,
atom maps, hashes and old scan points are frozen in
[paper_closure_inputs.json](../paper_closure_inputs.json).

## Closure decisions after these results

| Observation | Required next action | Acceptance evidence |
|---|---|---|
| Continuous scan maximum on the intended branch | Refine that geometry as an unconstrained B2PLYP saddle, then run its two path directions through the existing owner | Small stationary gradient, one reaction-mode imaginary frequency, both atom-mapped connections |
| Scan changes H attachments or electronic branch | Identify the first switch and compare the opposing trace before selecting a local constrained follow-up | A continuous intended branch; an energy jump is not a transition-state energy |
| Monotonic intended branch reaches the existing outer domain | Match the inner/outer geometry, states and asymptote before assigning the retained bottleneck treatment | Continuous sampled path and an applicable outer VRC surface; monotonic finite sampling alone is not a universal barrierlessness proof |
| Missing reverse IRC finishes | Map both endpoints to the existing spectrum/geometry library | Reuse matching spectra; calculate a frequency only for a genuinely new retained endpoint |
| Any task fails or exhausts its bound | Diagnose the native output and exact attempt before admitting another job | A retry must address an evidenced failure class and still change a consequential decision |

These are the initial jobs whose inputs can be fixed from current evidence.
A saddle geometry or branch-specific follow-up cannot be selected honestly
before its supporting scan exists. The owners preserve all results for that
next decision; they do not automatically promote a constrained maximum into MESS.

## Reused TS5 and fragment threshold ledger

Seven retained well/fragment records were independently matched to native
Gaussian inputs/outputs: theory, charge/multiplicity, correlated B2PLYP energy,
ZPE, frequencies, normal termination, and geometry distance matrices. Maximum
geometry discrepancy is below 0.000001 Å. See
[reused_fragment_verification.json](reused_fragment_verification.json).

| Quantity | kcal/mol relative to W4 |
|---|---:|
| W3 | 21.022629 |
| HC(=S)SH singlet + H doublet | 60.080360 |
| Paper-like B2PLYP TS5 | 64.916106 |
| CH2SS singlet + H doublet | 73.003054 |
| CH2S singlet + SH doublet | 43.875989 |

Convention: B2PLYP-D3(BJ)/cc-pVTZ electronic energies plus unscaled harmonic
ZPE, with the recorded spin-free states. TS5 lies **4.835746 kcal/mol above**
HC(=S)SH + H in this convention. Existing separate fragments supply that
threshold; the recaptured optimized IRC endpoint is not used for it. No new
fragment QC or CC path is selected. This ledger does not silently apply a
spin–orbit, rotor or anharmonic correction.

## Paper comparison: decisions available now

Level: B2PLYP-D3(BJ)/cc-pVTZ. F/R numbers are Sella path steps, not optimizer
iterations from an unrelated endpoint job. The prior export records converged
paths and positive endpoint spectra for these cases. See
[prior_l2_connectivity.json](prior_l2_connectivity.json) for individual receipts.

| Object | Existing result | Decision | New QC for this decision |
|---|---|---|---|
| TS1 | W1 ↔ HCS···H2S; 241/156 steps | Expected abstraction connection supported | None |
| TS2 paper seed | W2 arrangements on both sides; 158/119 | Nominal paper TS2 connection not reproduced by this B2PLYP search | None |
| Native `54_1` | W2 ↔ CH2S···SH complex A; 158/50 | A native B2PLYP route supplies W2–complex connectivity | None |
| TS3 paper seed | W2 ↔ CH2S···SH complex B; 208/91 | Nominal W3–W1 connection not recovered at B2PLYP | None to report that result |
| TS4 / native `2_1` | W4 ↔ W2; TS4 149/129 | Reuse one mapped representative | None |
| TS5 | W3 parent side; raw Gaussian product-side H separation at two path steps, followed by recapture during optimization | Retain an inner H-loss saddle; reject barrierless-only representation of this route | None to make this assignment |
| TS6 / compact `12_1` | W2 ↔ W3, 142/117 and 131/129 | Preserve distinct saddles after geometry/mode mapping | None |
| Lower `12_2` | W2 ↔ complex B; 146/107 | Remove from the W2–W3 barrier competition; retain its actual B2PLYP connection | None |
| CH2SH H loss | `scale_1p6` gives CH2SH ↔ CH2S···H, 80/123; `1p9` duplicates it | Include the inner saddle in the reaction model | None |
| CH3S `1_2`, `3_1` | Two loose CH2S···H arrangements | Reassign these observations rather than import their old reaction labels | None |
| CH3SS `2_2` | W2 ↔ CH2···HSS, 214/243 | Distinct from W4–W2 and CH2S···SH | None |

Paper-seeded TS2/TS3 non-recovery is a **method-specific comparison result**.
It neither requires successful CC calculations nor proves that the paper's CC
surface is wrong. TS5's recaptured endpoint is not the separated product.

## Diagnostic attempts to retire from the proposed retry list

| Work | Evidence | Disposition |
|---|---|---|
| CC TS2 optimization/curvature | Gradient about 0.04357 Eh/bohr; nonstationarity already independently confirmed | Reject this CC geometry. No further optimization, gradient or Hessian for the B2PLYP decision |
| CC TS3 electronic gate | Quadratic ROHF converged; subsequent CC did not | Keep the diagnostic history; no new CC iteration trial merely to finish this optional benchmark |
| CC compact `12_1` | Failed refinement and elevated T1/D1 diagnostics | Keep B2PLYP saddle/path evidence; no forced CC acceptance or retry |
| CC TS5 forward/reverse paths | Forward reaches a QSD limit; reverse ROHF fails | No retry needed to establish the already-supported inner H-loss saddle |
| CC W1 refinement | Incomplete failure evidence | No retry needed for B2PLYP complex identity or topology |
| TS1-F TZ/QZ corrections | TZ retry has a generated MAXIT syntax defect; QZ failure tail is incomplete | Shared syntax fix is useful software maintenance. It does not make these single points necessary for topology |
| Independent Gaussian checks of accepted Sella paths | Corrector failures or early-stop behavior | Retain diagnostics; stop repeating equivalent connectivity questions |
| Repeated `12_2` departure checks | Original/positive-time trials stop early; bounded departure hits its point limit | The earlier two-sided Sella assignment is already available; no additional departure wave |
| W2/W3/W4 and paper-like TS5/TS6 CC results | Valid stationary structures and T/Q energies | Reuse. No recalculation |
| Existing orbital reconnaissance | Human-selected spaces and native exports exist | Reuse after exact channel/state/basis binding; no blanket new screen |

Retiring a proposed retry means **no new submission**. It does not delete a
failed attempt, relabel it successful, cancel a running job, or change its
scientific thresholds. A future quantitative ANL/composite claim may need a
specific missing energy correction, but that need must be tied to the chosen
energy ledger rather than inferred from this diagnostic task's failure.

## Actual remaining R9 gaps

| Gap | First action using existing evidence | When an additional calculation becomes necessary |
|---|---|---|
| W3 → CH2S + SH inner route | Map all existing B2PLYP saddles/scans to W3 and encounter arrangements; the nominal TS3 search and `12_2` both lead elsewhere | If this direct route is retained and the inner path remains unsupported, a targeted B2PLYP parent-to-encounter investigation is needed. Outer CAS/VRC evidence cannot substitute for it |
| Direct CH3SS → CH2SS + H | Review retained parent-bound dissociation evidence; prior guessed saddles gave other chemistry or failed | If retained, establish its bottleneck at the production geometry level. Failed TS searches do not authorize PST |
| TS5 product threshold | Matched native isolated HC(=S)SH and H records; threshold 60.080360 kcal/mol above W4 in the stated B2PLYP convention | Resolved by reuse; zero new fragment calculations |
| Other PST candidates | Bind all retained entries in the 105-edge inventory to their actual inner/outer evidence | Calculate only a demonstrated missing bottleneck in the selected scope; a PST label is not an admission rule |
| Remaining 44 L2 cases | Reconcile the 13 completed directions, reuse identical validated endpoint spectra, map each candidate to retained edges | Only a missing retained direction or genuinely new endpoint spectrum; the 75 missing direction records are not automatically 75 required jobs |
| 24 L1 energy exclusions | Preserve 14 usable L2 energies; reconcile ten unresolved saddle refinements with retained scope and actual energy range | Only an exclusion needed for the final scope and unsupported by another valid bound. Do not repeat all ten failed searches automatically |
| Paper energy comparison | Use completed T/Q energies with explicit geometry, ZPE and reference conventions | A specific missing component only when needed for a retained quantitative claim; optional CC geometry replication is separate |
| SI Eckart depths | Keep literal SI and energy-consistent variants; two W1-side depths differ from the energy ledger by exactly 3 kcal/mol | No QC; this is a model-input comparison |
| Scope/branching | Fix topology, degeneracy, spin, energy distribution and observable first | Corrected kinetic sensitivity after MESS is intentionally resumed; no new QC merely to rerun old pruning percentages |
| Orbital/state binding | Bind existing WFU hashes, roots, weights, atom maps, asymptotes and domains | New orbital calculations only for a demonstrated uncovered domain/state |

An omission from a provisional model must be labeled as a scope limitation;
it must not become a claim of physical absence or negligible yield. Likewise,
a complete 105-edge inventory is not independent validation of every edge.

## Current captured status

- Native panel: 245 task records — 97 complete, 98 failed, 48 blocked, 2 staged.
  These counts include historical attempts and optional diagnostics.
- Separate remaining-L2 panel: 44 cases; 13/88 directions have path/optimized
  endpoint records; 75 lack those records in the export. None of those 13 has
  an endpoint-frequency record in this particular export. Existing matching
  endpoint spectra elsewhere must be checked before requesting new ones.
- CC thermochemistry retained: W2/W3/W4/TS5/TS6. Their QZ + DZ-harmonic-ZPE
  relative energies are 24.510, 21.091, 0, 64.629 and 60.098 kcal/mol; the
  comparable supplied SI values are 24.564, 20.197, 0, 64.643 and 61.220.
  This is a stated comparison convention, not exact protocol reproduction.
- Energy-exclusion audit: 14 usable L2 decisions span 151.24–170.98 kcal/mol;
  ten attempted saddle refinements remain unresolved. The 151.24 case is not
  a robust exclusion across a 140–155 kcal/mol study interval.
- The user reports no current jobs or matching validation owners from the
  requested queue/process check. The stale saved L2 owner record is superseded.
  A replay of the captured `1_1` case with the actual QC status reader requests
  only `saddle_IRC_R`, reusing its forward result.
- MESS remains paused. Corrected R9 launch inputs have not passed acceptance.

## Commands for this decision

Use [RECOVERY_HANDOFF.md](RECOVERY_HANDOFF.md). It stages the four frozen scans
and starts their existing dispatcher plus the selected Sella owner. The other
43 L2 cases and all failed CC records remain preserved. MESS remains paused.
The five initial allocations are deliberately scoped to missing consequential
evidence; a blanket restart would add work without a corresponding decision.
