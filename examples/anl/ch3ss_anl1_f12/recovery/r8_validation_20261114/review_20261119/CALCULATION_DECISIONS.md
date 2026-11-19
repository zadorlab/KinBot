# R8 calculation necessity decisions — 19 November 2026

**Use the B2PLYP evidence to make B2PLYP judgments now.** A failed diagnostic
CC calculation is not, by itself, a reason to spend another allocation.
This ledger supersedes the earlier recommendation to finish TS2/TS3 CC
refinement and all independent native paths before making those judgments.

**Current new-QC list: empty.** No CC retries, duplicate Gaussian IRCs, fresh
Hessians, new active-space screens, or whole-panel restart are selected. The
remaining R9 gaps below are real, but their existence does not yet define a
necessary, nonduplicating calculation. This does not approve R9 or declare the
remaining 44-case campaign finished.

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
| TS5 product threshold | Match existing isolated HC(=S)SH and H states, geometries and energy convention; do not use recaptured/loose complexes as separated products | Only a missing or invalid required fragment reference warrants a new calculation |
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
- The native owner has a finished receipt. The separate L2 owner's saved
  running record and empty log do not establish current process/queue state.
- MESS remains paused. Corrected R9 launch inputs have not passed acceptance.

## Commands for this decision

There is **no QC staging or retry command selected by this review**. In
particular, do not use a whole-panel launcher or the optional TS1-F CC retry as
a finalization step. The justified immediate operation is the consolidated
read-only capture in [RECOVERY_HANDOFF.md](RECOVERY_HANDOFF.md), after updating
the checkout. Its output resolves the remaining owner/queue uncertainty and
allows completed work to be adopted before selecting any missing jobs.

For any future proposed task, record all five items before adding it to a
staging command: exact retained channel; decision affected; evidence missing;
existing-result reuse check; bounded calculation with unchanged accepted
geometry level and electronic-state definitions. Current approved task IDs:
`[]`. This ledger has no controller-consumed approval field and cannot release
held work by itself.
