# R9 readiness work ledger

## Current acceptance criteria — 19 November 2026

The objective is to resolve the paper–KinBot agreements and disagreements from
native evidence, then define a scientifically justified `_r9` calculation and
kinetics model. Completion of the Blodgett validation jobs is the start of the
final evidence review. `_r9` is not ready until the gates below close.
The 19 November export has now been reviewed: see the
[full evidence review and closure sequence](R8_FINAL_EVIDENCE_REVIEW_20261119.md).
It covers all 245 native task records, the separate 44-case L2 campaign,
24 energy-exclusion decisions and 105 serialized network connections. The
remaining scientific and operational findings are explicit in that review.

**Necessity rule:** existing converged B2PLYP paths support B2PLYP connectivity
judgments even when an optional CC or redundant Gaussian check fails. No further
CC optimization/path campaign is required to make those judgments. The current
[new-calculation admission](../examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/review_20261119/CALCULATION_DECISIONS.md)
is empty; missing jobs must be selected only after retained-channel and reuse
review. R9 model closure is not the same as completing every diagnostic.

This section supersedes earlier readiness implications and unconditional
"PST → VRC" treatment recommendations in the scouting scope documents. The
historical scope tables remain inventories of the model that was calculated;
their channel labels and branching fractions are not approval of that model.

- [ ] **Paper comparison ledger:** for each disputed connection, record the
  paper claim, atom-mapped species/fragment states, native saddle and path
  evidence, matched energy reference/theory convention, and an explicit verdict:
  supported agreement, supported disagreement, or unresolved. Distinguish a
  method-dependent result from a failed calculation, incomplete path, or parser
  error. Bind each verdict to source files and attempt hashes.
- [ ] **Mechanism and bottleneck audit for every retained/candidate exit:**
  explicitly review beta scission and H loss as well as other dissociations.
  Resolve whether the evidence supports an inner conventional saddle, an outer
  loose-fragment separation, or a path containing both. A `nobar`, `hom_sci`,
  beta-scission, or PST label, a disconnected optimized endpoint, and failure to
  locate a saddle do not establish a barrierless path. Retain inner barriers
  where supported; assign PST/VRC only to an evidenced applicable outer/loose
  bottleneck. Replacing PST with VRC alone does not resolve a misclassified inner
  reaction. Record the evidence needed wherever classification remains open.
- [ ] **Stationary points and paths:** verify native geometry convergence,
  imaginary-mode character, continuous path progress and stopping criteria,
  both endpoint identities, and relevant step sensitivity. Keep raw IRC
  endpoints separate from subsequent optimization and fragment recapture.
- [ ] **Electronic/thermochemical consistency:** verify reference/state
  continuity, convergence, diagnostics, fragment multiplicities, ZPE, basis and
  triples conventions, and matched well/TS/asymptote energy references before
  comparing barriers or accepting high-level corrections.
- [ ] **Reassess the R9 scope:** after correcting topology or rate treatment,
  recompute the affected scouting/sensitivity comparisons before reusing their
  pruning recommendation. Include paper-relevant competing channels explicitly;
  previous lean/broader percentages alone cannot settle the corrected scope.
- [ ] **Exact retained-channel model map:** assign each justified bottleneck its
  thermochemistry and kinetics treatment. For selected multireference VRC work,
  bind the human-selected active space, orbital/state evidence, atom mapping,
  basis-specific assets and spin treatment to the exact channel. Reuse prior
  evidence only after checking applicability.
- [ ] **Generated kinetics input audit:** trace every retained MESS edge back
  to its accepted channel, inner/outer bottleneck and numerical inputs; verify
  energy references, degeneracy/symmetry/spin factors and reversible-pair
  consistency, and remove duplicate representations of the same pathway.
- [ ] **Independent verification pass and launch package:** perform a separate
  source-to-conclusion check without additional agents. Deliver the paper
  comparison ledger, corrected network, accepted/excluded/unresolved channel
  ledger, exact R9 inputs/reuse plan, and checks of the generated workflow and
  kinetics inputs. Keep each unresolved issue visible with its consequence and
  next evidence requirement; readiness requires closing issues that affect the
  selected R9 scope.

## Historical work ledger and scope evolution

R7 is retired and remains immutable. The user has now authorized staging all remaining necessary orbital checks
after examining the all-channel native export. Submission is a separate action;
no new diagnostic or production job was submitted locally.
Do not call the submission ready until every requirement below is verified.

- [x] Verify seed evidence, databases, terminal search history and old optimizer
  trajectories; attribute the original MRCC checkpoint (seed review report).
- [x] Durable accepted ASE geometry points and bounded dispatcher allocation
  retries, with local interruption/orchestration tests.
- [x] Periodic native candidate preservation, including hard-kill survival tests.
- [ ] General independent run import with source immutability, preserved terminal
  failures, compatible complete-result reuse, pending-submission reconciliation,
  path rebinding, interrupted-transaction recovery and repeated invocation.
- [ ] Cross-PES shared ownership under an explicit shared root, including
  simultaneous registration, completed/live work, and restart tests.
- [ ] Native checkpoint selection/restoration with phase/applicability checks;
  persistent bounded retries in all owning drivers. Native validation claims
  must be distinguished from local synthetic tests.
- [ ] Apply L1 150 kcal/mol pruning and matched-ZPE 0.5 kcal/mol vdW checks in R9,
  including reassessment of imported accepted paths without reopening failures.
- [ ] Analyze the existing closure panel; retain earlier choices and export an
  exact, evidence-linked selection for every known VRC request.
- [ ] Assemble independent R9 inputs and all basis-specific model/seed assets;
  fail submission preflight if any known required channel has no model.
- [ ] Reorganize CH3S/CH3SS examples and preserve their historical evidence and
  user edits; keep only the two case folders at examples/anl before the demo.
- [ ] Add minimal ethane ANL1-F12 demo with augmented basis, CASPT2-F12 T/Q CBS
  correction and CCSD(T)-F12b T/Q CBS spin splitting, ten exclusive nodes/eight
  processors; general opt-in single-state CAS(2,2) discovery for applicable
  simple doublet-radical pairs, with native evidence gates.
- [ ] Run complete simulated-site orchestration and example/preflight validation;
  explain each final CAS choice and its evidence in chat.

Closure evidence is local: archive SHA-256
`c213d108b7c27721684eed9921d1fb9d2a86bc8ad7c3e0ba4933cee60b83ebef`.
All 4,012 exported files passed their manifest hashes, with no collection issues;
all 90 closure probes have native completion. Scientific decisions are being
recorded from CI, orbital character and state correlations; completion alone is
not approval. The subsequent all-channel export is verified and the remaining verification
calculations are staged under `verification_20261027`; completion and scientific
approval remain separate.

The proposed CH + HSSH CAS(5,4) approval was withdrawn: the fixed-CI controls
do not validate the exact optimized partition. It is stored only as a blocked
proposal, with no selected model. R9 is not launch ready. The other assessments
must be checked individually against their native optimized orbital evidence;
the number of completed diagnostics cannot satisfy this requirement.

2026-10-27: All-channel native archive `15b6ef48db20e9aa97b901a3e15759003728bc3bea585865cffcef4c90f3a256` verified file by file (125,793 files, no issues). The staged verification inventory covers 29 retained/pending-L1 VRC requests, reuses three initial reviewed models, excludes five channels, aliases the identical F/R exit endpoint, and preserves the prior choices as candidates. It does not mark R9 ready.

2026-10-27 scope change: high-level production is deferred to R9. Retained channels will be chosen explicitly after the independent R8 L1/L2/PST scouting pass. Existing orbital verification evidence remains required only for those chosen channels.
