# R8 discovery and population audit — 2026-11-05

Scope: the supplied November 5, 07:48:50 PST monitor and local source/native
evidence. No native QC jobs were launched, stopped or retried for this audit.
The current full native export from the restarted run is not present locally;
the available seed archives predate these new HIR and homolysis results.

## Corrected in this update

- The additional homolysis checks and their terminal-inconclusive export hold
  were introduced on this branch, not inherited from master. Completed natural
  nonfindings no longer hold a worker or combined MESS export. Evidence remains
  intact and the capture-model assumption is stated in the generated model.
- An accepted tight path with the same configured endpoints supersedes the
  homolysis placeholder. Companion selection now agrees with the existing
  final writer, including rearrangement paths. A failed TS guess cannot replace
  a channel; a different electronic/stereochemical product is not a match.
- Unfinished companions still finish normally. Refresh adopts existing native
  work; it does not reopen exhausted searches or certify barrierlessness.
- Completed reverse campaigns with no retained reactions no longer need a
  fictitious integration step. All requested motifs and current terminal
  checkpoints must be covered. Accepted or incomplete reverse results remain
  distinct, and cannot disappear from the final model silently.
- `run_seed_audit --ordinary-native` can capture every ordinary worker and
  reverse worker, including conformers, HIR, models, reuse receipts and MESS.
  It starts no calculations and records changing/missing evidence explicitly.

## Current run: all identified outstanding checks

| Item | Evidence and disposition |
|---|---|
| Homolysis holds | The supplied monitor lists 7 CH3S and 34 CH3SS per-search holds. These are not 41 failed channels. An ordinary refresh applies the corrected terminal-nonfinding policy. Still-active searches remain active. |
| Live characterization | Fourteen allocations were running HIR work. CH3S had three ready workers; CH3SS had ten ready and three exploring workers. Preserve these calculations. Duration and old nonzero exit records alone do not justify retrying them. |
| Reverse abstractions | Automatic integration of accepted reverse saddles/entrance complexes into combined MESS is **not implemented**. Inspect the current `reverse_search_result.json` files. If any new saddle survives, match its native entrance/product evidence and integrate it before calling the revised combined PES complete. If all searches settle without accepted reactions, the new no-result exemption applies. |
| w7/w14/w15 | They have not been proven to be duplicate minima or a single thermodynamic population by the current evidence. The earlier native low-barrier TSs and modes are documented in `SOFT_MOTION_REASSESSMENT.md`. Compare the new selected geometries and HIR paths with those TS/IRCs on the same energy basis. No arbitrary rotamer merger is implemented. |
| w3 comparison | Evaluate the same soft-motion/HIR criteria; absence of a second discovered minimum is not proof of a different physical mechanism. No special-case conclusion is embedded in code. |
| rxn20 and mirror counting | Sulfur-axis selection and harmonic-midpoint logic are implemented. Whether the newly measured S–S rotor covers the relevant mirror domain, and whether that changes the counting factor, still requires the completed native rotor/model output. Detecting an axis alone does not establish coverage. |
| Shared minima and downstream work | Atom-mapped L2/HIR reuse is implemented and tested. It does not merge wells or population factors. Current reuse receipts and final conformer geometries are needed to verify which remote requests shared work; preserve independently running calculations even if later equivalence is found. |
| Final combined model | After these native outcomes are available, replay the owning writers and check all retained connections, mode replacement, population factors and any inverse/reverse routes. If a population merge is justified, every distinct entrance/exit must be remapped and internal flux must not be counted twice. The present monitor cannot certify this. |

Natural rejected saddle searches and historical QC errors are not added to a
retry queue by this audit. A concrete numerical/software defect needs its
specific native attempt evidence and applicable bounded recovery.

## Handoff

Use an ordinary `run_refresh RUN --wait-idle 600 --apply` after updating code.
Do not repeat `--check-all-homolysis` or `--reassess-thermochemistry` solely for
this policy correction: the previous restart already requested those searches
and characterization. Preserve the 20-job pool and existing native allocations.

Local regression coverage includes full simulated-site controller replacement,
submission/staging races, decision and sampling adoption, repeated invocation,
worker-to-combined MESS writing, and actual captured structural/mode fixtures.
The clean HEAD checkout with only this patch overlaid passed 551 tests and
105 subtests (2 skipped); unrelated local edits were excluded from that check.
The optional complete captured-PES replay was then run separately against the
November 1 native handoff and passed, writing both combined models twice with
QC submission forbidden and retained inputs/databases protected. The remaining
skip needs the Linux pidfd ABI. This older capture does not contain the current
restart's new HIR or reverse-discovery outcomes.
It is not verification of a remote restart or new Gaussian convergence. The
remaining population/integration questions above are explicit limitations;
neither test success nor removing the blanket hold resolves them by assumption.
