# Native verification staged for the retained R8 channels

Status: **staged, not submitted; no new model approved; R8 is not launch ready.**
This is a finite three-phase campaign. Its purpose is to resolve the remaining
scientific evidence gaps, not to turn normal termination into acceptance.

The all-channel archive has SHA-256
`15b6ef48db20e9aa97b901a3e15759003728bc3bea585865cffcef4c90f3a256`.
All 125,793 members (24,123,802,271 bytes) passed local size/hash verification,
the embedded and sidecar manifests agree, and no changed-source or export issue
was reported. `production_native_audit.json.gz` reconciles native inputs,
outputs, geometry graphs, CI/orbital tables and root energies by channel.
The archive intentionally omitted 9,345 intermediate WFU copies; it is not a
claim that every historical scratch file was preserved.

## Coverage and reuse

`plan.json` lists every one of the 34 known channels, exact basis-specific
checkpoint hashes, source files, atom maps, state averages, reactive roots,
trial indices, and the purpose of each comparison.
`coverage.md` gives the channel-by-channel table. `source_character.json.gz`
records 52 native starting-orbital tables with file hashes and line numbers;
all requested trial indices were printed. These coefficients are qualitative
character evidence, not orbital populations or approval of the optimized result.
`pre_submission_audit.json` maps eight unresolved scientific questions to tests
for every verification case, including the extra CH+HSSH and H-exit controls.
Its outcome is **coverage checked; native outcomes pending**. It cannot certify
that all 29 models will pass and receive immediate promotion. Thirteen separate
L1 shallow-well decisions and final model/seed binding also remain explicit.

- **29 channels remain in scope:** 16 retained by the current L1 audit and 13
  awaiting matched-ZPE shallow-well checks. Unknown channels are not discarded.
- The original CH3+S, CH3+S2 and CH3S+S models retain their previous native
  chemical, state, expansion and orientation evidence. They receive no repeated
  orbital panel. Their existing review documents state the quantitative limits.
- The other **26 request identities** are covered by **25 distinct cases**.
  The two shift labels for the identical HCS+H2S product endpoint share one set
  of calculations; their identical basis seeds and state definitions were
  checked. Other stereoisomers and different product geometries remain distinct.
- Five channels excluded by the saved L1 results receive no new calculation.
- Nineteen matched, completed CASPT2-F12 probes are reused directly: three for
  the CH2+SH expanded trial and sixteen for the CH2S+SH inner-approach trials.
  Their wavefunctions seed new return/angular branches. No old failed branch is
  relabeled successful.
- Four completed production geometries are reused for CH2+SH and CH+H2S.
  Their mid-distance checks use the explicitly recorded nearest completed grid
  point. The existing production radial scans remain evidence; this panel does
  not repeat those full scans. New fixed comparisons use the original common
  diagnostic fragment geometries, so relaxed production energies cannot be
  silently substituted for them.

## Exact remaining work

| Phase | Maximum new jobs | Purpose |
| --- | ---: | --- |
| Fixed geometries | 1,149 | Matched compact/control CASPT2-F12, outgoing and return paths, angular tests and separated references, in aDZ-F12 and aTZ-F12 |
| Conditional relaxation | 46 | Near and intermediate separation, only the interfragment anchor constrained; four other geometries reused |
| Conditional relaxed-geometry single points | 208 | Every compact/control trial evaluated at the same final geometry in both bases |

The 1,149 fixed jobs are rendered locally in 112 configurations. Phase-two
workflow specifications and phase-three input builders are staged now, but bind
the actual completed checkpoint hashes and final geometries only when those
prerequisites exist. They are not fictitious native results or guessed future
geometries. The launcher updates the existing pool to a maximum of 23 jobs, preserving its
node exclusions and recovery controls. One controller drives all three phases under the
existing detached **23-job shared cap**, with eight processors and exclusive
allocations. General recon/ANL drivers own submission and bounded recovery.
This is an upper bound; the existing node exclusions and actual Slurm
availability still apply. A capacity-only change adopts the same input files,
job records and completed work; a scientific change is refused.
Interrupted submission journals stop for reconciliation rather than duplicate
jobs. Failed scientific prerequisites are recorded; independent cases continue.
Repeated invocation adopts the existing calculations and graph identities.

All main comparisons preserve the chosen spin/root averages and reactive-flux
definitions. Correlation uses CASPT2-F12 with shift 0.1, IPEA 0.25, a common
state-averaged reference, and D/T extrapolation exponent 3.7. Each model and
orientation has its own separated reference. No energy window chooses orbitals.
The orbitals are selected from hash-bound native records and then optimized.

The regular fixed grid is the own near geometry, an intermediate point, 6, 20
and 30 angstrom, followed back to 6 and the near geometry. Two independent
fragment orientations at the intermediate separation have separate 30 angstrom
references. An atomic fragment is not redundantly rotated. Previously completed
forward points are reused, while backward tests remain new by design.

## Chemical questions being tested

| Cases in plan | Preserved candidate | Additional evidence sought |
| --- | --- | --- |
| 1: CH2+SH | CAS(5,4), 2Q/4D support | SH sigma/antibonding pair in CAS(7,6), matched angular and relaxed potentials |
| 2: CH+H2S | CAS(5,4), 4D | H2S antibonding direction alone, and its occupied partner: CAS(5,5)/CAS(7,6) |
| 6: HCS2+H2 | CAS(7,6), 2D | Compare CAS(5,4) without the identified H2 sigma/antibonding pair; establish whether inclusion is necessary |
| 7: CH2S+SH inner approach | CAS(7,6), 2D | Reuse existing correlated radial evidence; add orientations, return and relaxation checks against CAS(5,4) |
| 8, 10, 13: CH2S/HCSH+SH approaches | CAS(5,4), 2D | Relevant additional sigma pair or occupied sulfur direction; reject replacement of an SH pi component |
| 9, 21: CH2+HSS inner approaches | CAS(5,4), Q/D | Omitted occupied sulfur direction in CAS(7,5), preserving the necessary carbon frontier span |
| 11, 12, 24: CH2+HSS outer approaches | CAS(3,3), Q/D | Actual sulfur donor in CAS(5,4), including bound-S2 source 18 for case 24 |
| 14, 15, 16, 29: HCSSH+H outer exits | CAS(1,1), D1 | CAS(5,5) occupied-plus-external pair trial and independent MRCI+Q(1,1); not the previous ineffective one-hole CAS(3,2) test |
| 17, 18: HCSH+SH stereochemical approaches | CAS(5,4), 2D | Own-mapping CAS(7,5) occupied union, orbital retention and correlated interaction energies |
| 19: CH2+HSS | CAS(5,4), 2Q/2D | Molecular C-S sigma/antibonding pair in CAS(7,6) |
| 20: CH2SH+S | CAS(5,4), 3Q/3D | Bound-sulfur donor in CAS(7,5), retaining the complete departing atomic 3p shell |
| 22: stereochemical CH2SH+S | CAS(7,5), 3Q/3D | Next molecular occupied direction in CAS(9,6); retain established donor and all atomic 3p components |
| 25/26, 33: HCS+H2S outer exits | CAS(3,3), D1 | Actual H2S donor 19 in optimized CAS(5,4), rather than the previously tested HCS occupied direction |
| 27: CH+HSSH, unresolved | Proposed CAS(5,4), 2D | Explicit near source 18,20,21,22; matched CAS(3,3) and CAS(7,5); independent inward starts from both far HSSH occupied combinations |
| 28: HCSH+SH | CAS(5,4), 2D | Optimized occupied union with actual correlated potential; earlier frontier-substituted solutions cannot count as convergence |

Case 27 includes **24 independent inward jobs**. At separation, source 18 and
19 are mixed HSSH sulfur combinations; neither is mislabeled as a localized
bound-sulfur donor. Both starts are transported inward to test branch dependence.
No orbitals are frozen to force the desired answer.

For the four H exits, the two external starting orbitals are explicitly **trial
virtual directions**. Their ordinal positions do not prove valence character.
The exported AO coefficients, full XML and optimized subspace must establish
that character. The MRCI+Q control includes single/double excitations into the
external space and supplies independent pair-correlation evidence. A diffuse
substitution in the small trial is a failed test, not proof that CAS(1,1) is
converged. The general explicit `open_lowest_external` partition option only
constructs these requested controls; it never chooses production orbitals.

## Scientific acceptance

`analyze.py` checks native file identities and exports all MULTI state/CI/orbital
sections, labelled correlated energies, same-geometry complete-MO subspace
comparisons, own-asymptote interaction energies, return differences, root gaps,
20-to-30 angstrom changes, angular references and D/T extrapolations. It never
writes a selected model or changes a request status.

Each final choice requires chemically complete fragment orbitals and spin/spatial
states; stable optimized active character; consistent incoming/outgoing branches;
correct separated fragment/state correlation; unchanged fragment connectivity
in relaxed checks; and adequate matched potential/gap agreement against the
relevant control. Occupations near 2/0 and total-energy lowering alone cannot
establish minimality. Compare the coordinate-dependent changes and CI character.
There is no new universal numeric threshold that overrides the prior scientific
review: report residuals explicitly, including any larger than those established
for the original CH3+S reference.

A finite panel cannot guarantee that a proposed model will pass. A failure stays
unresolved and must not be hidden by an automatic larger CAS, root relabeling,
changed product state, or declaration that every model is approved. These are
all currently identified verification tasks; unexpected native or scientific
failures remain possible.

## Git-based handoff

Use the existing `composite` checkout and virtual environment on Blodgett.
The prior standalone-bundle instructions are superseded. No SCP of a runtime
bundle or alternate installation is needed. This panel is independent of the
unfinished examples-directory reorganization.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
test "$(git branch --show-current)" = composite
git pull --ff-only origin composite
bash examples/anl/ch3ss_anl1_f12/orbital_reviews/verification_20261027/run_blodgett.sh
BASH
```

The launcher sets the existing detached pool to at most 23 jobs, stages the
checks, then starts its controller. It prints the PID and log path. An existing
controller prevents a second owner; completed/submitted job identities survive
repeated staging and a capacity-only plan revision.

When the controller finishes, gather its native evidence with:

```bash
cd "$HOME/KinBot"
bash examples/anl/ch3ss_anl1_f12/orbital_reviews/verification_20261027/collect_blodgett.sh
```

The collector includes new inputs, outputs, XML, WFUs, optimizer trajectories
and execution/attempt records. Prior evidence remains bound to the already
transferred archive and original hashes; another all-R7 export is not needed.
The frozen `prior_selection_ledger.json` and `unresolved_ch_hssh.json` are the
provenance inputs used to generate this plan. Neither is a new model approval.

Local tests exercise the complete three-phase orchestration, immutable sources,
source-hash rejection, exact donor swaps, genuine external pair controls,
shared-pool ownership, failed prerequisites, interrupted submission journals and
repeated invocation. They are simulated-site tests, not native Molpro successes.
See `validation.json` for the exact tested command and result.
