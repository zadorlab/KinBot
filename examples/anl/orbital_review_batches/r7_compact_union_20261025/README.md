# Native review, driver fix, and ten targeted probes — 25 October 2026

**Completion update, 17:25:59 PDT:** all ten compact-union jobs are complete.
CH2/HSS supports CAS(5,4); C(SH2)/SH still shows an inactive SH-p direction.
See the [new result audit](union_audit.md) for the evidence and eleven prepared
selections. The combined publication remains held; do not treat the historical
probe hypotheses below as accepted models.

**Follow-up after the full XML comparison:** the last canonical-seed trial
selected occupied sigma character instead of the identified missing SH direction.
That [four-job natural-union panel](../r7_natural_union_20261025/README.md) is now
complete and archived. It did not resolve the last channel. Use the
[26 October diagnostic handoff](../r7_root_spectrum_20261026/README.md).
The ten-job launch below is historical and need not be repeated.

The 14-job continuation is complete. The transferred native archive has been
verified, and the VRC endpoint-discovery bug is fixed locally. Two active-space
questions remain; this folder prepares **ten CASSCF calculations for those two
questions only**. It does not publish a production model. Final publication of
all twelve remaining requests and the required controller refresh follow native
review of this finite panel. The existing controllers therefore retain their
current code until that final refresh.

## What the driver fix changes

`VTS.calculate_correction_potentials()` discovers fragment coordinates before
model selection and again before dispatch. The old vdW branch called
`make_extra_bond()`, which modified the endpoint bond graph. A later discovery
then found one connected fragment and stopped. The replacement computes the
closest cross-fragment atom pair directly from the mapped endpoint coordinates,
without adding a bond. It also respects atom-map reordering during product
matching. Repeated calls and reentry preserve the endpoint graph and geometry.

The two already published product exits retain anchors **[6,4] and [5,4]** in
the regression using their recorded L2 coordinates. Their orbital selections,
correction grids and COM sampling grids remain unchanged. The two other
previously published channels are already scanning in the supplied monitor.

No earlier PES search, IRC, L2 optimization, ANL calculation or CBH arithmetic
retry is indicated by this defect. No accepted native energy or geometry was
changed. The affected operation is VRC preparation; after loading fixed code,
the owning worker reconstructs its reaction objects from retained results and
repeats that preparation. A genuinely connected endpoint still fails the
fragment-identity check; the code does not cut a bond to force it through.

**Final refresh requirements:** load this code into the PES/KB controllers,
preserve live Slurm jobs and their adoption records, preserve policy/exclusions
and operator rejections, publish/reconcile the reviewed models through their
owners, and inspect actual scan creation. A case-level applied receipt alone
does not establish success. Check all three channels under
`791472110370000000002`, both channels under
`791671580700000000002-s36f1e2b924862f12`, all twelve newly reviewed channels,
and adoption of existing RotdPy campaigns. Any further failure is inspected
from native evidence before a targeted retry; there is no blanket rerun.

## Native evidence and the two remaining questions

The archive is retained locally at
`examples/anl/tools/compact-continuation-20261025T232718Z/orbital-evidence.tar.gz`.
All **4,056 files** matched their manifest hashes, including the WFU files.
Archive SHA-256:
`7db5f41aee419391295616af5fd029ed02ab36ff93f2a013548f2f611e8a63f1`.
The extracted text/XML was reviewed without changing source files or checkpoints.
See [provenance](evidence_provenance.json), the
[14-job results](../r7_compact_continuation_20261025/continuation_audit.md), and
[60-job inventory](../r7_compact_resolution_20261025/resolution_audit.md).

**C(SH2) + SH product exit.** The two CAS(7,5) solutions differ by approximately
20 kcal/mol at the same L2 complex geometry. A complete same-geometry MO
comparison now distinguishes an actual subspace difference from a misleading
list of leading AO coefficients:

| Singular values of active-space overlap | DZ | TZ |
| --- | --- | --- |
| Natural 18–22 in both solutions | .99949, .99839, .99627, .94103, **.02547** | .99939, .99831, .99614, .93923, **.02393** |
| Include inactive natural 17 with 18–22 | .99956, .99949, .99839, .99743, .99354, **.94102** | .99956, .99940, .99831, .99728, .99358, **.93922** |

The inward solution's inactive carbon-p orbital 17 has squared overlap
0.93970/0.94344 with the independent solution's active 20. Conversely, the
inward SH pair spans much of the independent solution's inactive SH-p orbital
17. Thus one of five active directions is almost orthogonal between the two
solutions. Adding one occupied direction accommodates that exchange; residual
orbital relaxation remains. The six-orbital comparison is a hypothesis to test,
not a claim that the minimum has already been established.

[compare_orbitals.py](compare_orbitals.py) reproduces these numbers from the
full native XML, with file hashes and results in
[orbital_comparison.json](orbital_comparison.json). It verifies identical
geometry, atom order and AO basis, then reconstructs the change of basis between
the two complete orthonormal MO sets. Relative reconstruction error is below
7e-15 and orthogonality error below 1.4e-13. No truncated-AO population estimate
or cross-geometry overlap is used for these numbers.

**CH2 + HSS C–S homolysis.** CAS(3,3) inward continuation reproduces the
independent near solution. This rules out a direct far-to-near initialization
jump as the sole explanation. The near active span contains one carbon p
direction and two sulfur-rich directions, while occupied carbon/sulfur valence
character remains outside it. The three far SOMOs alone have not demonstrated
the intended frontier span through the inner region. The next comparison adds
the specifically identified occupied S2 donor, far SCF 19, while retaining
20–22. This tests CAS(5,4); it does not add an arbitrary virtual window or the
other valence bonds. The extra 2.5 Å point brackets the observed change between
2.5 and 2 Å. The large near doublet/quartet gap and averaged 1.5/1/0.5
occupations are not, by themselves, the justification for expansion.

| Channel | Existing → comparison | Explicit source SCF indices, both bases | Points per basis; state average |
| --- | --- | --- | --- |
| `791671450790040000002_hom_sci_1_2` | (3,3) → **(5,4)** | 19,20,21,22: S2 occupied donor plus CH2 frontier pair/HSS radical | 2, 2.5, 30 Å; one Q and one D, total weights .5/.5 |
| `791871710370000000002_R_Addition_CSm_R_5_1_2` | (7,5) → **(9,6)** | 17,18,19,20,21,23: occupied S2 and C/S valence directions, both SH p components, carbon frontier/acceptor | Actual L2 product complex and 30 Å; two D roots, .5/.5 |

Each source is the same channel's already completed, hash-bound 30 Å canonical
checkpoint. Source indices are packed once; expected optimized active windows
are 19–22 for (5,4) and 17–22 for (9,6). Those windows are bookkeeping, not
predictions of the optimized orbital character. Full outputs/coefficients must
show what remains in them. Counts, changed indices and hypotheses are explicit
in the two recipes. Spin states, root counts, weights, fragment geometry,
orientation and convergence settings are retained from the preceding panel.
There are no CASPT2 jobs, production restarts, cancellations or auto-acceptance
inside this diagnostic launch.

Molpro's `DIAB` rotates within the active span; it does not prevent exchange
with inactive orbitals. See the
[MULTI documentation](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#diabatic_orbitals).
Another identical continuation would not resolve the demonstrated missing
direction. The preceding source files and all unsuccessful attempts remain
available, including the superseded VIRTORB/core-substitution evidence.

## Disposition of the other ten requests

Atom order is C1 S2 S3 H4 H5 H6. The source indices below belong to each
channel's own canonical record; natural indices belong to its completed
2140.2 checkpoint. D/Q counts are optimization support, not reactive flux.
The intended reactive contribution for these requests remains D1 only.

| Parent / suffix | Selected compact model and evidence | Further review QC |
| --- | --- | --- |
| `791871710370000000002 / hom_sci_1_2` | CH + H2SS: (5,4); source 17,18,21,22 → natural 19–22; sulfur donor plus CH frontier; two equal D. Consistent endpoint character in DZ/TZ. | None |
| Same / `hom_sci_2_3` | CHSH2 + S: (5,4); source 18,20,21,22 → natural 19–22; all three atomic-S p directions plus fragment radical; three D and three Q, equal total block weights. | None |
| Same / `r12_insertion_R_2_1_5` | HCS + H2S outer: (3,3); source 19,21,22 → natural 20–22; HCS pi/radical/pi-star, one D. Adding the tested occupied HCS valence orbital changes the endpoint relative CASSCF energy by only about .22/.21 kcal/mol. | None |
| `791671450790040000002 / r12_insertion_R_2_1_5` | CH2 + HSS outer: (3,3), source/natural 20–22; CH2 frontier pair and HSS radical, one D and one Q, equal total weights. Retained at both product-side endpoints. | None |
| Same / `r12_insertion_R_2_1_6` | Same compact contents and counts from this exit's own checkpoint and atom map. | None |
| Same / `hom_sci_2_3` | CH2S + SH homolysis: (5,4); source 18,19,21,22 → natural 19–22; CH2S pi/pi-star plus SH pair; two equal D. Alternative donor initialization converges to the same solution. | None |
| Same / `r12_insertion_R_2_3_4` | CH2S + SH outer: the same orbital types and counts from this product-complex checkpoint, with its distinct geometry. | None |
| Same / `r12_insertion_R_5_1_6` | HCSSH + H5 outer: formal minimum **(1,1)**, retain natural 21 (H SOMO); one D. The two tested (3,2) donors give the same energy and occupations exactly 2/1 at reported precision. | None |
| Same / `r12_insertion_R_6_1_2` | HCSSH + H6 outer: formal (1,1), its own natural 21 H SOMO; same two-donor evidence. | None |
| Same / `r12_insertion_R_6_1_5` | HCSSH + H6 outer: formal (1,1), its own natural 21 H SOMO; same evidence at this distinct endpoint. | None |

For the last three rows, (1,1) is an explicit compact model choice based on the
single-configuration endpoint results, not a claim that a separate (1,1)
CASSCF or correlated-potential comparison was run. The filled donor can remain
inactive and optimized; no evidence here establishes a need to make it active.
Production still performs its configured optimization and retains orbital,
root, connectivity and scan checks. All ten decisions await the combined
publication; none is installed by this document or by the launch below.

The existing and earlier five selections, evidence, weights and limitations
remain documented in [ALL_ACTIVE_SPACES.md](../r7_progress_20261024/ALL_ACTIVE_SPACES.md).
The C + SH3 reaction remains operator-rejected. No new reaction is discarded.

## Driver snapshot and remaining work

The [16:31:04 PDT monitor](production_monitor_163104.txt) reports cap **23**, policy
revision 14, **22 allocations** (12 CH3S, 10 CH3SS), zero diagnostic allocations,
and 16 applied live-control requests. Existing node exclusions remain in force.

| Area | Reported state | Next work |
| --- | --- | --- |
| CH3S exploration | Postprocessing; four wells; 9 retained and 20 omitted outer selections | Retained CBH-ANL and kinetics work continue |
| CH3S VRC | Four retained state surfaces each have 48 accepted scan entries including 30 Å | Four RotdPy campaigns running; final flux/kinetic validation remains |
| CH3S ANL | CC geometry and higher-excitation single points running; 1,115 completed tasks across the displayed graphs | Dependency-ready tasks continue through the owners; no simplification requested here |
| CH3SS exploration | **Five wells exploring and five at VRC**, 5/5 search slots | Search, product characterization and TS optimization remain; the final network is not known |
| Initial CH3SS | C–S 48 accepted entries and RotdPy running; S–S 22 accepted, 6.4 Å running | Continue scans/sampling |
| Previously published CH2SHS | C–S 1 accepted, 2.2 Å running; S–S 2 accepted, 2.4 Å running | Confirms these two reviewed models reached actual scan dispatch |
| Previously published CH2(SH)S parent | Endpoint-graph hold blocks preparation | Fixed code plus final controller refresh; verify three channels start |
| Twelve remaining requests | Ten compact decisions documented; two require this comparison | Review ten new native results, publish all twelve together, refresh/reconcile |
| CH3SS ANL/CBH | New work waits for outer L2 selection | Finalize network after exploration, then assess sharing/simplification/pruning as requested |
| MESS | Waiting; executable ready | Requires accepted thermochemistry and required kinetic fluxes; no completed solve claimed |

The RotdPy `0/116` fields count final surface output files, not individual
sample attempts. Live jobs are visible under Surface_38/46/83. This snapshot
alone is insufficient to declare either a stall or completed converged flux.
The old optimizer-limit and zero-imaginary-mode TS failures remain rejected
attempts. No current failure evidence justifies changing their validation,
cancelling unrelated work or rerunning them indiscriminately.

Reported CH3S expensive allocations have elapsed roughly 24 hours to 3 days
8 hours; elapsed time alone does not determine their remaining work. Related
compact near probes took 16–22 seconds in DZ and 155–192 seconds in TZ. New
spaces and scheduler admission can change that timing. No defensible finish
time for CH3SS exploration follows from a single snapshot while new wells
remain active.

The existing [ANL live controls and thermal contract](../../../../docs/anl_live_controls.md)
remain applicable: retained high-level geometry/harmonics in MESS, individually
projected TZ/QZ anharmonic shifts, separate total-ZPVE arithmetic, L2 rotor
data/geometry, harmonic-only vdW binding differences with independent product
thermochemistry, and dependency-aware live simplification/sharing. This fix
does not change those choices. Local regressions do not establish that every
native ANL result or final MESS file has already been scientifically accepted.

## Run the finite review panel on Blodgett

Use the user's SSH session. This launches only the ten new diagnostics, up to
three concurrently within the existing shared cap23. Completed reviews and
production allocations remain untouched. Preflight checks every seed hash and
deck before staging anything. `--resume` and the campaign lock prevent duplicate
staging/owners. No production controller restart is needed for this standalone
review; the final combined refresh is deliberately deferred until its review.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
test "$(git branch --show-current)" = composite
git pull --ff-only origin composite
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PLAN="$PWD/examples/anl/orbital_review_batches/r7_compact_union_20261025/screens.json"
KB_PANEL="$KB_RUN/orbital_recon/compact_union_20261025"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_RUN" "$KB_PLAN" "$KB_PANEL" --shared-pool --preflight --resume
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_RUN" "$KB_PLAN" "$KB_PANEL" --shared-pool --resume
KB_LOG="$KB_PANEL/run.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.requested_orbital_recon \
  run "$KB_PANEL" --report-failures >"$KB_LOG" 2>&1 </dev/null &
printf 'Review PID: %s\nLog: %s\n' "$!" "$KB_LOG"
BASH
```

Read-only progress and character extraction (writes only the report):

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PANEL="$KB_RUN/orbital_recon/compact_union_20261025"
KB_REPORT="$HOME/compact_union_character.$(date -u +%Y%m%dT%H%M%SZ).txt"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor "$KB_PANEL" >"$KB_REPORT"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon character "$KB_PANEL" >>"$KB_REPORT"
"$PWD/.venv/bin/python" -B -m kinbot.run_monitor "$KB_RUN" >>"$KB_REPORT"
printf 'Report: %s\n' "$KB_REPORT"
cat "$KB_REPORT"
BASH
```

After completion, export the ten new jobs with full WFU/XML/outputs and latest
driver metadata. The already verified source campaigns need not be transferred
again. This reads source files and writes a verified new export directory:

```bash
cd "$HOME/KinBot"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_EXPORT="$HOME/kinbot-evidence/compact-union-$(date -u +%Y%m%dT%H%M%SZ)"
"$PWD/.venv/bin/python" -B examples/anl/tools/export_orbital_evidence.py \
  examples/anl/orbital_review_batches/r7_compact_union_20261025/export.json \
  --repo "$PWD" --runs "$HOME/kinbot-runs" --production "$KB_RUN" \
  --output-directory "$KB_EXPORT"
"$PWD/.venv/bin/python" -B -m kinbot.run_monitor "$KB_RUN" --all --details \
  >"$KB_EXPORT/production_monitor.txt"
printf 'Evidence directory: %s\n' "$KB_EXPORT"
```

## Validation and limits

Local tests cover repeated coordinate discovery, the two published exit anchors,
both selection/dispatch passes, sampling reentry, retained files, and genuine
fragment-identity rejection. The actual ten-job plan is tested for whole-panel
preflight, native-source partitions, unchanged endpoint geometry/state weights,
bad-hash refusal, repeat staging and simulated scheduling retaining 22 other
allocations under cap23. Controller-refresh tests include submission/staging
races, live-job preservation, model/decision application, sampling adoption and
repeat invocation. MESS/vdW, ANL arithmetic/live controls, sharing and resource
policy regressions are also run. Test counts are recorded in the accompanying
handoff response. These are local software tests; a successful native retry,
remote refresh, new orbital acceptance or completed rate is not claimed.
