# r7 VRC channel and orbital handoff, 2026-10-21

This records the three initial reviewed channels, the applied CH2SH selection,
and the three additional orbital selections completed on October 21. Live
refresh adds accepted channels without replacing existing models. The two
provisionally screened CH3SS cleavages still await their accepted requests.

## Bounded scan recovery and live review, October 22

The C + SH3 relaxed scan changed chemical connectivity at **C–S = 2.8 A**:
S(2)–H(4) grew from 1.436 to 2.9422 A while C(1)–H(4) became 1.1099 A.
The separated, optimized 30 A reference has S(2)–H(4) = 1.4785 A.
This is an H transfer, not a parser error. Native optimization completion
does not make that geometry a valid C + SH3 correction point. Passing the
connectivity check at 2.0–2.6 A also does not prove absence of earlier drift.

The supplied L2 network already contains the corresponding well-to-well
H-transfer TS, with a 0.689 kcal/mol forward barrier and approximately
-747 cm^-1 imaginary frequency. The outer PES retains its reverse discovery
`471230610040000000002_intra_H_migration_2_5` (84.260 kcal/mol reverse
barrier). The duplicate forward searches are omitted. This supports a
competing tight-TS rearrangement; it neither proves nor disproves a separate
C + SH3 dissociation channel. A `hom_sci` label is a channel proposal, not
proof of the complete barrierless potential.

Production behavior after the controller refresh:

1. Keep the existing 47-point scan and dividing-surface grids, 2–20 A, and
   the **separate 30 A optimized fragment reference**. No distances, active
   spaces, methods or convergence thresholds are changed by this repair.
2. If a Molpro relaxed point breaks an internal fragment bond, preserve the
   rejected graph and try that distance once with the validated separated
   fragment geometry and its orbital checkpoint. This is a separate graph,
   limited to **100 optimizer steps and at most two hours**, or the shorter
   explicit user walltime. Automatic step extensions and further allocations
   cannot multiply this retry budget. It uses the existing shared pool,
   node exclusions and resource requests; its walltime ceiling takes precedence
   over a pool-wide walltime floor. Explicit later user retries can change
   that budget.
3. On success, use the accepted retry geometry/checkpoint as the next seed.
   If the retry fails, retain a named gap and attempt subsequent distances
   from the validated fragments. Never publish the rejected point, interpolate
   across the gap, or silently substitute a frozen-fragment energy.
4. Persist each point's original geometry/orbital seed before dispatch.
   Controller refresh and later repair of an earlier gap reattach to later
   accepted/allocated graphs using their original seeds. A failed automatic
   retry is not renewed merely because the driver or input is refreshed.
5. Any connectivity failure holds the **entire affected scan for review**,
   even if the fragment retry succeeds and the grid completes. Complete
   flagged grids can finish their correction calculations; their channel
   assembly/rates remain held. Incomplete grids retain later successful
   points but cannot become a correction curve. Other surfaces/channels,
   exploration and ANL continue independently. Final MESS completion still
   requires resolution of every required channel.

Evidence is written under `vrctst/scan_progress/<leaf>/<context-hash>/` and
`vrctst/scan_review/<leaf>/scan_review.json`, alongside original native
files, rejected-attempt markers, fallback graphs, accepted artifact hashes
and orbital checkpoints. Review history survives replay. The monitor keeps
this scan-path hold distinct from an already accepted orbital-model receipt.
ROTD generation and result handoff reject a held scan even when older
correction/rate files exist.
Correction-to-scan bindings use sidecar receipts, preserving the existing
correction bytes/hash when the science is unchanged. Publication is atomic.

The C + SH3 request in r7 is expected at:

```bash
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_SCAN_REVIEW="$KB_RUN/ch3s/471430450000000000002-s06a7b4c221beb748/vrctst/scan_review/471430450000000000002-s06a7b4c221beb748_hom_sci_1_2__single_m2_r1/scan_review.json"
"$HOME/KinBot/.venv/bin/python" -m kinbot.vrc_scan_recovery status "$KB_SCAN_REVIEW"
```

After this code is installed, ordinary users can request a targeted retry
through the live owner, without source edits or another controller replacement.
For example, **after reviewing an unresolved 2.8 A point**, save a schema-1
fallback JSON such as `{"schema":1,"optimizer":{"steps":200},"resources":{"walltime":"04:00:00"}}`
and request:

```bash
"$HOME/KinBot/.venv/bin/python" -m kinbot.vrc_scan_recovery retry \
  "$KB_SCAN_REVIEW" --distance 2.8 \
  --reason "Reviewed failed point; retry with the specified optimizer budget" \
  --fallback /path/to/reviewed-fallback.json
```

This queues only the last failed attempt at that gap. Existing per-attempt
fallbacks also support calculator keywords, commands, environment/setup and
resources within their validation rules. They retain the requested scientific
identity and parser checks; changing methods, active spaces or atom identity
requires the corresponding reviewed input migration. Duplicate pending
requests and attempts to replace an accepted point are rejected.
For a finished graph whose Slurm record has been purged, a live retry requires
a matching saved allocation ID, verified stage/output provenance and a
successful queue query confirming it is absent. Missing accounting alone
does not force a restart; ambiguous or unavailable scheduler evidence blocks
the retry without cancelling or resubmitting anything.

Only after **all points are accepted and the complete path has been reviewed**,
the `approve` action with a substantive `--reason` records approval bound to
the exact evidence hash. The live owner wakes, revalidates/reuses work and
continues; no new driver is started. Approval cannot waive a missing point,
changed artifacts or unresolved orbital review. An optimizer success alone
is not a sufficient reason to approve this H-transfer-prone path.

Use the [single controller-refresh command](ANL_COMPONENTS_V2_AND_R7.md#existing-r7-native-interface-repairs-october-22)
once to install this implementation into the current r7 controllers. It
retains running Slurm work. The new CH + H2S orbital review remains a separate
requirement; this recovery does not choose an active space for that channel.

## Additional orbital-only batch, October 21 evening

The plan `examples/anl/orbital_review_batches/r7_additional_20261021/batch.json`
stages these **three screens, eight jobs each**, under **one shared cap of
three allocations** (reduced if the current production review reserve is
smaller). Each uses independent ROHF/AVAS CASCI at 2, 3.5, 6 and 30 angstrom
in aDZ-F12 and aTZ-F12. There are no CCSD(T), CASPT2, geometry optimization,
or potential-validation calculations in this batch. Printed roots and
natural-orbital AO coefficients/occupations support subsequent human selection.

| Screen | Source | Diagnostic space | Averaged doublet/quartet roots |
|---|---|---|---|
| C + SH3, C–S cleavage | Actual new review for `471430450000000000002-s06a7b4c221beb748` | CASCI(13,11), C/S valence plus all H 1s targets | 3 / 3 |
| CH2 + HSS, C–S cleavage | Provisional search at `791671580700000000002-s36f1e2b924862f12` | CASCI(19,12), C 2s/2p and both S 3s/3p | 2 / 2 |
| CH2SH + S, S–S cleavage | Same provisional parent, separate bond/record | CASCI(19,12), same target classes, independently generated orbitals | 3 / 3 |

These larger **diagnostic spaces are not production active-space choices**.
The CH3SS spaces retain all valence electrons in heavy-atom AVAS targets;
they do not include independent H-antibonding targets. If the printed
character indicates those are needed, expand the diagnostic space before
selecting production orbitals. The roots are diagnostic averages, not an
instruction to count that many reactive surfaces. C(3P) and S(3P) motivate
three spatial components per coupled spin; CH2 + HSS uses one extra
diagnostic root per spin beyond its nondegenerate ground-fragment channel.

Preparation reads current requests and well inputs and writes only the
separate recon directory. A provisional channel must already be an active
or accepted `hom_sci` search in `kinbot_monitor.out`; rejected searches are
not staged. Known CH2SH orbital work is not in this plan. Channels already
supplied in the case input are skipped, including a published refresh that
has not yet removed its old review marker. Repeating preparation refuses an
existing campaign; repeating `run` resumes its manifests and validated
outputs. Manifest locks prevent a second controller submitting those jobs.
Independent screens continue if one screen fails. Startup/allocation retries,
node exclusions, exclusive jobs, local scratch and partition-maximum walltime
come from the existing recon execution path.

Run from the cluster checkout (no production-controller restart):

```bash
cd "$HOME/KinBot"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
KB_BATCH="$HOME/kinbot-runs/r7_additional_orbital_recon_20261021"
if [ ! -f "$KB_BATCH/review_launch.json" ]; then
  .venv/bin/python examples/anl/tools/requested_orbital_recon.py batch \
    "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7" \
    examples/anl/orbital_review_batches/r7_additional_20261021/batch.json "$KB_BATCH"
fi
.venv/bin/python -u examples/anl/tools/requested_orbital_recon.py run "$KB_BATCH"
```

Use the same helper's `monitor "$KB_BATCH"` and `evidence "$KB_BATCH"`.
When a provisional channel gets its actual review request, bind its completed
evidence before selecting orbitals:

```bash
.venv/bin/python examples/anl/tools/requested_orbital_recon.py bind \
  "$KB_BATCH/ch2shs_cs" \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/ch3ss/791671580700000000002-s36f1e2b924862f12/vrctst/orbital_review/791671580700000000002-s36f1e2b924862f12_hom_sci_1_2/request.json"
```

Binding checks parent, atom ordering/connectivity, bond, charge, spin and
separated fragment identities. A change in optimized coordinates is recorded
and both geometries retained. Binding does not approve an active space or
change production. A reviewed selection recipe still goes through `review
--apply`, which uses the existing live refresh. S–S binding uses the
`ch2shs_ss` child and `hom_sci_2_3` request. Native evidence remains separate
if a tentative reaction is later rejected.

## Completed additional screens and selected production spaces

All 24 orbital-only probes completed on October 21, as confirmed by the
subsequently imported native archive described below.
The last aTZ S-S far record confirms the same orbital character as aDZ.
The following selections use **each screen's own 2 Å natural record 2140.2**
in each basis. Orbital indices are source-file positions, never a transferable
label across different molecules or dissociation channels.

| Channel | Selected CAS and source indices | Frontier content | Orbital average | Reactive leaves |
|---|---|---|---|---|
| C + SH3 | (3,4), 12–15 | C three 2p orbitals / two electrons, SH3 radical / one electron | 3 doublets + 3 quartets | D1 |
| CH2 + HSS | (5,4), 19–22 | CH2 two frontier orbitals / two electrons, HSS frontier pair / three electrons | 2 doublets + 2 quartets | D1 |
| CH2SH + S | (5,4), 19–22 | CH2SH radical / one electron, S three 3p orbitals / four electrons | 3 doublets + 3 quartets | D1 |

The explicitly reviewed recipes are
`examples/anl/orbital_review_batches/r7_additional_20261021/*_selected.json`;
`selections.json` groups them by production case. Each root has equal weight
within its spin block (1/3 or 1/2); the two blocks have equal total weight.
These orbital-optimization weights are not electronic statistical weights.
Molpro's [MULTI documentation](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi)
describes the WF/STATE/WEIGHT declarations. D1 is the chosen ground-parent
reactive path; additional averaged roots support the orbitals and spin-gap
calculations. Degeneracy of an atomic fragment does not by itself prescribe
additional reactive parent paths. These new selections keep explicit Q1 for
the MR reference, as did the preceding CH2SH selection. The original three
channels retain their separately reviewed quartet selectors.

The evidence supplied for these selections is:

- **C + SH3:** far orbital 11 is carbon 2s, occupation 1.9559; it is closed.
  Orbital 12 has the SH3 S/H coefficients and occupation 1.0030. Orbitals
  13–15 span three carbon p directions with occupations about 0.6814 each.
  Near occupations of 12–15 are 1.2420, 0.6900, 0.6855, 0.4532. Both bases
  agree. The fractional sum in this larger correlated space is not an integer
  electron-count prescription; the compact CAS contains the three frontier
  electrons. Retaining C 2s would instead be a five-electron/five-orbital
  expansion, not the earlier CH3S CAS(5,4).
- **CH2 + HSS:** far 19/20 are on both sulfurs, about 1.5095/1.5050, while
  21/22 are carbon-centered with occupation 1 each. Near values are 1.7403,
  1.5040, 1.1602, 0.6081, consistently across bases. The 18/23 pair is largely
  filled/empty (far 1.958/0.049). The two diagnostic asymptotic roots differ
  by about 0.017 Eh: the averaged 1.5 occupations do **not** prove a doubly
  degenerate HSS ground state. The second root remains explicit diagnostic
  orbital support, not a second counted reaction.
- **CH2SH + S:** far 19–21 are entirely the departing S p directions with
  occupations 1.33333 each; 22 has carbon/bound-S radical character and
  occupation 1.0062. The last aTZ record agrees. Near aTZ occupations are
  1.57253, 1.53627, 1.36968, 0.60832. The 18/23 spectator pair is about
  1.977/0.027 far. This is a different fragment frontier space from the
  initial CH3SS -> CH3S + S channel's CAS(7,5).

This is orbital assignment from the larger-space CASCI evidence, **not** a
completed optimized compact-CAS, angular-potential or rate-accuracy
validation. Normal production performs the compact CASSCF/CASPT2 calculations.
Do not resume an extra CC or correlated validation campaign for this handoff.

### Geometry provenance and the live refresh

The previous CH2SH review used the actual accepted request's atom ordering
and parent geometry. Its live-refresh receipt was reported **applied** for
all three participating workers. The C + SH3 request was rechecked against
the new recon snapshot: identical parent, bond, symbols, charge, spin and
products, with zero change in pair distances. These requests contain the
accepted L2 parent geometry: B2PLYP-D3(BJ)/cc-pVTZ for r7. Recon rigidly moves
one connected fragment to 2, 3.5, 6 and 30 Å, preserving intrafragment
coordinates. It does not claim these are independently relaxed MR points.

The CH3SS screens used the saved well-input geometry before its two accepted
review requests existed. The 17:49 monitor and subsequent filesystem check
still had no such requests. **Keep those selections prepared.** The helper
will bind them only when actual requests appear and their parent, atom
ordering/connectivity, bond, fragment identities, charges and spins agree.
It saves both geometries and reports the maximum pair-distance change.
A changed molecular identity is refused. A changed geometry is evidence to
inspect, not proof that the electronic state must remain the same.

The original CH3S and two CH3SS seeds retain their own historical near/far,
angular and expansion evidence and separate DZ/TZ WFU inventories. They were
not generated from every current r7 coordinate. This is legitimate seed reuse:
production opens FILE,2 **before** declaring the current geometry, basis and
spin, so the WFU's stored geometry cannot replace the requested structure.
The shared renderer, ASE calculator and sampling regression tests cover that
ordering. The current production geometry is optimized with the reviewed MR
method; a recon geometry is never installed as a converged production result.

Run this after pulling the repair; it starts no driver or recon calculation:

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
.venv/bin/python examples/anl/tools/requested_orbital_recon.py review-batch \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7" \
  "$HOME/kinbot-runs/r7_additional_orbital_recon_20261021" \
  examples/anl/orbital_review_batches/r7_additional_20261021/selections.json \
  --apply
BASH
```

At the reported state, this adds C + SH3 to CH3S and reports CH3SS as
`waiting_for_accepted_channel`. Rerun the same command when those reviews
appear. It reuses all saved recon evidence, groups ready additions into one
revision per case, and reconciles repeated submissions. If an older refresh
is pending it reports `waiting_for_previous_refresh` without overwriting that
input or receipt. No existing model is replaced, no QC is cancelled, and the
input's current grids, basis policies, sampling accuracy and resource limits
remain intact. It uses the existing pinned driver's live-refresh protocol;
loading this command does not require replacing the running driver.

The refresh repair adds an atomic proposal API: the publisher validates the
expected input, owner, workers and prior revision while holding its existing
control/phase locks, **before** writing an edited input. A separate helper
lock serializes channel additions. Tests cover stale writers, pending
revisions, identical resubmissions, joint channel additions, provisional
binding and inheritance by existing workers without duplicate driver startup.

### Evidence retained during production

Scan/asymptote steps retain executed input, native output, actual geometry,
optimizer history, basis-specific WFU/checkpoint provenance and
`*.orbital_audit.json`. Generated MR correction and sampling decks also
request native orbital/CI diagnostics. Reports include OCC/CLOSED, spins,
roots, weights, selected energies, natural occupations, AO coefficients,
CI/reference evidence, geometry/basis records and input/output hashes.
See [native diagnostic details](MOLPRO_NATIVE_ORBITAL_DIAGNOSTICS.md).

Printed AO/CI coefficients have the configurable default magnitude cutoff
0.05. They are not a complete overlap matrix or automatic proof against root
flipping. Quantitative state/subspace continuity can need the retained WFU
and additional analysis. Solver/geometry convergence plus a successful audit
parser alone does not establish the scientifically intended electronic state.
Keep the recon directories: their hashed WFU files are runtime seeds.

## Evidence archive coverage and export

The documentation records selections, rationale, geometry provenance and the
native excerpts supplied in the conversation. It is **not a substitute for the
full native outputs**. The export captured on October 21 at 18:19:41 PDT was
subsequently imported into
`examples/anl/tools/kinbot-r7-orbitals-20261022T011938Z-2923454/`.
Its 197,101,264-byte `orbital-evidence.tar.gz` has SHA256
`0396956230d8dabe0b216bd9cf90aa12cce88dcfdbdf86821a26d2c3c5eeff0f`.
The archive checksum and **all 2,951 member sizes and checksums were verified**;
the internal and sidecar manifests match exactly, with no collection issues.
The large archive stays in that local directory, outside Git. The tracked
[native export inventory](../examples/anl/orbital_review_batches/r7_additional_20261021/native_export_inventory.json)
records its identity, campaign statuses, native-file hashes, installed seed
bindings and refresh acknowledgements.

| Evidence | Source below `$HOME/kinbot-runs` | Coverage verified in the archive |
|---|---|---|
| CH2SH -> CH2 + SH orbital screen | `r7_ch2sh_cs_orbital_recon` | Eight completed aDZ/aTZ orbital probes |
| Optional compact/expanded CH2SH panel | `r7_ch2sh_cs_orbital_validation` | Four executed probes; the other 18 were not submitted. The compact near-point native run completed but its optional CC energy reader rejected the output. Do not resume this panel merely to export it. |
| C + SH3 and the two prospective CH3SS cleavages | `r7_additional_orbital_recon_20261021` | All 24 orbital-only probes complete; prospective channels still require accepted-channel binding |

The new r7 campaign trees contain **36 native outputs, 36 WFU files and 108
XML files**: 32 orbital-only screens plus four executed optional validation
probes. They contain no Molden files. The 18 unsubmitted validation inputs
are retained as inputs, not counted as completed calculations. For the 35
jobs with recorded provenance, prepared-input, output, WFU and both orbital
XML hashes match the archived bytes. All ten WFU seeds referenced by the
installed production models also match their recorded hashes. The optional
compact probe's native output and failed validation record are both retained.
The orbital-screen manifests still say `coverage_refinement_required`;
preserving their evidence does not turn them into a completed potential
validation or waive the provisional CH3SS review requirement.

The original CH3S C-S and CH3SS C-S/S-S selections already have curated native
archives, individual WFU seeds, hash inventories and review records under
`examples/anl/{ch3s,ch3ss}_anl1_f12/orbital_reviews`. The two
`caspt2_native.evidence.tar.gz` archives contain the accepted native evidence;
their `native_inventory.json` and `orbital_inventory.json` retain original
locations and hashes. The original WFU sources were in
`minimum_active_space_molpro_confirmation_20261016_r2` (CH3SS C-S) and
`minimum_active_space_molpro_confirmation_20261016_r2_recovery` (CH3S C-S and
CH3SS S-S). The export uses the curated copies, so those older source
directories need not still exist. Any historical archives present under the
two review trees are included as well.

The standalone exporter uses only the Python standard library. It includes
**all regular files** in the three r7 campaign directories, including prepared
and executed decks, full native outputs, XML, Molden, WFU, orbital audits,
manifests, failed-attempt evidence and launch records. It also captures the
curated original evidence, selection recipes, documentation, production inputs,
saved review requests and live refresh receipts. It does not run an audit that
rewrites manifests, query Slurm, submit QC, acquire controller locks, or change
source files. Production QC output trees and node-local scratch are outside
this orbital-recon export's scope. Keep them separately for the final rate
calculation audit.

Run on Blodgett after pulling `composite`:

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
KB_EXPORT="/tmp/kinbot-r7-orbitals-$(date -u +%Y%m%dT%H%M%SZ)-$$"
nice -n 10 .venv/bin/python examples/anl/tools/export_orbital_evidence.py \
  examples/anl/orbital_review_batches/r7_additional_20261021/export.json \
  --repo "$HOME/KinBot" --runs "$HOME/kinbot-runs" \
  --production "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7" \
  --output-directory "$KB_EXPORT"
printf 'COPY THIS DIRECTORY: %s\n' "$KB_EXPORT"
BASH
```

This can be a multi-GB export because it preserves the wavefunctions. It
streams each source once, compresses at level 1 and then verifies all archive
members on the export disk. The output includes `orbital-evidence.tar.gz`, its
`.sha256` and a readable `.manifest.json` with every source/member/hash/size.
No existing export is overwritten. Missing required sources, skipped symlinks
and files changed during their read are explicitly reported; exit status 2
means inspect that manifest, not that a calculation failed. An interrupted
archive remains `.partial`. Live metadata is a set of per-file snapshots,
not an atomic snapshot of the whole running workflow. Transport verification
does not mean scientific validation or that unsubmitted jobs ran.

Copy the printed directory from a local terminal, then verify its checksum:

```bash
scp -r hram@blodgett.ca.sandia.gov:/tmp/kinbot-r7-orbitals-TIMESTAMP-PID ./
cd kinbot-r7-orbitals-TIMESTAMP-PID
shasum -a 256 -c orbital-evidence.tar.gz.sha256
```

### Observed refresh wait, October 21 at 18:07 PDT

The C + SH3 revision `38d137b545fa4efba054ebed93fa360f` was published as three
worker requests but remained unacknowledged in the monitor. The prior CH2SH
review was already applied: its first production point completed at 18:02:31.
That worker then yielded; the original CH3S worker remained in VRC with
retained job 52915984. Both new CH3SS selections still awaited actual accepted
requests. These are separate from the 750-step exploration failures displayed
in the same monitor.

The cause of this delay is a missing cooperative epoch check inside the
dispatcher polling loop. The outer VRC loop and new submissions checked for
refreshes, but an already submitted point could be polled until it finished
before the controller surrendered `drive.lock`. PES deliberately grants no
revision while any branch still owns its graph. The repair checks the epoch
before every dispatcher poll, including scheduler-outage and failed-task hold
loops. It releases only controller graph ownership; submitted job IDs and
native files remain intact for reattachment. Regression tests verify the
unchanged ledger and released lock. Standalone dispatch has no PES epoch and
continues as before. The read-only monitor now names the branch controllers
still awaited when `quiescing` is true, rather than only showing a queued count.

An already running pinned Python controller does not load this code merely
because the checkout is pulled. Without replacing it, the current queued
revision can proceed when the remaining point reaches its existing boundary.
No scan cancellation or recon resubmission is required to let that happen.

The transferred archive contains all six worker acknowledgements for the
CH2SH and C + SH3 revisions, each `applied` without an error. The 18:22:42 PDT
monitor confirms this and shows C + SH3's asymptotic geometry job 52916238
running, alongside CH2SH point 01 and original CH3S point 07. Thus this
observed refresh wait has cleared; a restart is not needed for these models.

### Serial scheduling observed before the parallel-channel update

The original CH3SS S–S channel already has its reviewed model. In the
18:22:42 PDT monitor it has not submitted a scan because the current
implementation processes channels **serially within each well's worker**.
`VTS.calculate_correction_potentials` calls `calculate_channel` synchronously
for each selected reaction; that function completes each electronic leaf's
scan and correction energies before returning. Consequently CH3SS C–S does
its asymptotic geometry, pivots, relaxed points and correction assembly before
S–S begins. CH3S D2 similarly follows D1. Other wells' VRC workers and PES
exploration can run concurrently under the shared resource cap.

This was an implementation limitation, not a scientific dependency between
the two bonds or a missing S–S active-space review. All 23 shared allocations
were occupied in this snapshot, but S–S was not yet a Slurm pending job.
The following update removes that serialization for shared-pool runs.

### Parallel channels and electronic surfaces

Reviewed MR channels now advance through a bounded pool of in-process surface
controllers. The default `vrc_tst_parallel_surfaces: 4` permits both initial
CH3SS bonds and both CH3S doublet surfaces to progress independently. This
controls local surface orchestration, **not additional Slurm allocations**.
The existing shared resource policy continues to enforce the common limit,
node exclusions and orbital-review reserve. Without a shared resource policy,
execution remains serial so independent graphs cannot multiply a legacy
per-graph node cap. Setting the input value to 1 also selects serial execution.

Each surface still runs its MR asymptotic geometry, fixed-geometry pivots,
ordered relaxed points and correction energies in dependency order. Every
point uses its own preceding geometry/WFU checkpoint. Review-held channels
are omitted until their model is accepted; approved channels can proceed.
A native failure is recorded per surface while independent surfaces continue.
The main controller assembles a channel's correction manifest only after all
its required surfaces succeed. ROTD generation/execution and plotting remain
on the main thread; this update does not change Monte Carlo or MESS ordering.

Graph identities, surface names, geometry inputs, methods, bases, roots and
statistical weights are unchanged. The in-place refresh adopts completed and
allocated graphs through their existing bindings. Identical pivot requests
share one graph controller; graph creation, database publication and correction
batch metadata have serialized access. A partially copied-back correction
batch now reattaches its remaining tasks to the matching existing graph,
including its Slurm IDs, instead of preparing a duplicate subset batch.
Live PES refresh waits for **every** surface controller to release its graph
lock before a revised model can be applied; it does not cancel QC allocations.

`vrctst/parallel_surfaces.json` records queued/running/complete/failed/paused
surface-controller states. Here `running` means the controller is advancing
the surface; its QC task may still be waiting for shared-pool capacity.
`run_progress` includes all current correction batches without double-counting
a batch reused by two surfaces. The ordinary monitor continues to report the
actual Slurm tasks and their states.

Install this code with one full **controller-only** refresh of r7:

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python -u -m kinbot.run_refresh \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7" --apply
BASH
```

This command uses the same run directories and preserves running/pending QC,
completed results, review models and the resource policy. It does not request
pending-job replacement or a scientific migration. The new source snapshot
applies to replacement PES/KB controllers; previously allocated QC keeps its
original verified worker code. After `phase: launched`, retained Slurm IDs
should still appear and the initial CH3SS S–S and CH3S D2 asymptotes become
eligible alongside the existing surfaces. A full shared pool can still delay
their allocations until other jobs finish.

Regression coverage includes simultaneous bonds/roots, independent coordinates,
per-surface failures, main-thread plotting, shared-pivot reuse, concurrent graph
creation, partially completed correction-batch adoption, unchanged allocated
ledgers on epoch refresh, and joining all controllers if thread startup fails.
Native electronic-structure inputs are unchanged; these are orchestration and
restart tests, not a new Blodgett native calculation.

## D1/D2 counting at the MESS boundary

The current single-bottleneck handoff already implements Union's operation:

`N(E) = sum_i weight_i * g_i * N_i(E - offset_i)`.

Each surface contributes once, with `g_i=2` for an unsplit doublet. Explicit
fine-structure levels, if supplied for a surface, replace that default spin
factor. State-averaging weights do not enter. Conserved fragment vibrations
and rotor factors are shared by these surfaces and convolved once after
the sum; MESS receives `ElectronicLevels 0 1` for this combined barrier.
The CH3S well's separate electronic-level partition function is retained.

MESS `Union` likewise sums constituent microcanonical state counts and
Boltzmann-shifted canonical weights: see the
[official manual, section 6.3](https://tcg.cse.anl.gov/papr/codes/mess/messmanual_3_23_16.pdf)
and `UnionSpecies::states` / `UnionSpecies::weight` in the
[official source](https://github.com/Auto-Mech/MESS/blob/master/src/libmess/model.cc).
An additional Union of the two leaves **and** the already summed table
would double count. Replacing this valid sum with a second aggregation is
not required. Tests integrate unequal-threshold doublet fluxes at four
temperatures and verify the analytic Union weights; identical doublets
contribute four spin-resolved states, not eight and not their average.
The manifest reader now also rejects duplicate spin/root entries or a
repeated scalar input within one bottleneck.

This equivalence assumes the configured common conserved-mode model and
the shared sampled energy range. Distinct state-specific fragment vibrational
models would require separate convolutions before summation. It does not
add distance-dependent spin-orbit or nonadiabatic coupling; those remain
limitations of the selected spin-free adiabatic model, not problems solved
by writing `Union`. Existing stereochemical barrier Unions remain separate:
they combine distinct spatial reaction pathways after their own duplicate
classification.

## Ground terms, surfaces and averaging are different quantities

CH3S has a ground electronic term X²E. At C3v and without spin–orbit
coupling, E has two degenerate orbital components; the superscript 2 denotes
spin multiplicity. This is an electronic-structure property, not a selectable
orbital count. Symmetry-lowering Jahn–Teller distortions can split the orbital
components. Spin–orbit coupling changes the physical levels and competes
with the distortion. Two exactly equal energies everywhere along a general
scan are not implied. See [Marenich and Boggs, 2004](https://doi.org/10.1021/jp040443f).

Treating those components as two independent, uncoupled spin-free adiabatic
VRC surfaces is a **model choice**. The production input selects doublet roots
1 and 2 for CH3S. It does not implement coupled spin-vibronic dynamics or a
distance-dependent spin–orbit Hamiltonian. The chosen roots must retain the
intended character; root numbers alone do not prove that they do.

CH3SS's ground-parent term is X²A″ in Cs, which is orbitally nondegenerate;
see the [photoelectron assignment](https://doi.org/10.1063/1.1343081).
Consequently its two chemical dissociations each select one ground-parent
doublet surface. The S–S products nevertheless include CH3S(²E) and S(³P),
so that channel needs a much larger *state average*. Product-manifold
degeneracy does not make every asymptotic component a separate reaction into
the ground-parent well. Selecting only the ground-parent adiabatic pathway
is not proof that excited-surface or nonadiabatic contributions are zero.

| Chemical channel | Reactive scalar relaxed paths | Active space | Averaged roots (doublet / quartet) | Relative input weight per root (doublet / quartet) |
|---|---|---|---|---|
| CH3S -> CH3 + S, C–S | D1 and D2 | CAS(5,4) | 3 / 3 | 1/3 / 1/3 |
| CH3SS -> CH3 + S2, C–S | D1 | CAS(7,5) | 1 / 1 | 1 / 1 |
| CH3SS -> CH3S + S, S–S | D1 | CAS(7,5) | 6 / 6 | 1/6 / 1/6 |
| CH2SH -> CH2 + SH, C–S, new review | D1 | CAS(5,4) | 4 / 2 | 1/4 / 1/2 |

D1/D2 mean roots within spin multiplicity 2, not the two spin projections of
one state. Every reactive leaf has spatial weight 1 and a common spin-free
ground-fragment asymptotic reference. Averaging weights normalize to a total
of 1/2 per spin block: each CH3S root has total weight 1/6, each CH3SS C–S root
1/2, each CH3SS S–S root 1/12, and each new CH2SH doublet/quartet 1/8 and 1/4.

State averaging optimizes one compromise orbital set across the specified
states. It reduces bias toward a single component and helps avoid orbital
optimization root flipping; it is not a guarantee against crossings or an
occupation-based proof of state identity. [Molpro's state-averaging guidance](https://www.molpro.com/manual/doku.php?id=quickstart#state-averaged_mcscf)
describes this distinction. These weights are not thermal populations,
spin-statistical weights or an instruction to sum all averaged roots as flux.
The VRC assembly separately combines the declared reactive leaves, includes
their spin factors, and writes its combined MESS barrier with electronic
degeneracy 1 to avoid multiplying that factor again.

The occupations below are eigenvalues of the state-averaged one-particle
density, not populations of individual CASPT2 states. A value such as 1.5
can arise from averaging complementary occupations. AO coefficients locate
an orbital on particular atoms and show its angular/bonding character, but
their squares are not population percentages in a nonorthogonal AO basis.
The overall orbital sign is arbitrary; relative phases and the complete
subspace matter. Occupations alone are therefore insufficient for selection.
Each index range below refers to its own reviewed WFU file and record 2140.2,
not to a universal ordering shared by the molecules or channels.

## CH3S C–S: five electrons in four orbitals

CH3 contributes its singly occupied carbon frontier orbital (one electron),
and atomic sulfur contributes its complete 3p manifold (four electrons in
three orbitals). All three sulfur p directions are retained because angular
sampling must not privilege one atomic component. S(³P) has three spatial
components; coupling its triplet spin with doublet methyl yields three
doublets and three quartets in the spin-free ground-fragment manifold.
The parent X²E components motivate two reactive doublet paths; the third
averaged doublet is not automatically another ground-parent flux path.

The reviewed DZ/TZ natural records use indices **11–14**. Larger-parent
coefficients place 11–13 on sulfur p and 14 on methyl carbon. Their confirmed
far occupations are `(1.33333, 1.33333, 1.33333, 1.00000)` in both bases. At
the inner TZ point they become approximately
`(1.49748, 1.49748, 1.47637, 0.52867)`, showing bonding redistribution inside
the retained space. C–H bonds and sulfur 3s are inactive.

Adding sulfur 3s in CAS(7,5) changed the tested doublet angular residual by
0.00039 kcal/mol and the doublet–quartet angular-gap residual by 0.01640.
The occupied-valence (13,8) comparison gave 0.10803 and 0.03367. These support
the minimum choice on the tested cuts, not a universal accuracy bound.
[CH3S evidence and native inventory](../examples/anl/ch3s_anl1_f12/orbital_reviews/ACTIVE_SPACE.md).

## CH3SS C–S: seven electrons in five orbitals

The fragments are methyl doublet and ground triplet S2. Keep methyl's
one-electron frontier orbital plus the complete S2 pi/pi-star manifold:
six electrons in four orbitals. S2's intact sigma bond is outside this
minimum CAS. Unlike atomic S(³P), the ground S2 term is an orbitally
nondegenerate Sigma term; there is one spatial ground-fragment state for
each of the doublet and quartet spin couplings, hence one root per spin.

The larger parent record requires original indices **17,18,20,21,22**.
Preparation swapped 17 and 19 so the saved production records contain the
chosen set contiguously at **18–22**. Do not apply that swap again. Bonding pi
coefficients have corresponding same-phase sulfur contributions; the pi-star
partners have opposite phases. The methyl orbital is carbon-dominated.
Confirmed far TZ occupations are
`(1.95828, 1.95828, 1.04172, 1.04172, 1.00000)`.

The omitted sigma/sigma-star pair is not exactly uncorrelated: its larger
parent occupations are about 1.94564/0.06127. The explicit CAS(9,7) expansion
changed the tested doublet angular residual by 0.00239 kcal/mol and the
angular spin-gap residual by 0.00095, supporting its exclusion for this
channel. [CH3SS evidence](../examples/anl/ch3ss_anl1_f12/orbital_reviews/ACTIVE_SPACE.md).

## CH3SS S–S: a different seven-electron, five-orbital space

Keep the CH3S sulfur-centered frontier pair (three electrons) and all three
departing atomic sulfur p orbitals (four electrons). The two CH3S orbitals
support complementary (2,1)/(1,2) occupations. Freezing one as a permanent
lone pair would remove the other hole configuration. This is the reason for
CAS(7,5), rather than atomic S plus only one CH3S radical orbital.

The two CH3S spatial components times the three atomic S components give
**six spatial roots per coupled spin**. That motivates six doublets and six
quartets in the average. They need not be all degenerate at a distorted
fragment geometry. Only the selected ground-parent D1 path contributes
reactive flux in this model.

The channel-specific natural seeds also use indices **18–22**, but their
content differs from the C–S channel. Indices 18/19 belong mainly to bound
sulfur, and 20–22 to departing sulfur. Both bases give far occupations
`(1.50000, 1.50000, 1.33333, 1.33333, 1.33333)`. At the inner TZ point they
become `(1.71816, 1.63491, 1.49536, 1.36269, 0.78888)`. Canonical orbitals can
reorder these fragment groups inside the same CAS without losing the span.

Adding the intact C–S sigma/sigma-star pair in CAS(9,7) changed tested angular
doublet/gap residuals by 0.04494/0.05703 kcal/mol. At 2.1 Å, corresponding
interaction/gap differences reached 0.114/0.316 kcal/mol. There is no claim
of uniformly sub-0.1 kcal/mol agreement. The tested active-span overlaps and
return calculations support continuity, not universal tracking at every
future geometry. [Detailed native evidence](../examples/anl/ch3ss_anl1_f12/orbital_reviews/ACTIVE_SPACE.md).

## New CH2SH C–S review: five electrons in four orbitals

Triplet CH2 has two frontier electrons in two orbitals; SH(²Pi) contributes
three electrons in its pi pair. This gives CAS(5,4), source indices **11–14**
in each basis's original CASCI record. SH contributes two spatial ground
components, yielding two ground-fragment doublets and two quartets. The
four-doublet/two-quartet average deliberately includes two extra diagnostic
doublets; their energies/occupations alone do not identify their asymptotes.
The selected reaction into the nondegenerate parent uses D1 only. The new
recipe retains explicit Q1 for the MR reference; the earlier three channels
use their reviewed lowest-correlated-quartet selector because their native
quartet root ordering was observed to change.

At 2 Å, the completed aDZ CASSCF(5,4) occupations are
`(1.70808, 1.48606, 1.02341, 0.78245)`, with the intended SH/CH2 AO character.
The CAS(7,6) checks at 2, 3.5 and 30 Å leave the added S–H sigma pair near
1.98/0.02 and preserve the frontier pattern. Both original CASCI bases give
consistent near/far character. This supports orbital selection without
completing the optional 22-job potential panel. Full two-basis optimized-CAS
and angular comparisons are not claimed. [Evidence, recipe and refresh command](CH2SH_ORBITAL_REVIEW.md).

## Shared production sequence and controls

Each reviewed reactive leaf uses its own chosen CAS, state average and
basis-specific hashed natural-orbital seed. The source inputs prescribe:

1. A constrained separated-supermolecule optimization at the 30 Å reference,
   using CASPT2-F12/aDZ-F12; extract the fragment geometries. This reuses the
   full channel CAS and requires no separate fragment active spaces.
2. Fixed-geometry L2 Gaussian B2PLYP-D3BJ/cc-pVTZ fragment orbital jobs for
   HOMO pivot directions. These do not replace MR fragment geometries or
   enter the MR interaction energy. Atomic fragments need no such job.
3. CASPT2-F12/aDZ-F12 relaxed paths and matched frozen-fragment calculations,
   retaining the existing distance constraints, without imposing new C3v/Cs
   geometric symmetry. The 47 distances are 2–8 Å by 0.2, 8.5–12 by 0.5,
   and 13–20 by 1. The separate asymptote remains 30 Å.
4. High-level reference curves from ROHF-UCCSD(T)-F12b on the quartet plus
   CASPT2-F12 doublet-minus-quartet gaps. Both the CC reference and MR gaps
   use aDZ/aTZ-F12 extrapolation with power 3.7. In Molpro, `rhf` is the
   restricted-open-shell HF reference for these open-shell calculations.
   For one gap method the construction is
   `E_target = E_CC,quartet,CBS + E_MR,target,CBS - E_MR,quartet,CBS`.
5. Assemble the relaxation/high-level correction and perform CASPT2-F12/aDZ
   orientation sampling with RotdPy at the same 47 dividing distances.
   The relative flux sampling target is 5%; configured sample bounds still
   apply. Several pivot alternatives may exist at a distance, so 47 is not
   necessarily the number of RotdPy Surface objects.

CASPT2 uses IPEA 0.25 Eh, real shift 0.1 Eh and the reviewed same-spin
state-average reference policy. The earlier three channels have separate
DZ/TZ WFU files already in their case directories. The new recipe points to
the completed reconnaissance's two hashed seeds. The source configuration
therefore had four reviewed chemical channels and five reactive leaves after
the CH2SH addition. C + SH3 adds one channel/leaf; the two provisional CH3SS
selections can each add another after their actual requests appear. This does
not describe every channel future exploration may find.

The live command preserves the actual input's grids, levels, numerical and
resource controls rather than reinstalling defaults. Existing scratch-first,
exclusive-node, partition-walltime and shared-admission policies are untouched.
Other reviewed paths continue; a new unspecified channel remains held for
review. ANL still waits for outer-PES L2 selection. Neither successful orbital
review nor a passing parser certifies eventual convergence of every native
point or the complete kinetics calculation.

### Review status during production

Older controllers retain `vrc_model_review.json: needs_review` until the
entire VRC branch finishes, even after individual
`vrctst/orbital_review/<reaction>/request.json` records become
`model_accepted` and scans start. The read-only monitor now checks those
channel receipts: newer, matching acceptances supersede the old aggregate
hold. Missing, unreadable, older or unaccepted receipts remain visible as
outstanding or unconfirmed reviews. The original aggregate status/error is
retained in `--json` and `--details`. Model acceptance does not certify QC
completion. This display correction takes effect on the next monitor
invocation after pulling; existing drivers and jobs need no refresh.

The overview labels each branch pair by its well and, when available, shows
that well's current PES stage. `not_reported_yet` means a retained branch
report predates the replacement controller; an `exploring` worker can
legitimately keep that label until it reaches post-reaction dispatch.
