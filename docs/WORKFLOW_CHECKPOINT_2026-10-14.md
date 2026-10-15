# Workflow checkpoint and pause — 14 October 2026

**Resumed by the user.** This is the historical stopping record. Development
has resumed; the current recovery and production-default changes are recorded
in `Blodgett_workflow_runbook.md`. The existing eight-basis campaign remains
unchanged by those updates. The subsequent conventional-5Z RHF completed,
using about 138 GB of scratch, and failed only at the unsupported Molden
export. The native XML/WFU recovery preserves that converged SCF evidence.

Work pauses at the user's request after this checkpoint. The user runs all
SSH. Existing and subsequently resumed Blodgett calculations may continue
while development is paused. This document supersedes the earlier
`WORKFLOW_PAUSE_2026-10-14.md` snapshot. It is a stopping record, not a claim
that every desired workflow feature or chemistry validation is complete.

## Published implementation and native evidence

The `composite` branch contains these checkpoints:

- `4793a0f`: general orbital campaign, state/term models, native-basis orbital
  probes, multistate scalar ROTD orchestration, MESS state sums, CBH/ANL
  corrections, scoped cleanup, tests and initial audit.
- `18abd0d`: explicit `/home/hram/KinBot` paths; removed the unsupported
  unguarded Slurm-accounting query from Blodgett instructions.
- `435f012`: corrected an extra ROHF closing brace; validates generated
  delimiters; sets each job's WFU repository with `-W`; first native probe
  gates wider submission; native execution/error records. Production example
  selections now restrict homolysis to the requested input atom pairs, retain
  other reaction families, and require reviewed surface models.
- `a8e1a21`: accepts both Molpro `HARTREE` and legacy `AU` energy labels in
  reconnaissance and ASE/sampling readers. ANL F12 validation recognizes both
  and accepts consistent assignment/SHOW duplicates while rejecting mismatches.
- `3f09546`: clears stale campaign errors when a controller resumes.
- The final checkpoint adds an archived, scheduler-checked failed-job retry
  and a guard against silently dropping required ROTD data in full-PES MESS
  assembly. Its hash is supplied in the handoff.

The replacement CH3S augmented-DZ native pilot succeeded with return code zero,
three doublet and three quartet roots, normal termination, both Molden files,
and a nonempty job-local WFU. The original replacement audit rejection was a
reader error (`HARTREE` versus `AU`), not a failed electronic calculation.
The most recent user-provided status reports **four validated CH3S pilots**,
conventional 5Z failing during integral scratch I/O, and three additional
pilots running at 19:44. The user subsequently reports no jobs remain active.
Re-audit the existing outputs rather than assuming which other jobs succeeded.
CH3SS pilots and all full scans had not started in that status snapshot.

The current campaign is:

```text
/home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525
```

`/home/hram/kinbot-runs/last_sulfur_orbital_campaign.txt` is the authoritative
pointer for commands, so the workflow does not depend on this timestamp.

## Scratch failure and narrow recovery

Molpro reported a failed FILE21 write at offset 569671680 eight-byte words
(about 4.56 GB), ending in `writew`/I/O error and exit 255. The message alone
cannot distinguish disk exhaustion, quota, file-size limits, or another I/O
failure. It is not a basis-availability or energy-parser error.

The compute-node inventory supplied by the user shows `/tmp` has about 37 GB
free, `/scratch2` about 172 GB, and `/scratch` about 5.4 TB. `/scratch` is a
compute-node mount, absent from the login-node inventory. The existing runtime
selector only requires a 4096 MB free-space floor and can select `/tmp`; it
does not forecast a conventional 5Z integral build's disk demand. The selected
scratch and initial free capacity are recorded per native execution.

Use the existing generic runtime input:

```bash
export KINBOT_MOLPRO_SCRATCH="/scratch/$USER/kinbot-molpro"
```

Set it on the login node before starting the controller. No login-node
`mkdir /scratch` is needed: the allocated compute-node worker creates and
checks its own directory. Submission explicitly exports the controller's
environment. Each native job creates an isolated subdirectory and cleans
only that subdirectory. This selects a filesystem; it changes no electronic
method, basis, auxiliary choice, or previously completed result.

`python -m kinbot.orbital_campaign retry-failed ROOT` requires the campaign
controller lock and successful live scheduler queries. It refuses a manifest
with live allocations or uncertain submissions. It re-audits outputs, preserves
completed probes, and archives known failed native attempts (including inputs,
outputs, WFU and execution records) under each job's `attempts/NUMBER`. Only
then does it clear that failed attempt's old scheduler ID and make it eligible
for resubmission. Original input hashes and predecessor bindings remain
unchanged. Parser-only/unknown failures are not automatically recalculated.

After selecting scratch and marking failed attempts retry-ready, resume the
same campaign. Do not archive the whole campaign or manually submit individual
`.sbatch` files. Successful pilots are reused. If another native failure
occurs, retain its diagnostics; no unattended unlimited retry is implemented.

## What is implemented for these examples

The examples use native augmented conventional D/T/Q/5Z and native augmented
F12 D/T/Q/5Z orbital bases. **Default Molpro F12 definitions**, no rev2 override,
no generated/vendored basis data or auxiliary substitution. They treat CH3S
C–S and CH3SS C–S/S–S as separate channels. Spins follow input species and
fragment models; spatial-root counts follow supplied terms rather than
molecule-name branches. Physical state data and atom pairs are example inputs.

The initial full grids comprise 248 CH3S and 600 CH3SS probes, plus 24 native
pilots: 872 jobs before adaptive additions. All system pilots precede full
scans. System controllers run sequentially and wait for allocations to leave,
so the campaign cap is ten one-node jobs total. Radial probes span near-bond
geometries through an initial 30 A, with rigid-fragment orientation probes.
Every averaged root, energy steps, occupations and asymptotic plateaus inform
refinement. Exhausted refinement is an explicit unresolved status.

Reports expose root-labelled energies, natural populations, occupation-based
orbital candidates, Molden/CI-character evidence and hash-bound WFU records.
They provide diagnostic active spaces, not an automatically certified minimal
CAS or proof of diabatic root continuity. `select` can create reviewed,
basis-specific orbital-record models after chemical review.

Single-well multistate VRC orchestration generates independent native paths,
corrections and scalar ROTD runs for reviewed physical spin/root contributions.
Stage-specific active spaces, averaged roots and theory controls are supported.
Physical fluxes are aligned by explicit energy offsets/levels and summed once;
correction variants are not counted as extra electronic states. MESS receives
unit electronic degeneracy for the summed flux. Existing single or two-TS
statistical/minimum models are used without changing the pinned ROTD_py code.

Legacy insertion of a zero-flux origin is retained. Negative native thresholds
are preserved, and the MESS core ground shift is included in explicit RRHO
ZeroEnergy values. No rule requires the original ROTD grid to contain zero.

## Actual NH2O supporting-file verification

All ten publisher TXT supporting files were obtained and checksum verified.
The actual MESS network and nine Ne tables were inspected. The derived report
`docs/NH2O_SI_verification.json` records the source hashes, four barrier
contracts, state/level counts, negative thresholds and an isolated canonical
electronic-population check with relative errors below 0.05%. The licensed
full SI is not vendored into the project. Local source files remain under
`/private/tmp/nh2o_si`; the verification entry point is
`tests/verify_nh2o_si.py`. This is not a reproduction of the paper's full rates.
The publisher-numbered Ne files could not be unambiguously matched to all
original filenames referenced in the MESS network.

NH2O's retained entrance complexes and state-specific normalization must not
be imposed on sulfur inputs. The current user-approved scope defers automatic
homolytic vdW-complex discovery. Existing inherited optional-complex policies
remain separate from electronic-state accounting.

## Required work on resume

1. **Finish native recon and scratch recovery.** Pull the checkpoint, select
   the compute-node scratch through the environment, retry known failed
   attempts, and review the native pilots before bulk scans. Investigate any
   native method/basis/root/record errors from their actual outputs. Keep
   completed WFU records on Blodgett. The default scratch capacity policy still
   needs a general, resource-aware improvement beyond its 4 GB floor; standard
   mount discovery, disk-demand estimates and quota-aware reporting are future
   operational improvements. A density-fitted diagnostic mode is an optional
   alternative requiring explicit theory labeling and its own native pilot;
   no silent DF substitution has been made.
2. **Full-PES MESS staging is unfinished.** `MESS._rotd_model` used to return
   a deferred PST model for `pes=1` even when required/production ROTD data were
   requested. Required/multistate full-PES configurations now fail explicitly
   before a misleading handoff. The sulfur examples use the supported `pes=0`
   workflow. To implement full-PES support, factor the direct renderer/stager
   into a shared routine; save hash-bound per-worker/per-UQ context (frequencies,
   rotors, stoichiometry, symmetry, physical electronic counting and fragment
   reference shifts); and render final Core Rotd blocks after global ROTD
   execution using the final network energy reference. Stage Ne files and
   provenance in the top-level `me` tree. Retain legacy auto/PST behavior when
   no required ROTD model was requested. Test scalar/multistate, UQ shifts,
   optional wells, mirrors and relative paths through actual PES assembly.
3. **Validate the multistate contract more deeply.** Before full-run handoff,
   also audit required-ROTD selection: the global production mode currently
   reaches retained optional vdW exits as well as selected homolytic channels;
   unselected exits need an explicit, consistent requirement/fallback policy.
   Audit VRC convolution against retained multi-conformer fragment unions;
   the current direct renderer uses selected fragment frequencies/rotors and
   must not be advertised as complete multi-conformer VRC support. The supplied
   examples use the default single-conformer TST model. Finally, the scalar
   multistate leaf currently chooses reference root 1 for spin splitting;
   preserve or expose an expert-reviewed reference-root override rather than
   assuming that choice for every possible excited-state input.
 Local fixture tests cover
   independent native/scalar templates, high-spin direct CBS references,
   selected roots, aligned flux sums, explicit levels, provenance changes and
   MESS rendering. Licensed native correction -> ROTD -> MESSPF/MESS execution
   is still required per sulfur channel. Check energy-zero conventions and
   physical-state/path counting against those native results. No coupled
   nonadiabatic dynamics has been implemented or assumed necessary.
4. **Review the chemistry from recon evidence.** Confirm parent/fragment
   geometries, conformers/orientations, active bonding/antibonding and sulfur
   nonbonding manifolds, larger-space and basis convergence, and crossings.
   CAS(13,8) / CAS(9,7) are provisional production requests. Select consistent
   WFU records and averaged roots for each basis/channel. Establish which roots
   actually contribute to reactive flux and their energy offsets; averaged
   roots are not automatically independent reactive channels. Resolved sulfur
   atomic levels are supplied; molecular CH3S JT/SO/vibronic levels, S2 fine
   structure and molecular SO corrections still require physical evidence.
   A missing correction is a labeled spin-free approximation, not an exact zero.
5. **Complete production accuracy checks.** Validate native MR-F12 energy
   selection, numerical gradients, state/reference quality, CASPT2 intruders,
   MRCI+Q size consistency and separated-fragment gaps. Converge CBS/theory,
   pivots/surfaces, Monte Carlo statistics and energy/angular grids. Validate
   collision/energy-transfer inputs and rotor/harmonic fallbacks. The example
   collider/Lennard-Jones values are user inputs, not fitted sulfur constants.
6. **Retain the compatibility audit's limits.** The unchanged origin/master
   tests produced 602 passes, 41 failures and 369 passing subtests; the older
   fixture assumptions are classified in `ANL_workflow_audit.md`. Do not call
   that suite passing or claim every master regression has been excluded.
   Current replacement tests require actual artifacts and bound scheduler
   generations. Preserve legacy scalar inputs, arithmetic/distance-switch gap
   policies, unknown uncertainties and basis naming behavior.
7. **Acetaldehyde diagnosis remains incomplete.** Its cluster restart killed
   jobs/controllers and left stale branch labels. HIR harmonic fallbacks and
   optional vdW failure were nonfatal. The earlier relaxed-scan record identifies
   Molpro exit 255 on the first force evaluation, but the actual native logs for
   that attempt have not been supplied. Do not attribute it to the recon brace,
   unit-label or scratch failure without its own native evidence. Its scoped
   stop/archive instructions preserve those files; no shared scratch cleanup
   or user-wide kill is authorized by inference.

## Verification at the stopping point

Before the final retry/guard changes, the full suite passed 1228 tests and
389 subtests, with the existing deprecated `imp` import and Sella external
cache-permission warnings. The focused retry/guard/campaign/MESS checks pass
72 tests and 41 subtests. The final full suite passes **1231 tests and 391 subtests**; all 158
non-template Python modules parse, and `git diff --check` passes. These tests simulate native schedulers/calculators except for
the actual Blodgett evidence explicitly described above. No licensed Molpro,
Slurm or MESS execution is available on the development machine.

## Procedure before the full CH3S/CH3SS runs

Finish/recover native pilots -> adaptive radial/orientation reconnaissance ->
collect and review evidence -> select active spaces, basis-specific records,
reactive roots and physical energy/state model -> native correction/ROTD/MESS
validation per channel -> convergence checks -> full single-well campaigns.
Do not run the production example JSONs before their reviewed state/surface
model files exist. Do not run required full-PES ROTD/MESS until item 2 is
implemented and tested. Automatic homolytic vdW-complex characterization and
full nonadiabatic dynamics remain explicitly deferred future scope.

No SSH or cluster job was launched by this development session. Further
implementation waits for the user's request to resume; the user's cluster
campaign can continue independently.
