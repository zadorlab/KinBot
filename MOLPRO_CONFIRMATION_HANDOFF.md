# Molpro confirmation: targeted resolution handoff

**11:42:40 PDT CABS=3 monitor:** CH3S far is validated; S–S far is still
running with recent output and no reported failure. This is successful
completion of the six-root CH3S point according to the controller; the
coupled-CABS energies still require independent native review. The S–S
twelve-root point remains unresolved. Preserve it and the recovery results.

**Latest native report:** both CABS=2 far points failed after converged MRCI.
The proposed perturbational remedy is not accepted for production. The two
isolated tests of the documented coupled treatment, CABS=3, have been submitted;
S–S native completion remains pending. Preserve the recovery calculations
and all completed results, including MAXDAV20. See the updated
[resolution record](docs/MOLPRO_CONFIRMATION_RESOLUTION_2026-10-16.md).

**10:12:32 PDT monitor:** 19/26 validated, two failed, four waiting, one
running. The active point is S–S aTZ shared-reference CASPT2 `return`; its
`tilted` successor is waiting. The other three waiting points depend on the
failed CABS trials. This accounting indicates MAXDAV20 finished; its native
energies still need export review. The monitor validates completion, not the
scientific acceptance of all corrected interaction energies.

The minimum supported spaces remain **CH3S C–S (5,4), CH3SS C–S (7,5),
CH3SS S–S (7,5)**. The uploaded `_r2` run has **35 validated, two failed and two
waiting points**. All 833 exported files passed checksum verification. The
completed outputs and basis-specific orbital orderings have been independently
reviewed. Full KinBot production has not started.

The existing CH3S control resolves the large CASPT2 asymptotic splitting:
**equal-weight same-spin reference averaging reduces it from 1.1748 to
0.000144 kcal/mol, retaining IPEA=0.25 and real shift=0.3.** The running supplement
tests this policy on the current geometries/bases and on S–S, including the
requirement to preserve physical fragment excitation gaps. See the
[evidence, causal analysis and literature](docs/MOLPRO_CONFIRMATION_RESOLUTION_2026-10-16.md).

## Existing CABS submission (already launched)

The user handles SSH. The launch commands below are retained for reproducibility.
**Do not rerun them while the CABS tests are running.** Use the monitor below.

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/orbital_reviews/molpro_confirmation/probe_cabs3_on_blodgett.sh
```

This validates the source against the reviewed output/input/WFU hashes and
refuses a live source controller or allocation. It cancels nothing, deletes
nothing, and does not resubmit any old validated job. It prepares a sibling
`minimum_active_space_molpro_confirmation_20261016_r2_cabs3_gates` directory.
If the source changed since the reviewed export, the script stops before
preparing a new run so that newer results are not silently ignored.

The two jobs are CH3S C–S and CH3SS S–S at 30 Å, aTZ-F12, respectively
six and twelve correlated roots. Their only method change relative to the
failed CABS=2 far-point trials is `CABS_SINGLES=3`. Original seeds, orbitals,
references, convergence thresholds and relaxed-Davidson selection remain.
No center/tilted calculations are queued by this plan. The coupled option
is documented, but the manual does not promise it avoids the internal error.
The exact `strange..` guard cannot be diagnosed from the available output.

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py monitor /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_cabs3_gates --roots'
```

`--roots` streams native output to show the number of selected correlated
energies printed, the current call's declared multiplicity/root and the latest
iteration, even when the program header lies outside the normal log tail.
It ignores inherited RESTART variables and does not treat pre-CABS MRCI or
Davidson numbers as completed selected energies. These are progress markers,
not a convergence verdict or time-to-finish estimate. Output scans are optional;
the default monitor keeps its bounded log-tail reads.

Every detailed monitor now gives per-chain completion counts. To list every
successful repair point, recorded root count and waiting predecessor:

```bash
/home/hram/KinBot/.venv/bin/python \
  examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py monitor \
  /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_recovery --all --roots
```

This command reads files and scheduler state only; it changes no running job
or result. Completed-point summaries come from manifests without rescanning
all their native outputs. To inspect the last lines directly while S–S runs:

```bash
tail -n 40 -F /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_cabs3_gates/ch3ss_ss_atz_mrci_cabs3_gate/atz__far/atz__far.out
```

The recovery's successful native results must still be exported using `collect`
below and transferred locally for review. At this update only the original
`_r2` native export is present locally. Required review covers all repaired
CASPT2 asymptotes, angular/radial energies, same-geometry return orbital spaces,
and the MAXDAV20 MRCI energy selection. A successful job count alone does not
answer those scientific checks. Collecting one stopped run is independent of
the other CABS run's live allocations.

For native review after the two tests stop:

```bash
/home/hram/KinBot/.venv/bin/python \
  examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py collect \
  /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_cabs3_gates \
  --output /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_cabs3_gates.evidence.tar.gz
```

The two-test cap is separate from the running recovery. Its two failed CABS
chains cannot submit dependents, leaving at most five independent chains
progressing there. No new global scheduler/controller feature is introduced.
If the gates succeed, review the final energies and then prepare consistent
interaction checks; do not mix CABS=1 and CABS=3 dissociation zeros. If they
fail, retain the native reproducer for Molpro implementation diagnosis.

## Existing supplement: leave it running

**26 supplemental points, seven independent chains, already launched:**

- CH3S tilted MRCI: increase external Davidson capacity from 4 to 20; retain
  the method, CABS=1 and tolerances. Seed from the completed center WFU.
- S–S MRCI CABS=2: far failed; center and tilted remain dependency-blocked.
- CH3S MRCI CABS=2 control: far failed; center remains dependency-blocked.
- CH3S and S–S CASPT2: five points in each of aDZ/aTZ with the shared-reference
  policy. Include far, center, inner, exact return and tilted points. No IPEA
  series or enlarged-space campaign. The CH3SS C–S policy is unchanged and
  its completed checks are reused.

All choices are explicit in
[recovery_plan.json](examples/anl/orbital_reviews/molpro_confirmation/recovery_plan.json).
The original catalogue is historical evidence; its selected-root CASPT2 policy
is superseded by this remedy. The default CABS setting in general KinBot is
unchanged. The alternatives are confined to explicitly configured checks.
The latest 10:12:32 monitor counted 19 validated, 1 running, 4 waiting and
2 failed. Keep the last CASPT2 chain running; do not relaunch this supplement.

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py monitor /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_recovery'
```

For a stopped completed/failed supplement, export with:

```bash
/home/hram/KinBot/.venv/bin/python \
  examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py collect \
  /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_recovery \
  --output /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_r2_recovery.evidence.tar.gz
```

An existing export is never overwritten. Transfer it yourself for the native
closure review. The original `_r2` monitor will retain its historical failures;
use the recovery monitor for the new work. Successful new outputs do not
silently replace old method-labelled energies or automatically approve production.

Eight ranks retain 1000 MW/rank and a **104000 MB/node reservation**, below the
reported 126000 MB node limit. This is a reservation/budget check, not a measured
peak-memory guarantee. The seven-chain concurrency cap and all known node
exclusions are explicit. Infrastructure retries remain bounded and recorded.
CASSCF MAXIT stays 999; correlation MAXIT/MAXITI stay 2000, MAXVI=60,
THRDEN=1e-8 and THRVAR=1e-10. Only the stalled external solve adds MAXDAV=20.
`CLUSTER=0; ENERGD(1)` still selects relaxed Davidson after CABS completion.

## Deliverable status

| # | Deliverable | Current status |
|---|---|---|
| 1 | Literature-grounded minimum active spaces | Local choices complete; unchanged by this method diagnosis |
| 2 | ORCA chemical/radial/angular/return checks | Complete; 47 curated points retained |
| 3 | aDZ/aTZ native Molpro confirmation | 35/39 original points reviewed; 19/26 recovery points in last recovery snapshot; CABS3 CH3S complete/S–S running; new native review pending |
| 4 | General JSON controls | Orbital/state controls retained; explicit MAXDAV and MRCI CABS controls added/tested |
| 5 | General orbital reconnaissance | Procedure updated; chemistry remains example data |
| 6 | Complete production inputs | Counts configured; final orbital/surface/reference installation follows native closure |
| 7 | User-run end-to-end production support | Still to follow, including bug fixes through completion |
| 8 | Generality/integration audit | Existing controls retained; 26 recovery and 2 CABS3 native-input transformations tested locally |
| 9 | Reproducible evidence and conclusions | Original outputs/WFUs retained; hashed review and orbital map committed |
| 10 | Repository cleanup inventory | Recorded; actual cleanup remains deferred |
| 11 | Consolidated handoff | Targeted supplement above; production approval is not yet claimed |
| 12 | Preserve earlier work | All original runs retained; no cancellation/deletion in this supplement |

The [maintained checklist](docs/ORBITAL_RECON_DELIVERABLES.md) remains the full
scope. The [general procedure](docs/ORBITAL_REVIEW_WORKFLOW.md) can produce the
same declared checks; these scripts are example adapters, not molecule-specific
logic in the general runtime. After the remedies are demonstrated and the
corrected surfaces reviewed, install the accepted native settings and proceed
with the ordinary KinBot inputs through CBH-ANL, RotdPy and MESS.

## Artifacts

- [Native review](examples/anl/orbital_reviews/molpro_confirmation/native_r2_review.json)
  and [basis-specific orbital records](examples/anl/orbital_reviews/molpro_confirmation/native_orbital_orderings.json).
- The repository-root `molpro_confirmation_transfer_20261016.tar.gz` and SHA
  sidecar provide an optional transfer snapshot. Pulling the repository is
  sufficient; do not extract over the running tree.
- The full new Molpro evidence archive remains in `examples/anl/orbital_reviews`,
  SHA-256 `4c80871bfb7ddb8f20d930a076bfabe8906bdba0df10b8d76c168f47ac2625fb`.
- The full ORCA archive is retained under
  `local_runs/qchem_casscf_20261015/reproducibility/orca_minimum_space_final_20261016.tar.gz`,
  SHA-256 `13762a1124092834313eee3bf463383f9c5658535351715dc4efdf981feca523`.

Local validation cannot establish native success of the CABS=3 tests. The
new recovery monitor indicates MAXDAV completion, but its native energies
remain to be reviewed. The previous CH3S control establishes the reference-policy remedy
at its geometry; transfer across the current points is what the supplement
checks. These distinctions are part of the handoff, not a promise of zero errors.
