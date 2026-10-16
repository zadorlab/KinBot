# Molpro confirmation: targeted resolution handoff

The minimum supported spaces remain **CH3S C–S (5,4), CH3SS C–S (7,5),
CH3SS S–S (7,5)**. The uploaded `_r2` run has **35 validated, two failed and two
waiting points**. All 833 exported files passed checksum verification. The
completed outputs and basis-specific orbital orderings have been independently
reviewed. Full KinBot production has not started.

The existing CH3S control resolves the large CASPT2 asymptotic splitting:
**equal-weight same-spin reference averaging reduces it from 1.1748 to
0.000144 kcal/mol, retaining IPEA=0.25 and real shift=0.3.** The new run confirms
this reference policy on the current geometries/bases and on S–S. It preserves
physical fragment excitation gaps. See the
[evidence, causal analysis and literature](docs/MOLPRO_CONFIRMATION_RESOLUTION_2026-10-16.md).

## Run the supplement on Blodgett

The user handles SSH. On the existing `composite` checkout:

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/orbital_reviews/molpro_confirmation/retry_unvalidated_on_blodgett.sh
```

This validates the source against the reviewed output/input/WFU hashes and
refuses a live source controller or allocation. It cancels nothing, deletes
nothing, and does not resubmit any old validated job. It prepares a sibling
`minimum_active_space_molpro_confirmation_20261016_r2_recovery` directory.
If the source changed since the reviewed export, the script stops before
preparing a new run so that newer results are not silently ignored.

**26 supplemental points, seven independent chains:**

- CH3S tilted MRCI: increase external Davidson capacity from 4 to 20; retain
  the method, CABS=1 and tolerances. Seed from the completed center WFU.
- S–S MRCI far→center→tilted: use documented perturbational CABS=2 consistently.
  The two old waiting points were blocked by the failed far predecessor.
- CH3S MRCI far/center CABS=2 controls: compare relative energies/spin gaps
  against completed CABS=1 results before accepting the approximation.
- CH3S and S–S CASPT2: five points in each of aDZ/aTZ with the shared-reference
  policy. Include far, center, inner, exact return and tilted points. No IPEA
  series or enlarged-space campaign. The CH3SS C–S policy is unchanged and
  its completed checks are reused.

All choices are explicit in
[recovery_plan.json](examples/anl/orbital_reviews/molpro_confirmation/recovery_plan.json).
The original catalogue is historical evidence; its selected-root CASPT2 policy
is superseded by this remedy. The default CABS setting in general KinBot is
unchanged. The alternative is confined to explicitly configured checks.

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
| 3 | aDZ/aTZ native Molpro confirmation | 35/39 reviewed; targeted remedies and changed-reference checks ready |
| 4 | General JSON controls | Orbital/state controls retained; explicit MAXDAV and MRCI CABS controls added/tested |
| 5 | General orbital reconnaissance | Procedure updated; chemistry remains example data |
| 6 | Complete production inputs | Counts configured; final orbital/surface/reference installation follows native closure |
| 7 | User-run end-to-end production support | Still to follow, including bug fixes through completion |
| 8 | Generality/integration audit | New controls and all 26 native-input transformation paths tested locally |
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

Local validation cannot establish native success of the unrun MAXDAV/CABS
remedies. The previous CH3S control does establish the reference-policy remedy
at its geometry; transfer across the current points is what the supplement
checks. These distinctions are part of the handoff, not a promise of zero errors.
