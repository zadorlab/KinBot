# First production attempt: input repair and clean replacement

The first production attempt (commit `c340cc89`) stopped during both initial
L2 refinements at 18:35 PDT on 16 October 2026. No production ANL, VRC flux,
or MESS calculation had completed. This was a generated-input defect, not
evidence that node20 was unhealthy. Both initial UMA optimizations and the
first conformer search had progressed before the Gaussian failure.

## Native evidence and cause

The user supplied the generated `.com` files and Gaussian 16 Rev. B.01 logs
for `470831050000000000002_well_high` and
`791471490990000000002_well_high`. Both contained:

```text
scf(xqc,(xqc,conver=10))
QPErr --- A syntax error was detected in the input line.
Last state= "SCF2"
Error termination via Lnk1e
```

The ANL convergence-profile integration in `13a28bba7` supplied the SCF
value `(xqc,conver=10)`, including outer parentheses. The older QC option
merger split on commas without removing those parentheses. Combining it
with the base `xqc` produced the invalid anonymous nested option group.
The stderr segmentation message followed this explicit Link 1 parse error;
it does not establish a separate hardware fault.

The generated route is now `scf(xqc,conver=10)`. The generic Gaussian merger
accepts bare or outer-parenthesized lists, splits only at top-level commas,
and overrides suboptions case-insensitively. The ASE writer uses the same
normalization, covering direct profile calculators as well as generated QC
workers. Method/basis, checkpoint and other scalar fields remain scalar
overrides. Invalid anonymous groups/unbalanced lists raise before Gaussian
is invoked. Both ANL profile generators now supply bare SCF option lists.
No energy model, basis, active-space, spin-root, or convergence value changed.

## Other defects addressed

- PES previously treated a dead worker process as “Finished” and attempted
  postprocessing with a missing high-level energy. Workers and PES now record
  start/final status, PID, host, exit code and failure reason in
  `kinbot_status.json`/`pes_status.json`. PES requires the matching successful
  worker record before accepting that worker. Child processes are reaped
  with `Popen.poll()` rather than being mistaken for running zombie PIDs.
- The old monitor mostly printed a license banner from `driver.log`. The
  general [run monitor](RUN_MONITOR.md) reads actual PES/worker/native logs,
  scheduler resources/exclusions, QC exit markers, reaction/optimization
  progress, dispatch dependencies and downstream execution records. It
  works when Slurm accounting is disabled. Failed scheduler queries are
  reported as unknown, not an empty queue or successful completion.
- Production consists of two marked PES directories under an unmarked
  parent. Clean replacement inventories both, stops their processes and
  allocations, rechecks both, and only then archives the whole parent.
  It does not pass the unmarked parent to single-run cleanup. Replacement
  preparation and preflight happen before the first stop. A failure leaves
  remaining old files in place and prevents the new launch.
- The full suite caught two interactions in the new repair: PES wrote its
  status into an incompatible legacy directory before rejecting it, and
  rendering the example previews leaked `SBATCH_EXCLUDE` into later callers.
  Compatibility now precedes the first status write, and preview construction
  restores its incoming environment. Both regression checks pass.
- A separate, older `thread_kinbot` launcher contained a developer's absolute
  installation path, Python-2-style process-output handling, and undrained
  output pipes. It now launches the current Python modules, writes logs and
  polls its child processes. The production PES launcher also uses its own
  Python interpreter without a literal shell `&` argument. These are general
  launch repairs, not molecule-specific behavior.
- A whole-tree syntax check found Python 2 exception syntax in four legacy
  NWChem templates. Those clauses now use Python 3 syntax. A new regression
  renders and parses all 33 Python worker templates; all 170 non-template
  runtime Python modules also parse. These checks do not validate native
  NWChem execution, which is outside the selected production workflow.

## Scope of the audit

The due-diligence pass combined source inspection of the runtime interfaces
below, searches for example/site-specific implementation, and the **entire
discoverable `tests/` suite**. Final result: **1,748 tests and 424 subtests
passed; zero failures or skips** (80.29 seconds). Native MESS and MESSPF
executables were supplied, so their optional checks ran. The suite also ran
the actual RotdPy free-pair integrator check. Two warnings were emitted.
The reproducible command and source hashes are recorded in
[repair_validation.json](../examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/production_checks/repair_validation.json).

| Workflow area | Checks and relevant coverage |
| --- | --- |
| JSON/theory routing and Gaussian | Explicit values and profile overrides; generated well/TS/IRC/rotor/scan routes; force/frequency calls; scratch and worker exit handling (`test_parameters`, `test_theory_profiles`, `test_gaussian_input_routes`, `test_gaussian_runtime`, `test_qc_jobs`). |
| PES, conformers and reaction lifecycle | Worker completion contract; stereochemical identity; conformer/counting records; optimization and Hessian provenance; reaction selection, IRC endpoints and hindered-rotor restart/acceptance tests. |
| CBH/ANL | Task graphs, component parsing and energy provenance, CBS formulae, CBH/ATcT handoff, scheduler dispatch, formation enthalpies and thermochemistry/MESS integration (`test_anl_*`). |
| Orbital reconnaissance and Molpro | General state-model validation, immutable seed hashes and basis records, full declared same-spin references, correlated-energy selection and solver settings; campaign/retry and rendered-input tests. Existing native confirmations remain the electronic evidence. |
| VRC corrections | Per-reaction model selection, fixed selected-root gradient contract, high-spin anchor and gap construction, configured CBS exponent, asymptote handling and branch dispatch (`test_vrc_*`, `test_kb_vrc_selection`, `test_molpro_ase`). |
| RotdPy/MESS | Profile/table routing, result provenance and completion gates, electronic-level forwarding, native analytic RotdPy normalization and MESS/MESSPF parser/partition checks. Optional two-bottleneck hooks are tested but disabled in these examples. |
| Resources, monitoring and retirement | Allocation admission and exclusions, unavailable scheduler behavior, execution exit markers, bounded log following, scoped process/job cleanup and both-PES archive sequencing (`test_slurm_retry`, `test_run_*`, `test_caspt2_handoff`, `test_driver_status`). |

The general Python/template source search found no CH3S/CH3SS names, example
chemids, user-home paths, Blodgett hostname or bad-node list after the legacy
launcher repair. Inspection of active-space/state handling found the example
choices supplied by models/checkpoints, not a molecule-specific dispatch.
Example geometries, orbital indices and site defaults intentionally remain in
`examples/`. Physics constants, native record-number conventions, conservative
resource guards and documented configurable presets are still constants;
this is not a claim that all numerical literals have been eliminated.

No scientific default was changed in this repair. The production models retain
the reviewed (5,4)/(7,5)/(7,5) spaces, IPEA=0.25, real shift=0.1, gap CBS power
3.7, CASPT2-only correction, ordinary single-bottleneck VRC and 1000 kcal/mol
search cutoff. Stale parameter comments describing an automatic full-valence
choice and singlet/triplet-only correction were corrected to match the code.

Regression coverage renders real Gaussian input, not only intermediate
parameter dictionaries: all four automatic ANL tiers; well, TS, IRC, rotor,
rigid rotor, scan and Gaussian VRC argument variants; force and fixed-Hessian
calls; and the generated high-level worker script. Tests also cover the
actual two production input files and clean-restart success/failure paths.

Cluster preflight now writes 36 input previews under `input_preflight/`
using the actual two input/profile combinations (six argument variants,
energy/force/Hessian), saving their hashes in `preflight.json`. This is input
rendering, not execution of Gaussian or optimization on the login node.
The existing small native RotdPy/MESS/MESSPF preflight still runs.

Gaussian is not installed in the local test environment. Native acceptance
of the repaired route must therefore occur on Blodgett. Passing the local
suite does not certify every unexecuted backend, convergence, scientific
result, or future cluster behavior. The first attempted full production
run exposed a missing rendered-route regression; it was not valid evidence
of end-to-end success.

The clean restart runs from a fresh `_r2` directory. Old native evidence is
archived, with no selective reuse of the early results. The detached Python
PES/worker controllers run on the host executing the launch command, which
should be the login host. They request no Slurm allocation. The existing
10-allocation-per-PES admission budgets, 8-core QC requests and eleven-node
exclusion list remain in effect.
