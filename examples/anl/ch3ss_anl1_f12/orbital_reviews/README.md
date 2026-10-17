# CH3SS orbital review and production input

Use [../ch3ss.json](../ch3ss.json) as the single KinBot input. All chemistry
settings, per-bond CAS dimensions and indices, basis/state/solver choices,
resource limits and the **1% RotdPy flux target** are explicit there.
The 4 native aDZ/aTZ checkpoint files under `orbitals/` are tracked data;
no state-model JSON, installer or configuration in another folder is needed.
Copying this example directory preserves all production dependencies.

- [Active-space justification, occupations, ordering and native evidence](ACTIVE_SPACE.md).
- [Every explicit parameter and its provenance](PARAMETERS.md).
- [Workflow, recovery, Blodgett restart and monitoring](RUNNING.md).
- [Reviewed seed inventory](orbital_inventory.json) and [native-file inventory](native_inventory.json).
- [Candidate reaction inventory](initial_reaction_candidates.json), which is
  an enumeration, not a claim of accepted reaction paths.

Only the main JSON is at the example root. The old reconnaissance recipe
and historical reviews are retained here as evidence, not read by production.
The shared two-PES utility is [../../tools/run_examples.py](../../tools/run_examples.py);
it stages the declared inputs and binary data, checks dependencies, and
handles scoped cleanup/launch. It does not select chemistry or add hidden settings.

The current preset is Gaussian B3LYP-D3(BJ)/cc-pVDZ at L1 and
B2PLYP-D3(BJ)/cc-pVTZ at L2, with tight force criteria and bounded minimum
recovery. The full restarted native workflow, ANL thermochemistry, RotdPy
flux convergence and MESS network still require execution and review.
