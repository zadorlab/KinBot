# CASPT2-only KinBot production handoff

The handoff preparation is complete. The two end-to-end PES calculations
have **not** been run locally or certified in advance. The user handles SSH
and launches them on Blodgett; subsequent failures/results will be reviewed
through completion. Repository cleanup remains deferred.

## Accepted inputs and scope

| Channel | CAS | Active natural orbitals | SA spatial roots per spin (D/Q) |
| --- | --- | --- | --- |
| CH3S C–S | (5,4) | 11–14 | 3/3 |
| CH3SS C–S | (7,5) | 18–22: S2 pi/pi-star + methyl radical | 1/1 |
| CH3SS S–S | (7,5) | 18–22: CH3S sulfur pair + departing S 3p | 6/6 |

These one-based indices refer to six supplied, hash-bound, basis-specific
natural-orbital records. The two CH3SS channels have different subspaces.
The retained sparse radial/angular/return evidence supports these minimum
spaces; it does not prove that no uncomputed geometry can change state character.

Production uses ordinary **single-bottleneck VRC**, CASPT2-F12 sampling at
aDZ-F12, a constrained relaxed CASPT2-F12/aTZ-F12 target-root scan, and
aDZ/aTZ spin gaps extrapolated with **exponent 3.7** and anchored to ground
high-spin ROHF-UCCSD(T)-F12b/aTZ-F12. **IPEA=0.25, real SHIFT=0.1 Eh** apply
to all new CASPT2 stages. MRCI is excluded. The search cutoff is 1000 kcal/mol.
ANL1-F12/CBH thermochemistry is a separate, higher-level workflow.

Every CASPT2 root uses the full declared same-spin reference list and common
weights. The general `energy_selection: {"4": "lowest"}` option evaluates
all declared quartet roots for the correction and selects their minimum:
it does not assume CASSCF root 1 is always the lowest correlated root.
All roots must succeed. It does not provide universal root-character tracking,
add XMS mixing, or guarantee a correct ground-state CC solution. Numerical
geometry gradients remain bound to explicitly selected target roots.

## Documentation and evidence

- [Full workflow, expected candidate chemistry, resources and commands](examples/anl/caspt2_production/README.md).
- [Every explicit parameter, nested control, consequential default and master/branch provenance](examples/anl/caspt2_production/PARAMETERS.md).
- [Maintained 12-item status checklist](docs/ORBITAL_RECON_DELIVERABLES.md).
- [Orbital occupations, chemical character and native evidence](docs/SULFUR_ACTIVE_ORBITAL_EVIDENCE.md).
- [Matched shift-0.1 native acceptance](docs/CASPT2_SHIFT01_NATIVE_REVIEW_2026-10-16.md).
- [General orbital-reconnaissance procedure](docs/ORBITAL_REVIEW_WORKFLOW.md).
- [Deferred cleanup inventory](docs/REPOSITORY_CLEANUP_PLAN.md).

Local verification: **571 tests + 31 subtests**, plus **five separate native
RotdPy/MESS/MESSPF checks**. The latter retain native files in
[the backend evidence package](examples/anl/caspt2_production/backend_smoke_review.json).
The new minimum-selector wrapper is regression-tested; its individual CASPT2
root calls have native confirmation, but this wrapper itself has not yet
executed in Molpro. No full sulfur VRC flux or MESS network is claimed complete.

Atomic sulfur fine-structure levels are forwarded to MESS/PF using its existing
electronic-level capability. The KinBot JSON forwarding and explicit term-
centroid energy hook were added on this branch; neither invents a VRC SOC
potential. No inner/outer stages or automatic equilibration are enabled.
The full reference explains the separate energy-zero and partition-function
uses and their tests. Transport parameters remain provisional example estimates.

## Blodgett launch

The helper prepares and checks the new inputs first, including small native
backend tests. Only then does it stop/archive the explicitly named old orbital
campaigns and launch the two PES controllers. It preserves native files and
records archived paths; unrelated jobs are untouched. Missing dependencies
or a failed preflight stop the transition before retirement.

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
export PATH=/home/hram/KinBot/.venv/bin:$PATH
export PYTHONPATH=/home/hram/KinBot${PYTHONPATH:+:$PYTHONPATH}
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
mkdir -p "$KINBOT_MOLPRO_SCRATCH"
/home/hram/KinBot/.venv/bin/python examples/anl/caspt2_production/handoff.py \
  launch /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016 \
  --retire-parent /home/hram/kinbot-runs
```

Each PES has one active worker and a 10-allocation cap: at most 20 allocations
across the two runs. All eleven known bad nodes are excluded:
`node[1,9-10,13-14,16,21,23-24,30,33]`. Actual scheduler exclusions remain a
runtime check. Preparation saves complete resolved JSON and input hashes.
The helper refuses duplicate launch; diagnose/resume an existing run rather
than deleting its results or repeating broad cancellation.

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/caspt2_production/handoff.py monitor /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016'
```

After launch, review optimized structures/IRC, accepted channels, ANL component
and formation-enthalpy closure, relaxed corrections/asymptotes, root/subspace
behavior, Monte Carlo convergence, grid-boundary minima and final MESS rates.
The 5% Monte Carlo target is statistical precision, not a total kinetic-error bound.
