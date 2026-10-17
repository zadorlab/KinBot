# CASPT2-only production handoff

This package prepares the two ordinary KinBot PES runs. It submits no more
orbital-reconnaissance or MRCI jobs. The user runs the commands on Blodgett;
the helper does not use SSH. Full production results are still to be obtained.

The first production attempt failed on an invalid generated Gaussian SCF
option list. That defect, worker failure propagation and monitoring are
repaired. The [repair audit](../../../docs/GAUSSIAN_PRODUCTION_REPAIR_2026-10-16.md)
records the native diagnosis and test coverage. Use the clean-restart command
below, which archives the failed attempt and starts a fresh `_r2` directory.

## Accepted electronic model

| Reaction | CAS | Natural active indices in the supplied records | CASSCF spatial roots |
| --- | --- | --- | --- |
| CH3S → CH3 + S | (5,4) | 11–14 | 3 doublets + 3 quartets |
| CH3SS → CH3 + S2 | (7,5) | 18–22 | 1 doublet + 1 quartet |
| CH3SS → CH3S + S | (7,5) | 18–22, different subspace from C–S | 6 doublets + 6 quartets |

These are one-based indices within each supplied natural-orbital checkpoint,
not transferable labels for arbitrary SCF orbitals. Six independent aDZ/aTZ
records are bound by hash and atom order. The C–S seed permutation is already
applied. Do not repeat it. The evidence is in
[orbital occupations and character](../../../docs/SULFUR_ACTIVE_ORBITAL_EVIDENCE.md),
[minimum-space decisions](../../../docs/SULFUR_MINIMUM_ACTIVE_SPACE_DECISIONS.md),
and [native recovery review](../../../docs/MOLPRO_RECOVERY_REVIEW_2026-10-16.md).

Use **IPEA 0.25 and real SHIFT 0.1 Eh** throughout new CASPT2 calculations.
The [matched native check](../../../docs/CASPT2_SHIFT01_NATIVE_REVIEW_2026-10-16.md)
completed both bases and all 24 requested energies. Its extrapolated spin-gap
change is 0.267019 kcal/mol using **exponent 3.7**. It supports the chosen
settings at the sensitive tested geometry; it is not a rate-error estimate.
New asymptotes are computed at the new settings. Old shift-0.3 energies are
retained as evidence, not mixed into production interaction potentials.

Every CASPT2 call uses the declared full same-spin `REFSTATE` and equal
`WEIGHT` values. This retains the native-confirmed shared-reference treatment
which removed the artificial sulfur-component splitting. CASSCF weights
are equal within each spin and equal between spin blocks. These orbital
weights are not statistical populations of reactive surfaces.

The quartet correction reference evaluates every declared correlated root
and selects the lowest energy for each channel. For CH3S this is necessary because the native center point has
root 3 lowest while other points have root 1. The general JSON option is
`energy_selection: {"4": "lowest"}`. It keeps all root energies in `KBGAPS`,
requires shared references, and rejects a geometry-gradient request for the
minimum policy. The one-root C–S channel adds no work; the six-root S–S
reference adds correlation calls only to the correction cut, not to angular
sampling. The relaxed target scan continues to differentiate its
explicit doublet root. This fixes a root-number assumption; it does not
guarantee universal orbital/state tracking or CC reference quality.

## What these runs do

1. Explore each PES with UMA-s-1p2/OMol at L1; refine accepted structures,
   transition states, modes and rotors with the B2PLYP-D3(BJ)/cc-pVTZ L2
   profile. Conformer, TS and IRC checks determine accepted chemistry.
   The search barrier cutoff is **1000 kcal/mol**, with all reaction families
   enabled. Only the three explicitly reviewed initial-well homolytic bonds
   receive VRC treatment; newly found isomers do not silently reuse their
   orbital records merely because the element counts match.
2. Perform automatic CBH/ATcT/ANL1-F12 thermochemistry on the accepted network.
   The ANL reference uses augmented QZ-F12/5Z-F12 bases and its existing
   higher-order correlation, core-valence, scalar-relativistic, DBOC,
   state-specific spin-orbit and zero-point/anharmonic components. The
   existing size policy and component validation apply. This is separate
   from the cheaper VRC aDZ/aTZ electronic model. ANL1 is substantial work;
   completing the fixed-geometry confirmation did not complete these tasks.
3. For each reactive doublet component, run a constrained relaxed
   CASPT2-F12/aTZ-F12 cut, using numerical gradients. Evaluate matching
   frozen-fragment aDZ sampling energies and the trusted relaxed energies:

   `E_trusted(D_i) = E_ROHF-UCCSD(T)-F12b/aTZ(Q_ground) + Δ_CBS(D_i-Q_ground)`

   `Δ_CBS = (3^3.7 Δ_aTZ - 2^3.7 Δ_aDZ)/(3^3.7 - 2^3.7)`.

   Retain the aTZ-only gap and CBS-minus-aTZ difference as diagnostics.
   The relaxed/frozen difference includes geometry relaxation at the stated
   CASPT2 scan level. It is not a gradient optimization of the composite
   trusted energy. The one-dimensional correction is added to angular
   sampling; it cannot repair arbitrary angular errors.
4. Run **ordinary single-bottleneck VRC-TST**, with reactive-atom pivot
   surfaces over the configured distance grid. CH3S includes its two
   reactive ground-2E spatial components; each CH3SS bond uses its ground
   doublet component. Extra SA roots and quartet reference roots are not
   all counted as reactive flux. This spin-free adiabatic assignment is a
   physical model assumption. It does not assert nonadiabatic equilibration.
   The inputs contain no inner/outer stages, new complex well, or automatic
   distance-dependent spin-orbit model. The optional general two-bottleneck
   hook is unused. If the resulting rate-controlling geometry requires a
   different kinetic model, that requires a separate evidence-based review.
5. Integrate the corrected aDZ CASPT2 surface with RotdPy to a requested
   **5% Monte Carlo statistical target**, subject to its explicit sample
   ceilings. Apply the existing 0.85 local recrossing correction once.
   Inspect convergence and distance-grid boundary minima; neither a sample
   ceiling nor the preset establishes total kinetic uncertainty.
6. Feed validated state-count tables, accepted species/TS thermochemistry,
   electronic levels and rotor models into MESS/PF and MESS. Conventional
   accepted saddle pathways keep their ordinary RRHO/HR/TST treatment.
   VRC failure is not silently replaced by a PST success. The actual network
   and rates depend on completed optimization, IRC, ANL and sampling results.

The read-only initial search finds eight candidate searches for CH3S and
nineteen for CH3SS. They include H migrations, insertion/addition and shift
families, H2 elimination, and the requested homolytic scissions. CH3S
isomerization toward CH2SH and decomposition toward CH3+S are among the
chemically relevant possibilities. CH3SS includes both CH3+S2 and CH3S+S.
The [complete candidate family/motif inventory](initial_reaction_candidates.json)
is discovery evidence, not accepted product or barrier predictions.

## Sulfur levels and parameter provenance

MESS already supports electronic levels. The branch adds the general KinBot
JSON matching/forwarding hook. The atomic S list `(0,5), (396.055,3),
(573.640,1)` cm^-1 gives thermal populations. Its explicitly declared
term-centroid convention supplies **-195.756111 cm^-1** to the atomic ground-J
energy in ANL. These are two separate uses of the data, with duplicate
component protection. The code does not add a new SOC potential to VRC.
Molecular spin-vibronic spectra are not invented from term degeneracies.

[PARAMETERS.md](PARAMETERS.md) explains every explicitly supplied parameter,
every nested state/orbital/surface control, consequential defaults, and which
keys existed in master versus this branch. Transport epsilon/sigma and the
default energy-relaxation model are provisional example choices, not measured
or fitted evidence for these systems. They matter for pressure-dependent rates.
The 1000 kcal/mol search cutoff is separate from MESS's 400 kcal/mol model
energy limit. See also the [backend/interface audit](../../../docs/CASPT2_PRODUCTION_TRANSITION_2026-10-16.md).

## Prepare and launch on Blodgett

These commands archive native files and restart cleanly. Preparation checks/copies
all six checkpoints before any old work is stopped. Preflight renders 36
Gaussian route previews from the actual input profiles, checks installed
dependencies, and executes small native RotdPy/MESS/MESSPF tests, without QC
submissions. A missing dependency or failed check stops the transition.
The environment must include `pytest` to execute those small preflight checks.
Preparation also saves and hashes each complete `resolved_parameters.json`.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/caspt2_production/restart_on_blodgett.sh
BASH
```

The script calls `handoff.py launch NEW --replace-production OLD` with the
original production directory and a fresh `_r2` destination. It inventories
both old PES trees, stops their detached controllers and scoped Slurm jobs,
waits for teardown, then archives the common parent as `.stopped.*`.
It launches fresh detached controllers on the login host. It never uses
user-wide cancellation or deletes shared scratch. Retirement JSON records
old/new locations. Historical absolute manifest paths remain provenance;
these archives are not in-place resumable after the move. Seeds are copied
independently. The old orbital campaigns were already archived on the first
launch; they are not processed again. `prepare`/`preflight` actions allow
preparation/checks separately, without retirement or submission.

Each PES has one active KinBot worker and a shared **10-allocation cap**,
giving at most **20 allocations across these two runs**, excluding unrelated
user work. Eight processors are requested per QC/sampling allocation;
Rotd jobs request exclusive nodes. All known bad nodes are excluded:
`node[1,9-10,13-14,16,21,23-24,30,33]`.
The controller can use fewer allocations while dependencies are pending.
The actual scheduler's `ExcNodeList` remains the runtime verification.

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/caspt2_production/handoff.py monitor /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r2'
```

```bash
/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/caspt2_production/handoff.py \
  monitor /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r2 --follow
```

The helper refuses a duplicate launch after `launch.json` exists. `--all`
expands task/reaction detail and `--json` exports a machine-readable snapshot.
The monitor reads actual PES/worker/native logs and execution records, rather
than relying on driver stdout or Slurm accounting. See the
[general monitor guide](../../../docs/RUN_MONITOR.md).
Shells/processes on Blodgett remain under the user's control.

## Acceptance after launch

The local checks and native electronic evidence justify starting production;
they do not certify completed kinetics. Review the actual ANL component and
formation-enthalpy closure, accepted TS/IRC topology, relaxed cuts and their
asymptotes, root/subspace behavior, sampling convergence, boundary minima,
MESS/PF thermochemistry and final rate tables. A current-file warning is
diagnostic evidence, not automatically a terminal failure. Conversely, an
empty queue or successful parser is not scientific acceptance.

Local handoff verification: **571 tests and 31 subtests passed**, covering
input routing, native input rendering, gradients/energy selectors, CBS gaps,
checkpoint installation, admission limits, scoped retirement, thermochemistry
and MESS integration. Separately, **five native RotdPy/MESS/MESSPF tests passed**:
analytic zero-potential state counts, MESS single/statistical/minimum table
handling, and atomic sulfur fine-structure partition populations. The optional
two-bottleneck checks exercise an existing hook; these examples use only single.
Native outputs and hashes are in [backend_smoke_review.json](backend_smoke_review.json)
and `backend_smoke.evidence.tar.gz`. These are backend checks, not sulfur kinetics.
The lowest-manifold wrapper has rendered-input regression coverage; its
constituent CASPT2 root calls were confirmed natively, but the new wrapper has
not yet run in Molpro. The actual Blodgett dependencies are checked again at launch.

The subsequent full production-repair audit passed **1,748 tests and 424
subtests**, with zero failures/skips and the native backend checks enabled.
See [repair_validation.json](repair_validation.json) and the
[audit report](../../../docs/GAUSSIAN_PRODUCTION_REPAIR_2026-10-16.md).
The earlier focused counts above are historical, not additional unique coverage.

The maintained [12-item checklist](../../../docs/ORBITAL_RECON_DELIVERABLES.md)
distinguishes completed handoff work from the user's production launch and
the deferred repository cleanup. The latter remains deferred until both
examples complete, as requested.
