# ANL component version 2 and fresh r7 handoff

This version changes new ANL graphs. Existing version-1 graphs keep their
original definitions and are not relabeled as version 2. The r7 handoff starts
fresh chemistry and copies only the reviewed VRC orbital seeds/settings.

## Complete component matrix

`VnZ` means cc-pVnZ; `VnZ-F12` means cc-pVnZ-F12; `wCVnZ` means
cc-pwCVnZ. Prefix `a` means `aug-`, including hydrogen. `CVTZ-DK` means
cc-pCVTZ-DK. The literal input knob is `anl_augmented` (boolean), not
`augmented-ANL`. `null` inherits the older `anl_f12_augmented` option.

| Variant | anl_augmented | Geometry and harmonic surface | Reference electronic energy | Core-valence pair | Scalar-relativistic bases | DBOC | Higher excitations |
|---|---|---|---|---|---|---|---|
| ANL0 | false | CCSD(T)/VTZ | CCSD(T) CBS(VQZ,V5Z) | CCSD(T), wCVTZ/wCVQZ | CVTZ-DK | HF/VDZ | MRCC VDZ |
| ANL0 | true | CCSD(T)/aVTZ | CCSD(T) CBS(aVQZ,aV5Z) | CCSD(T), awCVTZ/awCVQZ | aCVTZ-DK | HF/aVDZ | MRCC aVDZ |
| ANL0-F12 | false | CCSD(T)-F12b/VTZ-F12 | mean of F12b CBS(VTZ-F12,VQZ-F12) and conventional CBS(VQZ,V5Z) | F12b, wCVTZ/wCVQZ | CVTZ-DK | CCSD/VDZ | MRCC VDZ |
| ANL0-F12 | true | CCSD(T)-F12b/aVTZ-F12 | mean of F12b CBS(aVTZ-F12,aVQZ-F12) and conventional CBS(aVQZ,aV5Z) | F12b, awCVTZ/awCVQZ | aCVTZ-DK | CCSD/aVDZ | MRCC aVDZ |
| ANL1 | false | CCSD(T)/VTZ then VQZ; TZ/QZ harmonic extrapolation | CCSD(T) CBS(V5Z,V6Z) | CCSD(T), wCVQZ/wCV5Z | wCVQZ-DK/wCV5Z-DK | HF/VTZ | MRCC VDZ/VTZ |
| ANL1 | true | CCSD(T)/aVTZ then aVQZ; TZ/QZ harmonic extrapolation | CCSD(T) CBS(aV5Z,aV6Z) | CCSD(T), awCVQZ/awCV5Z | awCVQZ-DK/awCV5Z-DK | HF/aVTZ | MRCC aVDZ/aVTZ |
| ANL1-F12 | false | CCSD(T)-F12b/VTZ-F12 then VQZ-F12; TZ/QZ harmonic extrapolation | mean of F12b CBS(VQZ-F12,V5Z-F12) and conventional CBS(V5Z,V6Z) | F12b, wCVQZ/wCV5Z | wCVQZ-DK/wCV5Z-DK | CCSD/VTZ | MRCC VDZ/VTZ |
| ANL1-F12 | true | CCSD(T)-F12b/aVTZ-F12 then aVQZ-F12; TZ/QZ harmonic extrapolation | mean of F12b CBS(aVQZ-F12,aV5Z-F12) and conventional CBS(aV5Z,aV6Z) | F12b, awCVQZ/awCV5Z | awCVQZ-DK/awCV5Z-DK | CCSD/aVTZ | MRCC aVDZ/aVTZ |

The reference average retains equal weights, independent extrapolations and
the existing exponent 3.7. F12b uses `scale_trip=1` consistently. Core-valence
is CBS(all-electron) minus CBS(frozen-core), with identical orbital/auxiliary
bases and geometry on the two arms. It uses weighted core-valence **orbital**
bases; it does not invent a `pwCVnZ-F12` basis family.

For F12 core-valence the explicit fitting sets are
`aug-cc-pwCVnZ/mp2fit` for `df_basis`, `aug-cc-pV5Z/jkfit` for
`df_basis_exch`, and `aug-cc-pwCV5Z/mp2fit` for `ri_basis` (large CABS),
with n=T/Q/5 as appropriate. These are shared between AE/FC arms. Override
them with `anl_cv_f12_auxiliary`, keyed by `T`, `Q`, `5`, each containing
`df_basis`, `df_basis_exch`, and/or `ri_basis`. There is no silently guessed
weighted-CV OptRI alias. Available orbital and fitting functions are
element/program-version dependent; a missing basis is reported, never
silently replaced by a lower cardinal.

Scalar relativity is **all-electron conventional CCSD(T)** with DKH2 minus
nonrelativistic CCSD(T), at the same DK basis on both arms. For ANL1 the QZ/5Z
differences are extrapolated with exponent 3.7. The native decks put
`set,dkho=2` or `set,dkroll=0` after the basis declaration and before RHF;
the latter prevents a DK basis from automatically enabling DKH on the NR arm.
This conventional component is necessary: Molpro explicitly disallows F12
with DKH/X2C. See the [F12 manual](https://www.molpro.net/manual/doku.php?id=explicitly_correlated_methods)
and [relativistic manual](https://www.molpro.net/manual/doku.php?id=relativistic_corrections).

Higher excitations retain the MRCC unrestricted coupled-cluster equations
with restricted HF references, and conventional matched orbital bases:
ANL0 adds CCSDT(Q)−CCSD(T) at DZ; ANL1 adds that difference at TZ and
CCSDTQ(P)−CCSDT(Q) at DZ. These terms are not F12. Experimental spin-orbit
corrections and their reference conventions are unchanged. One/two-electron
rank-exact reductions remain explicit, including zero higher excitations.

Closed-shell references are RHF; open-shell Molpro/MRCC CC uses ROHF-UCC
(`rhf` is Molpro's open-shell restricted alias). HF DBOC uses RHF/ROHF.
**The sole CC-reference exception is open-shell CCSD DBOC, which uses UHF**
because CFOUR does not support ROHF CCSD DBOC. Closed-shell CCSD DBOC uses
RHF. This is the specifically authorized exception, supported by the
[CFOUR DBOC manual](https://cfour.uni-mainz.de/cfour/index.php?n=Main.CalculationOfDiagonalBorn-OppenheimerCorrections).

L1/L2 Gaussian DFT remains RKS for closed shells and UKS for open shells,
as requested to avoid a separate RO→UKS VPT2 geometry optimization. The
Gaussian VPT2 correction remains frequency-only on the accepted matching L2
surface, with the existing symmetry recovery. No extra UKS optimization is
added. VPT2 and DBOC are not CCSD(T)-F12b implementations.

## Parsing, geometry, and compatibility

- Molpro→Sella now reads `show[1,f28.16],kb_geom_energy`, immediately after
  the undisplaced energy evaluation. The earlier eight-decimal `SETTING`
  value quantized energy changes and could spoil trust-radius decisions near
  convergence. The native constant-format check supplied by the user passed.
  The parser rejects conflicting/missing precision for new force calls and
  continues to read archived older outputs.
- F12 geometries use numerical gradients of the selected F12b energy.
  F12 harmonic tasks validate F12b energy evidence even when Molpro's
  FREQUENCIES heading only says UCCSD(T). Scalar-relativistic readers validate
  the explicit Hamiltonian. CFOUR selects the named HF/CCSD DBOC summary and
  checks the requested reference. Thus parser coverage was extended where
  required; the native validation during r7 is not claimed to have already run.
- ANL1 QZ geometry depends on accepted TZ geometry; QZ-dependent energies
  retain strict geometry hashes. Resource edits do not define a new
  scientific calculation. F12 components cannot be relabeled as a completed
  conventional ANL fallback by dropping only the F12 reference arm.
- All generic method, parser, scheduler and scratch code is independent of
  account/molecule/cluster. The sulfur examples and their site JSON are
  intentionally example-specific. Native WFU seeds are content-checked.

## Fresh r7: what the command does

`examples/anl/tools/fresh_run.json` is the editable campaign configuration.
Each prepared case also has a complete `input.json` containing the scientific
settings and live-policy defaults. The default campaign configuration:

- enables fully augmented ANL1-F12 v2 and vdW detection threshold 0;
- retains the example's 400 kcal/mol L2 search ceiling (10 kcal/mol extra
  at L1), temperatures through 5000 K, and the ordinary downstream-well
  threshold propagation. This conservative screening value is not a universal
  proof of negligible kinetic flux;
- permits five concurrent well searches per PES, sharing one QC node pool;
- selects six healthy partition nodes to exclude, preferring the named
  historically problematic nodes if healthy; excludes down/drained/unknown
  nodes and fixes the initial pool cap to healthy-count minus six;
- uses exclusive one-node QC allocations, a seven-day minimum clipped to
  the partition limit, one timeout retry with factor two, and two extra
  optimization-step chunks under the existing saved-geometry recovery;
- retains the reviewed CASPT2 active spaces for CH3S C–S (two doublet roots)
  and CH3SS C–S/S–S, MR supermolecule asymptotes, aDZ-F12 relaxation,
  L2 fixed-geometry orbital pivots, and aDZ/aTZ-F12 extrapolation of both
  CC reference energies and MR spin gaps;
- enables live controls without holding every scientifically rejected TS
  indefinitely. Unreviewed VRC channels have their own review hold, independent
  of this failed-search option. `hold_task_failures: true` keeps failed ANL/VRC
  dispatcher tasks available for manual fallback/release while independent
  tasks continue; `hold_failures: false` lets ordinary rejected searches finish.

The helper prepares inputs, checks reviewed WFU hashes, renders native inputs,
and checks the launch environment before cancelling anything. On `--apply`
it stops only r6's two marked PES trees, verifies their allocations/processes
have gone, renames r6 to a timestamped archive, performs the authorized cache
cleanup below, and launches r7 from a pinned source snapshot.
Integrals/scratch are not copied and no r6 result is imported.
Code pulls therefore do not change running r7 controllers.

The handoff journal survives preparation, retirement, archive-move, cleanup and partial
launch interruptions. Repeating the same command never launches a case with
an existing launch receipt again. An uncertain launch without a receipt is
reported and stopped for inspection rather than guessed safe. `launched`
means processes were started; confirm their subsequent health in the monitor.

### Cleanup between r6 and r7

The campaign configuration explicitly enables `cleanup.empty_home_wfu` and
`cleanup.remove_verified_home_scratch`. The first **deletes all contents of
the executing account's `~/wfu`**, as authorized for this handoff. The second
removes only owned `kinbot-*` directories immediately below
`~/.kinbot_scratch` that archived runtime records identify on this host.
Unknown directories, other-host paths, explicit retention requests and known
publication failures are retained and listed in `cache_cleanup.json` in r7.

Cleanup requires an empty **account-wide** Slurm queue, including pending and
completing allocations, and no detected local QC/controller or open cache
users. If unrelated work is present, the handoff stops before deleting caches
or launching r7; it does not cancel that work. Rerun the same command after it
finishes. Scheduler/proc inspection failures block deletion. Do not start
separate jobs during this short handoff: the guard is not a Slurm reservation.

Each target is moved to a journaled sibling before deletion, so a repeated
handoff never empties a newly created replacement directory. Symlink roots,
mounted directories and overlaps with run/source data are refused. Internal
symlinks are removed without following them. Reviewed r7 WFU hashes are checked
before and after cleanup. Archived native results and the repository's orbital
reviews stay intact. Package caches, virtual environments, `~/.cache`, conda,
mamba, basis-check reports, and node-local `/scratch/<account>` on compute nodes
are not swept by this login-host operation. There is no filesystem/account
name hardcoded in the cleanup implementation.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
export PATH="$PWD/.venv/bin:$PATH"
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
unset KINBOT_MOLPRO_SCRATCH KINBOT_MOLPRO_SCRATCH_FALLBACKS
.venv/bin/python -u examples/anl/tools/fresh_run.py \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7" \
  --replace "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage" \
  --apply
BASH
```

Omit `--apply` to prepare and preflight without changing r6. The same command
with `--apply` then uses that plan. Inspect `fresh_run.json` under r7 for the
selected node exclusions and exact cap. This is a campaign admission policy,
not a Slurm hardware reservation: KinBot leaves the six excluded nodes alone,
but cannot guarantee that other cluster users leave them idle. Node health is
not an automatic permission to remove later learned exclusions.

```bash
watch -n 15 "$HOME/KinBot/.venv/bin/python -m kinbot.run_monitor $HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
```

## Exploration, review, and continuation

Each well completes its L1/L2 exploration and publishes accepted L2 evidence.
The outer PES waits for **all discovered wells**, applies its L2 pathway
selection (including reverse discoveries), and only then schedules ANL for
selected states plus required CBH references. VRC can run reviewed channels
after their own well finishes. Unknown active spaces produce a review request
without blocking exploration or the later ANL network. See
[the ordering and selection audit](OUTER_PES_ANL_SELECTION_2026-10-21.md).

With `orbital_review_reserve: 3`, a current `needs_review` request reduces the
shared campaign ceiling by three. Running allocations finish; owned pending
jobs are held as needed. Accepted review models restore the original ceiling
and release only holds made by this policy. This frees campaign capacity for
recon; it does not bind three named nodes or cancel jobs to make room instantly.

At review, place the reviewed per-channel model and chosen theory/basis in
the case input's `vrc_tst_surface_models`, retaining existing reviewed entries.
Resubmit that same input, including while other wells are still exploring:

```bash
"$HOME/KinBot/.venv/bin/python" -m kinbot.run_control resubmit \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/ch3ss/input.json" --apply
```

This queues a revision to the existing owner, quiesces branch planning, adds
the new model, and continues without a duplicate driver. It does not permit
arbitrary changes to the science of running tasks. Resource rules, capacity,
exclusions, pending walltimes, selected retries, and network/rejected-channel
revisits use the existing `kinbot.run_control` commands. Records and prior QC
are preserved; dependent derived networks are regenerated after selection.
Searching/newly discovered wells inherit the accepted model at their ready
boundary. Existing reviewed channels remain unchanged. Changing the vdW cutoff
through file resubmission still requires the global postprocessing phase.

Manual numerical fallbacks now accept an embedded JSON proposal with native
input edits or generated-worker keyword overrides, executable commands,
allocation-only setup/environment, optimizer settings and per-attempt resources.
Use one exact target or repeat targets for an explicit batch. Requests validate
before cancellation, archive old attempts, retain sibling work, and remain
fenced to their observed generation. These controls can also patch an unsent
staged task; they cannot silently rewrite completed dependent results.
See [the complete live-control commands and fallback schema](LIVE_DRIVER_CONTROLS.md).

## Validation and limits

Tests cover all eight ANL variants for closed/open shells, basis/method/parser
agreement, QZ/5Z correction assembly, strict same-basis DKH/NR pairing,
high-precision energies, VRC extrapolation/orbital staging, shared capacity
drain/release, review reserves, and interrupted/partial fresh handoffs.
The full repository suite, including cleanup and live controls, passed with
**2,603 passed, 7 skipped, 429 subtests passed** (`.venv/bin/python -m pytest -q`,
106.72 seconds). The one warning is the existing test import of deprecated
`imp`. Skipped native/environment-dependent tests are not native QC evidence.

The audited launch sequence is:

1. Prepare complete inputs and reviewed seeds; validate code and environment.
2. Stop r6 controllers/allocations, verify termination, and archive r6.
3. Verify account inactivity, clean the explicitly enabled home caches, and
   reverify seed hashes. Stop here if cleanup is blocked.
4. Launch one pinned PES controller per case using one shared admission pool.
5. Explore discovered wells at L1/L2, with configured parallel well slots.
   Reviewed VRC channels can progress after their own well is ready.
6. Select the outer L2 PES after all discovered wells finish; only then plan
   ANL on selected states and required CBH references. TZ geometry feeds QZ;
   dependent single points require the correct accepted geometry.
7. Hold new VRC channels for reviewed active spaces, keep other work moving,
   and accept the existing live input-refresh procedure for continuation.
8. Assemble thermochemistry and correction potentials, run ROTD sampling,
   and run MESS when their required evidence is complete.

The existing SH symmetry recovery, VPT2 symmetry recovery, node-local scratch
verification, geometry-step and timeout recovery, outer-PES selection, and
live refresh tests remain in the suite. Native ANL v2 calculations are deferred
to r7 as authorized. Missing native functionality/bases or scientifically
invalid stationary points remain reported failures; they are not accepted by
loosening identity or curvature checks.

The last supplied r6 geometry evidence showed active CPU/I/O for all four
long jobs. One minimum had 29/36 curvature probes; the other three were still
far from their force targets. That established activity, not eventual success.
The precision correction addresses a real interface defect, and the F12
prescription changes r7's work, but neither supplies a reliable completion ETA.

Discovery is not proof of reproducing every stationary point in the supplied
paper diagram. VdW threshold zero admits only complexes with a positive
computed depth; it does not manufacture one or force the paper's connectivity.
Multiple distinct saddles can connect the same wells. The outer kinetic
selection and exact-identity reuse avoid unnecessary ANL work while retaining
raw evidence. Native IRC/curvature evidence and comparison with the paper are
still needed to establish completeness and quantitative accuracy.

## Additional audit: ordering, recovery, and resource limits

The executable paths audited are `pes_phase`/outer selection, `kb` branch
dispatch, `input_refresh`/`search_control`, shared `run_policy`/`pool_drain`,
dispatcher/ordinary retry transactions, and the fresh-run/cache-cleanup journal.
The review found and fixed these additional live-control gaps:

- Additive orbital reviews were unnecessarily restricted to global
  postprocessing; they now reach ready VRC branches during exploration and
  propagate to later wells at their safe boundary.
- A well becoming ready between proposal and consumption is included before
  publication. Competing revisions block the proposal without partially
  exposing it to workers. External input paths preserve their original base.
- Manual fallbacks are embedded, archived and generation-fenced. Successive
  attempts use the effective worker and retain explicitly selected resources.
  Native setup runs only in the allocation, not the login-host preflight.
- Approved task repairs preserve their logical VRC identity across later
  review refreshes. Unapproved subsequent manifest edits fail verification.
- Task-only failure holds keep expensive branch drivers available without
  freezing ordinary rejected TS guesses and preventing PES completion.

| Condition | Automated handling and acceptance boundary | Live operator route |
|---|---|---|
| Gaussian QC-SCF failure | Bounded alternative SCF algorithm; retain the requested convergence threshold | Supply reviewed `kwargs.scf`, guess/grid/symmetry keywords for one failed worker |
| Diagnosed SH/diatomic density oscillation | Fresh guess and automatic symmetry recovery, with force-frame handling; never accept energy-only apparent convergence | Targeted retry; native SH comparison supplied by the user supports this recovery |
| Gaussian frequency/VPT2 framework symmetry failure | Existing bounded symmetry recovery retained; matching accepted geometry required | Task fallback/retry; no false success from an error-terminated log |
| Optimization step budget exhausted | Additional chunks retain the current optimizer/Hessian; bounded later allocation recovery can seed saved geometry | Change steps/optimizer parameters and retry; a new optimizer does not promise the old Hessian |
| Sella geometry-update ODE failure | Existing bounded Cartesian recovery from verified geometry | Explicit Sella keyword fallback |
| Incorrect minimum/TS order or disconnected geometry | Curvature/connectivity rejection remains enforced | Retry a reviewed numerical strategy; never waive stationary-point validation |
| I/O, scratch or permission failures | Classified bounded recovery and verified local scratch; never home fallback | Repair setup/environment and retry the affected allocation |
| Slurm time limit | Confirmed timeout, archived attempt and bounded larger allocation within partition limits | New walltime/resources; native partial restart only where verified support exists |
| Native CC/MR failure | Native diagnostic and convergence checks; no generic success or unlimited retry | Exact native keyword/executable fallback, with unchanged method/state/geometry |
| Uncertain submission or completion racing cancellation | Retain ownership, reconcile receipts, do not submit a guessed duplicate | Inspect the blocked/uncertain receipt before proceeding |
| Missing VRC active space | Hold that channel for orbital review, proceed with approved work, reduce shared capacity for recon | Add a reviewed model and resubmit the same input |

The pool admits ready work from both PES trees under one serialized ceiling.
Reducing it drains running work and holds excess owned pending allocations;
raising it permits more work on the next poll. It does not create ready tasks
when scientific dependencies block them. Molpro uses requested MPI ranks with
one OpenMP thread per rank; threaded backends receive matching CPU allocations.
CPU and native memory changes are coupled, including Molpro per-rank overhead.
Controllers use bounded numerical-library threading. Exclusive nodes and
`--mem=0`/explicit memory are separate Slurm settings, as documented in the
[Slurm sbatch manual](https://slurm.schedmd.com/sbatch.html).

Utilization is dependency-limited, not guaranteed maximal: each relaxed scan
uses the previous geometry/orbitals, and leaves within a VRC worker are still
visited serially. Independent wells, ANL graphs and the two PES cases can run
concurrently. Numerical QZ/F12 gradients and curvature probes can dominate
elapsed time. Raising CPUs without a native scaling measurement is not assumed
to improve them. Seven-day requests favor completion but can reduce backfill
opportunities; the live walltime controls permit deliberate adjustment.

The manual/literature recheck preserves the intentionally revised ANL policy.
Molpro's [F12 manual](https://www.molpro.net/manual/doku.php?id=explicitly_correlated_methods)
supports restricted-reference UCCSD(T)-F12 and configurable auxiliary sets but
excludes DKH/X2C. Its [relativity manual](https://www.molpro.net/manual/doku.php?id=relativistic_corrections)
requires overriding Hamiltonian settings after a DK basis declaration.
The [CFOUR DBOC manual](https://cfour.uni-mainz.de/cfour/index.php?n=Main.CalculationOfDiagonalBorn-OppenheimerCorrections)
supports CCSD with RHF/UHF but ROHF only at HF, hence the authorized component
exception. An [MRCC developer clarification](https://mrcc.hu/index.php/forum/running-mrcc/404-ccmaxit-keyword-not-properly-applied)
states that `ccmaxit` applies to the CCSD program, while the higher-order MRCC
program can use a compiled iteration limit. A manual keyword retry therefore
cannot be advertised as curing every higher-order iteration exhaustion; a
compatible rebuilt executable may be needed and can be selected per attempt.

The [Veillet et al. paper](https://pubs.rsc.org/en/content/articlepdf/2026/cp/d4cp04873c)
uses F12/CBS energies with smaller-basis F12 geometries and identifies
multireference/submerged-barrier and vdW effects in CH2S+SH. It supports treating
those pathways carefully, not imposing its network as an automatic acceptance
criterion. The chosen ANL v2 and reviewed VRC models are the user's modified
prescription; native convergence, active-space adequacy and PES completeness
remain scientific checks during r7.
