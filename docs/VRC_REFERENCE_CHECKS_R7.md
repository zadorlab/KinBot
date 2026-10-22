# r7 native CC/MR reference checks, 2026-10-22

This panel runs new native calculations before the proposed CH + H2S refresh.
It does not accept a model, change production inputs, replace a controller,
cancel work, or clear any scan-review hold. Native results are still pending.

## Native export repair, October 22 morning

The first panel returned nine native failures at MATROP and five unsubmitted
candidate points. These were diagnostic-input failures, not failed chemical
state comparisons. MULTI `NATORB,SPIN=3` correctly selects a quartet, but the
general density selector interprets `SPIN=3` as S=3 and requests MS2=6.
The export now explicitly uses `MS2=3`. See the separate conventions in the
[MULTI manual](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi)
and [density selection manual](https://www.molpro.com/manual/doku.php?id=general_program_structure).
Old prepared decks receive paired diagnostic-card corrections in
their executed `.native.inp`; the execution receipt records the repaired
selectors and native-input hash. Prepared inputs retain their original hashes.

The next native attempt confirmed MS2=3 but again failed to find the density.
The first repair was incomplete: per-root NATORB cards do not save each root's
AO density to its own orbital record. In state-averaged MULTI, explicit `DM`
is needed to save the individual densities, and they share the record selected
by `ORBITAL`. The corrected diagnostic block uses `ORBITAL,3000.2; DM` and
loads every root from 3000.2 with explicit STATE/MS2 selectors. The normal
SA canonical/natural orbital records remain 2141.2/2140.2. Both prior generated
input versions are repaired in the executed copy, preserving original inputs
and archiving failed attempts. This is a diagnostic export change only.

Before panel fan-out, `run` now executes one existing singleton production
probe under a one-job cap. All other panel manifests count against capacity
but receive no submissions. Full native energies, all requested density
exports, hashes, MO orthonormality, electron traces and ROHF occupations must
validate before the remaining panel starts. A preflight failure stops there.
Successful execution does not automatically approve the chemical state match.
This native preflight is still pending on Blodgett; local tests use fixtures.

The four missing CH3SS points were skipped because the independent VRC branch
had no ANL `network_request.json` yet. Staging can now use the hash-bound
separated fragment geometry/maps and the current explicit electronic-state
models. It uses production's own perception and state-model functions, requires
exact fragment-map agreement and refuses undeclared states or charged
fragmentation. It does not infer a quartet merely from electron parity.

For the already-staged panel, run these before restarting its controller:

```bash
python examples/anl/tools/requested_orbital_recon.py retry-reference-exports PANEL
python examples/anl/tools/requested_orbital_recon.py production-checks PRODUCTION PANEL --extend
```

Both commands require an idle panel and hold its campaign/manifest locks.
The first archives only known legacy density-export failures via the existing
retry machinery. The second fills missing interior/asymptote roles, preserving
all existing point selections and attempts. Review the returned skipped list;
a missing point never counts as verified. Resume using the usual `run PANEL`
command. No production refresh is required for this repair. The corrected
native export still needs confirmation on Blodgett; local tests do not run Molpro.

## Production scope audit after the export failure

Audited at `c56e150c`: the failing `MATROP LOAD ... DEN` diagnostic path is
called only by `molpro_reconnaissance`. Neither the production single-point
generator (`kinbot/molpro.py`) nor the ASE/Sella gradient generator
(`kinbot/ase_modules/calculators/molpro.py`) loads these individual densities.
The missing DM directive therefore does not invalidate their existing energies
or gradients. No production controller replacement or QC archival is required
for either diagnostic export defect.

The production generators specify MULTI/RS2 wavefunctions with `WF` spin=2S,
use the selected root explicitly (or the declared lowest-energy selector for
energy-only reference calculations), read the SA canonical orbitals from
2141.2 for CASPT2 and propagate the packed SA natural orbitals from 2140.2.
The separate ROHF-UCCSD(T)-F12b reference uses its explicit WF multiplicity;
the pivot/cube path selects the alpha orbital type for UKS, without a density
SPIN selector. The dispatch SPIN placeholder is also multiplicity minus one.

Validation: 167 existing tests passed for Molpro ASE, VRC corrections,
multistate routing, orbital pivots and scan dispatch. A separate renderer audit
checked 112 inputs using the repository's six current scalar-leaf models plus
the proposed CH + H2S model, in both aDZ-F12 and aTZ-F12, for initial and packed
continuation records: target-doublet gradients, target/reference MR energies,
and quartet CC energies. Synthetic geometry/checkpoint bytes were used solely
to exercise generation; no native calculation or live production write was
performed. This establishes the input/selector scope, not native state identity
or absence of unrelated scientific errors. The panel remains the native check
at the actual hash-bound production geometries; existing scan holds remain.

## What runs

* CH + H2S: six fixed-geometry jobs, at 2, 3.5 and 30 A in aDZ-F12 and aTZ-F12.
  Each uses the proposed CAS(5,4), four doublets plus one quartet, all five
  root-resolved CASPT2-F12 energies, and a quartet ROHF-UCCSD(T)-F12b energy.
  Each basis starts from its verified 2 A CASCI checkpoint, exchanges source
  8 with 11 once, and retains 12–14. These test the proposed initialization;
  they are not a replacement for a relaxed scan.
* Existing scalar production leaves: two fixed-geometry jobs per leaf, in
  the **actual saved scan basis** (aDZ-F12 in r7). Use the accepted bonding-region
  point nearest 3.5 A and the separate optimized 30 A fragment reference.
  If the scan has not reached 3.5 A, use the nearest available accepted point.
  If it has not yet recorded its final far grid point, use the hash-bound
  fragment geometry/checkpoint already saved in its scan context.
* All six current leaves give twelve additional jobs: CH3S -> CH3 + S D1/D2,
  CH2SH -> CH2 + SH D1, C + SH3 D1, CH3SS C–S D1 and S–S D1. Together with the
  new-channel six, expect **18 jobs, at most three simultaneous allocations**.
  Discovery is from the actual run, so any additional eligible leaf is listed
  explicitly; absent/unverifiable leaves are reported as skipped, not passed.

Every production probe retains its saved geometry, atom order, CAS, state
average, selected roots, correlation settings and final packed natural
checkpoint. The original source-orbital permutation is not applied again.
Fresh ROHF, CASSCF and root-resolved CASPT2 calculations are diagnostic
single points; no gradients or geometry optimizations are requested. Native
reference-density saves and orbital exports are added without changing the
state-average objective. Convergence controls match the saved production
profile. Full MR root energies are recorded even where production selects
the lowest quartet, so the report can identify which root that selector used.
These are not byte-for-byte replays of numerical-gradient inputs.

The staging preflight verifies native geometry and artifact hashes, verifies
that the current correction contract agrees with the saved scan, and checks
the checkpoint before use. An inconsistent point is reported as skipped.
Rejected geometries are never substituted. If any expected leaf is skipped,
its identity remains unchecked; inspect the staging report.

The controller shares one cap across every manifest, honors the current
reserved slots and node exclusions, and continues independent probes if one
fails. Staging and launch use one campaign lock. Append the production checks
before starting the campaign. Repeating the append command is idempotent;
repeating `run` resumes the same manifests. Keep this panel outside r7.

## What identity is actually tested

MULTI `DM` saves **all individual CASSCF reference-root densities** to the
explicit diagnostic `ORBITAL,3000.2` record. Immediately after the
independent ROHF calculation, before CC, MATROP exports the AO overlap, SA-CAS
orbitals, ROHF orbitals/density and the state-specific MR densities. Native
matrix files are hashed in `execution.json`. The commands use Molpro's
[individual density saving and selection](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi)
and [MATROP matrix export](https://www.molpro.net/manual/doku.php?id=matrix_operations).

The read-only report checks MO orthonormality, electron traces and the ROHF
occupation pattern, then computes:

1. Squared singular values of the overlap between the ROHF singly occupied
   subspace and the MR active subspace. These are invariant to rotations/signs
   within either subspace; a value below 0.90 is flagged for inspection.
2. Each MR root's density in the ROHF MO basis, its occupations on the ROHF
   singly occupied orbitals, and its Frobenius distance from the ROHF density.
   It flags when the production-selected MR reference is not the closest
   density. With a `lowest` selector, the correlated energies determine the
   root actually selected for this comparison.
3. Native CI coefficients, active AO characters/occupations, reference ROHF
   orbitals and printed CC T1/D1/amplitude diagnostics, alongside all labelled
   energies. XML orbital exports and the WFU remain available for inspection.

These test **reference correspondence**, not a many-electron CC/MR overlap.
The MR density here is the CASSCF reference density underlying CASPT2; it is
not a relaxed CASPT2 density. Degenerate roots can rotate and produce an
ambiguous nearest-root assignment. The 90% screen is a diagnostic threshold,
not a scientific acceptance rule. No pass automatically approves a model,
proves CC accuracy, or establishes continuity at untested geometries. Missing
matrices are reported as unavailable. Converged energy alone is insufficient.

## Energy correction audit

At each common geometry and separately for each basis, the construction is

`E_best,target = E_CC,reference + E_MR,target - E_MR,reference`.

Each energy term is extrapolated using aDZ/aTZ-F12 and power 3.7. The relative
curve subtracts the **whole assembled expression at the same 30 A reference**:

`V_target(R) = E_best,target(R) - E_best,target(30 A)`.

The frozen sampling potential is referenced separately; the ROTD correction
is the trusted relative curve minus the frozen sampling relative curve.
The three energy terms at any point must use the same geometry. No spin
multiplicity factor is inserted into this arithmetic. Existing channels'
assembly and finite-reference subtraction were already correct; regression
checks exercise both the original triplet/doublet fragments and the new
singlet/doublet fragments, including a nonzero far target–reference gap.

The bug was in **reference eligibility/routing**: it restricted the auxiliary
reference to the ground-fragment coupled spins and could silently choose a
direct curve. The fix honors an explicit compatible reference multiplicity,
preserves inherited CC/gap settings, and continues to validate the reactive
target against the accepted products. It does not invent or approve an
excited reference automatically.

For CH + H2S, a quartet cannot be a ground-fragment reactive coupling, but
CH(a4Sigma-) + H2S is a possible auxiliary reference. CH ground/excited terms
are tabulated by [NIST](https://webbook.nist.gov/cgi/cbook.cgi?ID=C3315375&Mask=1000).
Its nonzero asymptotic excitation is removed consistently by the MR gap and
full reference subtraction. Higher spin alone is not proof of better CC
accuracy: the panel checks the determinant/MR reference correspondence and
provides CC diagnostics before that choice is accepted.

## Active spaces, averaging and flips

There were four CH3S scalar paths before the new CH + H2S proposal, from
three chemical channels. The counts below are spatial roots within each
spin, not spin projections, rate degeneracies or additional reactive flux.

| Channel / reactive leaves | CAS and chemical content | SA doublet / quartet roots | Evidence and remaining concern |
|---|---|---|---|
| CH3S -> CH3 + S, D1 and D2 | (5,4): methyl SOMO plus all three S 3p orbitals | 3 / 3 | Source 11–14; far S occupations about 4/3 each and methyl about 1. S(3P) motivates three spatial components per spin. D1/D2 represent the bound CH3S 2E components; finite distortions may split them. |
| CH2SH -> CH2 + SH, D1 | (5,4): CH2 frontier pair, 2e, plus SH pi pair, 3e | 4 / 2 | Source 11–14; far SH 1.5/1.5, CH2 about 1.206/0.794. Two roots per spin support SH(2Pi); two extra doublets retain the screened low states. Compact CASSCF keeps the intended characters; expanded (7,6) leaves extra S–H sigma/sigma* about 1.98/0.02. |
| C + SH3, D1 | (3,4): C 2p set, 2e, plus SH3 SOMO, 1e | 3 / 3 | Source 12–15; far SH3 about 1.003 and three C p occupations about 0.681 each. C 2s source11 about 1.956 stays inactive. C(3P) motivates three components. The actual scan has H-transfer failures and remains held regardless of this panel. |
| Proposed CH + H2S, D1 | (5,4): CH frontier set, 3e/3 orbitals, plus S donor, 2e/1 orbital | 4 / 1 | Near source8 and12–14; far donor is source7. Four doublets retain the screen; one quartet tests the auxiliary Sigma reference. New mixed-spin optimization is not validated by the earlier doublet-only screen. |
| CH3SS -> CH3 + S2, D1 | (7,5): methyl SOMO plus S2 pi/pi* set, 6e/4 orbitals | 1 / 1 | S2 ground Sigma has one spatial component. Already packed active18–22; do not repeat the original17↔19 source swap. |
| CH3SS -> CH3S + S, D1 | (7,5): CH3S sulfur frontier pair, 3e, plus atomic S 3p set, 4e | 6 / 6 | Two CH3S orbital components times three atomic S components motivate six per spin. Bound CH3SS has one selected parent component; averaging six roots does not add six rate paths. |

The original-channel evidence and validation limits are recorded in
[VRC_CHANNEL_HANDOFF_R7.md](VRC_CHANNEL_HANDOFF_R7.md#ground-terms-surfaces-and-averaging-are-different-quantities).
The two provisional CH2SHS cleavage models still await actual accepted review
requests; this panel does not promote them or fabricate production geometries.

Each current average gives equal **total** weight to the two spin blocks.
For example 4D+2Q gives each doublet 1/8 and each quartet 1/4; 4D+1Q gives
1/8 and 1/2. These are orbital-optimization choices, not thermal weights.
Extra roots reduce single-state bias and retain nearby components, but do
not certify a converged active space or prevent all crossings.

The CH + H2S 8-to-7 change is an ordering change between independent source
calculations: identify the donor by character, select near orbital8 once,
repack it into11, and propagate the whole selected span11–14. Applying far
index8 would select the wrong sulfur orbital. Source15 has S–H antibonding
occupation about0.086 at2 A, so omitting the S–H pair remains a compact-model
approximation that needs monitoring.

For other channels, checkpoint continuation also stabilizes the active span
but is not a quantitative root tracker. The initial three channel models
select the lowest correlated quartet; native quartet ordering changes were
already observed. Taking the lowest energy is **not** proof of character
continuity. CH2SH and C+SH3 currently use Q1. Current reviewed models disable
active-orbital alignment; no claim is made that temporal root/orbital flips
are impossible. The new invariant subspace/density checks supply additional
native evidence at selected points, leaving the full scan audit necessary.

## Commands

Run on Blodgett. The directory is intentionally new, so an existing campaign
cannot be overwritten. No production refresh belongs in this block.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
panel="$HOME/kinbot-runs/r7_vrc_reference_checks_20261022"
.venv/bin/python examples/anl/tools/requested_orbital_recon.py candidates \
  "$HOME/kinbot-runs/r7_ch_h2s_orbital_recon_20261022" \
  examples/anl/ch3s_anl1_f12/orbital_reviews/ch_h2s_cs_reference_check.json "$panel"
.venv/bin/python examples/anl/tools/requested_orbital_recon.py production-checks \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7" "$panel"
nohup env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  "$HOME/KinBot/.venv/bin/python" -u \
  "$HOME/KinBot/examples/anl/tools/requested_orbital_recon.py" run "$panel" \
  > "$panel/controller.log" 2>&1 < /dev/null &
echo "Panel controller PID: $!"
BASH
```

Monitor and report (the latter prints available completed evidence):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
"$HOME/KinBot/.venv/bin/python" "$HOME/KinBot/examples/anl/tools/requested_orbital_recon.py" \
  monitor "$HOME/kinbot-runs/r7_vrc_reference_checks_20261022"

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
"$HOME/KinBot/.venv/bin/python" "$HOME/KinBot/examples/anl/tools/requested_orbital_recon.py" \
  reference-report "$HOME/kinbot-runs/r7_vrc_reference_checks_20261022" \
  > "$HOME/kinbot-runs/r7_vrc_reference_checks_20261022/reference-report.txt"
cat "$HOME/kinbot-runs/r7_vrc_reference_checks_20261022/reference-report.txt"
```

The report independently reconstructs the new-channel CBS terms and anchored
interaction energies, with the direct MR curve alongside. Existing-channel
DZ checks are labelled as such; they cannot establish TZ/CBS or whole-scan
state continuity. Retain the entire panel directory: inputs, native outputs,
execution receipts, matrix exports, orbital XML files and checkpoints.

After results are reviewed, the existing r7 snapshot needs an in-place
controller refresh to load the eligibility fix, followed by publishing the
new channel through the existing live-review command. This is a code upgrade,
not a new PES run. Future runs started with this version can specify an
auxiliary reference through normal reviewed live input updates. Existing
channels' models, grids, sampling uncertainty and QC work remain unchanged.
