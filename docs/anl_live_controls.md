# Live ANL protocols, cross-PES sharing, and TZ/QZ anharmonicity

These controls apply to ordinary automatic ANL graphs. They do not depend on
particular molecules, production directories, or reference answers. Scientific
changes are applied by the owning driver. The CLI previews the selected graphs
and affected task attempts before it publishes requests.

## Cross-PES sharing

`python -m kinbot.anl.control share ROOT --source PES_A --target PES_B`
previews a directed sharing route. Add `--apply` to enable it. Paths are absolute
or relative to ROOT, and both PES owners must be within that live-control root.
The route persists for future graphs and does not require a controller restart.

The receiving driver can adopt a canonical graph before the receiver starts
native work. A running canonical graph can be joined; a completed one can be
consumed. Its original owner, Slurm IDs, native files, and attempt history stay
in place. No completed flags or arbitrary individual energy components are
copied between graphs. A receiver with native work already started is retained.

Eligibility checks include the accepted molecule and electronic state, L2
theory, geometry, energies, vibrational evidence, connectivity/TS endpoints,
requested protocol, and actual graph-to-source binding. The cross-PES geometry
match uses a 1e-4 angstrom tolerance with atom permutation and proper rotation;
it is not a match based solely on a species name or SMILES. The underlying
accepted-state comparison also checks its energy and frequency tolerances.
Different conformers, states, theories, or selected protocols are not merged.
Routes must be acyclic and aliases remain flat.

The receipt is `ANL_BASE/cross_pes_anl_reuse.json`; downstream consumption is
recorded in `pes_reuse_consumed.json`. Identity is rechecked when used. A changed
scientific source blocks consumption for review instead of silently accepting
new results. Resource-only changes do not change scientific identity. A live
protocol edit on a receiver detaches that receiver and archives its old alias;
the canonical allocation continues. Editing a canonical protocol with registered
receivers is held until those receivers have been explicitly detached/revised.
This is conservative whole-graph reuse, not a general geometry-transfer cache.

## Live calculation choices

```bash
python -m kinbot.anl.control protocol ROOT \
  --graph 'CASE/WORKER/anl_automatic/states/STATE/composite' \
  --geometry tz --harmonic tz --reason 'Retain accepted TZ geometry and Hessian'
```

Add `--apply` to queue that exact selection. Repeated `--graph` patterns select
subsets, `--species` selects exact accepted SMILES (including their declared
states), and `--all` selects all currently existing automatic graphs under ROOT.
The preview lists the resolved graphs, retained geometry/harmonic sources,
omitted tasks, and affected job IDs. It does not change native calculations.
`--all` is a selection of current graphs, not a rule to simplify future species.

| Choice | Effect |
|---|---|
| `--geometry tz --harmonic tz` | Retain the declared TZ CC geometry/Hessian; omit QZ geometry/Hessian; evaluate retained higher-level electronic terms at TZ geometry. |
| `--geometry qz --harmonic cbs` | Use QZ geometry and the declared TZ/QZ harmonic-ZPE pair. MESS retains QZ harmonics. |
| `--geometry l2 --harmonic l2` | Use accepted L2 geometry/harmonics for the retained electronic protocol. |
| `--omit-task TASK_ID` | Omit a specific task and any consumers that cannot run without it. The preview identifies the closure. |
| `--omit-component dboc` | Omit the DBOC contribution. |
| `--omit-component higher_order_tz` | Omit the paired ANL1 TZ higher-order increment; the DZ increment remains independently eligible. |
| `--omit-component higher_order_dz` | Omit the paired ANL1 DZ higher-order increment; the TZ increment remains independently eligible. |
| `--omit-component higher_order` | Omit all higher-order increments. |
| `--omit-component core_valence` / `scalar_relativistic` | Omit the named complete correction family. |
| `--tier ANL0-F12` | Rebuild the selected graph at that tier, retaining equivalent accepted tasks and their original IDs. |
| `--tier ANL1-F12`, `ANL1`, `ANL0`, `L3`, `L3-SP@L2`, `L2` | Other supported tier selections. |
| `--no-anharmonic-qz` | Disable the QZ double-hybrid branch and use the accepted TZ VPT2 treatment where applicable. |
| `--anharmonic-qz` | Enable that branch for an eligible ANL1/ANL1-F12 graph. |
| `--clear-omissions` | Clear task and component omission lists. |
| `--geometry default --harmonic default` | Restore the selected tier's normal sources. |

Options compose with the graph's previously applied choices. Omission lists
replace the corresponding earlier list; repeating the option supplies multiple
entries. A tier change does not implicitly clear explicit geometry or omission
choices. Include the desired resets when changing them.

Removing one member of an energy difference or CBS pair never produces half a
correction. Its complete contribution is omitted or the higher-tier recipe is
unavailable. Records label the actual geometry, harmonic treatment, and omitted
families. CBH reactions and barrier differences still require a common compatible
expression; changing one species can therefore change the common expression
available to that reaction. The code does not silently mix incompatible species
protocols into a nominally uniform ANL result.

Changing a geometry invalidates its dependent single points. Calculations at
the old geometry remain archived evidence; they are not relabeled as calculations
at the new geometry. Equivalent tasks, unrelated completed results, and unrelated
running allocations are preserved. Existing site launch commands and minimum
recovery settings are retained during tier changes where their correspondence
is unambiguous.

The default is to wait for affected allocations. `--cancel-running --apply`
explicitly authorizes the owner to cancel only the previewed affected attempts,
after verifying scheduler ownership, stage files, and attempt generations.
Cancellation sends one job-specific `scancel` and waits for a terminal Slurm
state; `COMPLETING` does not release the graph for replacement. If completion
wins the race, the superseded result is archived with its attempt rather than
being relabeled at the new geometry. Retained allocations continue to count
against the shared pool. Protocol changes do not rewrite the pool's capacity,
exclusions, walltime policy, fairness settings, or revision; subsequent work
uses the same admission controller. Freed slots are available to all cooperating
branches, not reserved exclusively for the edited graph.
Ambiguous submissions block a change. The manifest transaction is journalled
and restartable; superseded task directories precede replacement manifests in
the archive. Old workers without manifest-lock support are allowed to finish
before incompatible manifest edits. Newly enabled branches receive site and
scheduler preflight before submission.

Already handed-off graphs are reopened through the existing live PES revision
owner, which archives derived assembly before rebuilding it. The protocol fence
is published first, so no new handoff can race it. This requires a live
revision-aware PES owner; a completely exited standalone calculation must first
be resumed under its owning workflow. Ordinary live choices do not replace a
controller. `python -m kinbot.anl.control wake-owners ROOT` retries a failed owner
wakeup while retaining its queued fence.

Inspect receipts with `python -m kinbot.anl.control status ROOT`. Full plans are
in `.kinbot_controls/anl_protocol_requests/*.json`; affected attempts and the
prior manifest are under each graph's `attempts/live_*/`. `protocol_preflight.json`
records site validation; a queued/applied request is not native QC completion.

## Shared ANL1/ANL1-F12 anharmonic extension

For polyatomic minima with accepted L2 B2PLYP-D3BJ/cc-pVTZ data, both ANL1 and
ANL1-F12 now add the following L3-owned branch:

1. `l3_vpt2_qz_geometry`: B2PLYP-D3BJ/cc-pVQZ optimization seeded from accepted
   L2 TZ coordinates, with force criterion 1e-5 eV/angstrom and Gaussian
   SCF convergence at least 1e-10.
2. `l3_vpt2_qz_harmonic`: separate harmonic stationary-point validation at that
   QZ geometry.
3. `gaussian_vpt2_qz`: VPT2 at that accepted QZ minimum.

Open-shell species use the Gaussian open-shell double-hybrid treatment. The
original accepted L2 geometry, Hessian, VPT2, and rotor scans are unchanged.
The QZ branch is independent of the expensive CC geometry branch and does not
wait for its completion. TSs and atoms do not receive a minimum VPT2 branch;
custom L2 methods/bases are not silently converted into this pair. Existing open
eligible graphs acquire the extension through their owner with no reset of
existing tasks. A normal resumed owner reopens its graphs before applying the
extension. A closed graph still in a completed live workflow can instead receive
an explicit `--anharmonic-qz` protocol request and owner revision.

The scalar expression is

```
ZPE = ZPE_harmonic(selected CC treatment) + delta_ZPE_anharmonic(TZ,QZ)
X_CBS = X_QZ + alpha * (X_QZ - X_TZ)
alpha = 3^3.7 / (4^3.7 - 3^3.7)
```

The anharmonic quantity being extrapolated is the native **total anharmonic ZPE
minus harmonic ZPE**, including the terms in that native result. It is not the
total double-hybrid ZPE and is not half the sum of fundamental shifts. CC
harmonic ZPVE is separately extrapolated when both recipe-matched TZ/QZ Hessians
are accepted and selected. Selecting a single retained harmonic result disables
that harmonic extrapolation. Disabling QZ VPT2 disables its extrapolation.
Native hashes, geometry/method/basis identity, and VPT2 quality checks are required
for both members. A failed or unaccepted upper result is not a completed CBS pair.

The exponent 3.7 is the existing profiled ANL extrapolation convention, now also
used for the requested double-hybrid anharmonic extension. The resulting
`VPT2-CBS` label identifies a profiled extension; this is not a claim that this
exponent or this double-hybrid treatment is separately benchmarked for every
species. Existing bounded native-failure recovery and explicit quality-review
requirements remain applicable.

## MESS frequencies, geometry, and rotors

MESS receives the selected accepted high-level geometry and harmonic spectrum:
normally QZ, or the explicitly retained TZ pair. Individual CC harmonic
frequencies are **not** CBS-extrapolated. The total harmonic ZPVE used in the
ANL energy remains separately extrapolated as described above. MESS geometry,
frequencies, and zero-point energy are independent inputs; they need not be
derived by summing the same list of effective fundamentals.

For each retained CC mode i, the TZ and QZ double-hybrid modes are independently
aligned and assigned using mass-weighted normal-mode overlaps. Then

```
delta_i(TZ) = nu_i_anharmonic(TZ) - nu_i_harmonic(TZ)
delta_i(QZ) = nu_i_anharmonic(QZ) - nu_i_harmonic(QZ)
nu_i(MESS) = nu_i(CC, retained) + delta_i(QZ)
             + alpha * (delta_i(QZ) - delta_i(TZ))
```

The two basis calculations use their own optimized geometries and mode vectors.
Matching is not by mode number or nearest frequency alone. The default checks
require squared overlap >= 0.5, harmonic mismatch <= 100 cm-1, and no competing
assignment within 0.1 overlap with a materially different (>1 cm-1) shift.
Wrong dimensions, bad geometry/state provenance, ambiguous correspondence, and
nonpositive corrected stable frequencies require review. Source maps, overlaps,
both shift sets, and the extrapolated shifts are recorded for later inspection.
An explicitly retained L2 harmonic spectrum can also receive these independently
projected paired shifts when its native vectors and frequencies verify.

Accepted L2 hindered-rotor identities are projected out of the retained Hessian
for MESS. Rotor potentials and their geometry remain the L2 scan data, including
the L2 geometry in the rotor block. If MESS demotes a rotor, mode projection and
both sets of anharmonic shifts are recalculated from the stored sources.

CBH-ANL 0 K energies continue to use the full harmonic ZPVE plus the full scalar
VPT2 correction. They do not simply delete the harmonic torsional ZPVE. A
consistent quantum hindered-rotor replacement would require its ground-state
energy and corresponding treatment of torsion-vibration coupling; the rotor
ground energy cannot generally be assumed zero. Retaining the vibrational
treatment is the explicit approximation here, while MESS uses the L2 rotor
partition functions for thermal/statistical mechanics. Large-amplitude modes
remain a limitation of perturbative anharmonicity.

The existing vdW binding path remains harmonic-only for both the complex and
its fragments in the binding-energy difference. Independently computed fragment
formation thermochemistry can retain VPT2. This change does not introduce
covalent CBH fragmentation of a vdW complex, nor VPT2 corrections into that
binding difference. Existing vdW/MESS mode and inner/outer path handling remains
in force.

The choice to retain QZ harmonics is consistent with the QZ-F12 rovibrational
treatment in [Glarborg et al., 2025](https://doi.org/10.1016/j.proci.2025.105942).
The [HFREQ2014 benchmark](https://doi.org/10.1021/ct500174q) studies CBS frequency
accuracy using pointwise force-field energies; it does not establish that this
same scalar ZPVE formula improves every directly extrapolated mode. The
[MESS manual](https://tcg.cse.anl.gov/papr/codes/mess/messmanual_3_23_16.pdf)
defines geometry/rotational constants, frequencies, and zero energy separately.

## Blodgett deployment and subsequent use

Run the deployment commands from the accompanying handoff in your SSH session.
The `git pull --ff-only` plus `merge-base --is-ancestor COMMIT HEAD` check selects
the tested code. The one-time `kinbot.run_refresh --apply --wait-idle 600` pins
that code and replaces the local PES/KB controllers while retaining their
recorded Slurm allocations and QC directories. Its idle wait concerns safe
staging/submission boundaries, not completion of long-running CC calculations.
It neither prunes reactions nor chooses lower protocols. The normal pool is 23.
Pending jobs are retained as well: do not add `--replace-pending` for this update.
Independent RotdPy campaigns retain their owners and allocations; these ANL
features do not require replacing their sampling workers.

After both PES explorations are accepted and the desired sharing scope is
reviewed, an example opt-in route is:

```bash
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
"$HOME/KinBot/.venv/bin/python" -B -m kinbot.anl.control share "$KB_RUN" \
  --source ch3s --target ch3ss
# Repeat with --apply when selecting this route.
```

Examples of later protocol previews (replace the graph placeholder with a
resolved path from the preview/inventory):

```bash
"$HOME/KinBot/.venv/bin/python" -B -m kinbot.anl.control protocol "$KB_RUN" \
  --graph 'ch3ss/WORKER/anl_automatic/states/STATE/composite' \
  --geometry tz --harmonic tz --reason 'Retain accepted TZ pair'

"$HOME/KinBot/.venv/bin/python" -B -m kinbot.anl.control protocol "$KB_RUN" \
  --species '[SH]' --omit-task ccsdtqp_dz --reason 'Omit this higher-order increment'

"$HOME/KinBot/.venv/bin/python" -B -m kinbot.anl.control protocol "$KB_RUN" \
  --graph 'ch3ss/*/anl_automatic/states/*/composite' \
  --tier ANL0-F12 --geometry default --harmonic default --clear-omissions \
  --reason 'Selected lower ANL protocol'
```

None of these simplifications or sharing routes is enabled merely by installing
the code. Keep the current production choices until the full L1/L2 networks are
available for the planned review. After an applied request, inspect the protocol
receipt and `kinbot.run_monitor ROOT`; implementation tests and a launched
controller receipt do not establish successful native results on Blodgett.

## Validation of this change

Local validation used native-output parser fixtures and simulated scheduler and
controller fixtures; no remote QC calculation was submitted by these tests.

* ANL, live refinement, controller refresh, RotdPy campaign and worker-code
  suites: **664 passed, 1 skipped**.
* Input refresh and CASPT2 handoff suites, in a separate process: **34 passed,
  1 deselected**. The deselected example renderer test already fails on the
  unchanged baseline.
* The full-suite audit reached the same 25 failures on both the working tree and
  a clean archive of baseline `ea154bb1`. These concern existing VRC example
  expectations, branch fixtures, and suite-order-dependent PES-owner fixtures.
  The repository-wide suite is therefore not reported as fully passing.

The exercised paths include native scalar correction extrapolation, independent
TZ/QZ mode permutations, rotor demotion, retained-TZ recipe assembly in both ANL1
variants, independent higher-order increments, sharing provenance rejection,
receiver detachment, ambiguous submissions, generation changes, scheduler-owned
cancellation, interrupted transaction recovery, closed-owner wakeup, old-worker
guards, site preflight, controller replacement, sampling adoption, and repeat
invocation. The full-capacity cancellation regression exercises 23 allocations,
both pending and running targets, the `COMPLETING` interval, cancellation versus
natural-completion races, admission-claim reconciliation, unchanged policy bytes,
and admission of exactly one replacement into the freed slot. Successful native
TZ/QZ VPT2 results and a successful Blodgett refresh
remain to be verified from their actual receipts.
