# Downstream repairs after the r5 overnight run

## Native evidence and diagnosis

Run: `/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r5_dft`.
Evidence: user's 17 October 08:16 PDT monitor and subsequent compact native
reports. The saved-results restart moves obsolete attempts and controller logs
intact to an external `.repair_archive.*` sibling. These are software/interface
failures; production thermochemistry, rates and fluxes are not validated yet.

1. **CH3S discovery passed L1/L2 and accepted C-S homolysis.** Five reaction
   searches were accepted and three rejected/failed. One insertion search
   explicitly failed SCF convergence even after its bounded retry. The
   other two rejection reasons are not established by the monitor alone.
   Neither an empty queue nor an accepted channel certifies the full PES.
2. **VRC reached its Molpro relaxed scan in both PES.** At requested C-S
   separation 2.5 A, `relaxed_scan_step_0001.xyz` instead contained the seed's
   approximately 30 A geometry. CH3S output gave energy -437.413169043280 Eh;
   CH3SS output gave -835.235261489150 Eh. The Python coordinate check rightly
   refused to assign these gradients to the requested geometry. Native
   files are under each `vrctst/molpro_scan/*/tasks/relaxed_scan/` (and retained
   `attempts/relaxed_scan/1/`). **Do not reuse those energies as 2.5 A points.**
3. **CH3S CBH-ANL ran reference L2 and began composite tasks.** Reference SH,
   `state_20691a101e49c36b/composite/tasks/gaussian_vpt2/vpt2.log`, terminated
   normally at 00:09:17. Its route requested `Symmetry=(PG=C1)` but Gaussian
   retained nuclear framework `C*V`, with Abelian subgroup C2V. KinBot wrongly
   rejected that complete VPT2 result solely for not printing framework C1.
4. Reference H2S, `state_ed133ec78c51e793/composite/tasks/gaussian_vpt2/vpt2.log`,
   did fail natively at 00:07:57: Link 717 reported `New: CS - Old: C1`, despite
   the PG=C1 route. This is distinct from the SH parser rejection.
5. **CH3SS CBH reference L2 completed HSS** at 03:00:56, but model lookup
   failed. Requested `[S]S`, charge 0, doublet; exported `S=[SH]`, charge 0,
   doublet, symbols S/S/H. The accepted positions were
   `[1.0782044661,0.0110778729,-0.1098514366]`,
   `[3.0591896588,0.0966331206,-0.0405757860]`,
   `[3.4015152586,-0.9200311831,-0.8637874889]` A. These are alternative
   localized bond/radical representations of the same H-S-S connectivity.
   The source is `anl_automatic/reference_l2/state_6aa3d2b553427a20/`.

## Repairs and documentation basis

- Open the staged WFU **before** current geometry, basis, charge, spin and
  thresholds. FILE restores saved information, not just orbital coefficients.
  The generated ASE gradient, VRC correction single-point and RotdPy sampling
  decks now all follow title/memory/FILE/current-settings order. The method
  still selects the same record/indices and uses DIAB alignment. The force
  coordinate check remains strict. The reconnaissance deck already placed
  FILE before geometry; this repair does not invalidate its separate native
  active-space confirmation results.
  [Molpro quickstart, restarting calculations](https://www.molpro.net/manual/doku.php?id=quickstart),
  [Molpro restart control](https://www.molpro.net/manual/doku.php?id=program_control).
- Use documented `Symmetry=None` for newly generated **VPT2** tasks. PG=C1
  limits symmetry; it did not bypass the framework identification that failed
  in the H2S anharmonic analysis. None bypasses symmetry processing. The
  explicit framework-error recovery command is updated to the same option
  and archives its old input/native output. This does not change atomic-SCF
  handling, L1/L2 optimization routes, method, grid, or Conver10.
  [Gaussian Symmetry documentation, vendor mirror](https://www.conflex.co.jp/gaussian_support/symmetry.php).
- VPT2 parsing records the requested group cap and observed native framework
  labels separately. It no longer assumes they must match. Normal native
  termination, requested method/basis/dispersion, complete fundamental bands
  and internally consistent anharmonic ZPE remain required. The H2S error
  cannot pass that validation. Successful old SH VPT2 is reusable in principle.
- Reference model export retains requested SMILES after checking explicit-H
  connectivity, element/isotope labels, total charge and radical count.
  Spin is independently checked on import. Perceived SMILES remains in source
  provenance. Legacy exported records undergo the same restricted comparison.
  Hydrogen transfer/loss, different connectivity, charge, radical count and
  incompatible stereo are rejected. This does not merge ATcT identities or
  introduce a sulfur-only alias.
- VRC failures now include the task's native execution error rather than
  `unknown task failure`. The two-PES dashboard summarizes routine records
  and dependency waits and places both case overviews first. Full records
  and log tails remain available through `--all --details` and `--json`.

## Pruning decision

Keep the currently requested broad reaction search for this repair. The
1000 kcal/mol cutoff is applied **after** a successful TS search, before IRC;
it cannot prevent an unsuccessful TS search. It can allow extra expensive
IRC/product/recursive work. The CH3S root only attempted eight searches and
did reach both downstream branches, so these specific blockers are not caused
by that cutoff. Homolysis has a separate product-energy threshold derived
from it, so lowering it blindly could remove a desired channel. After native
barrier energies are available, an informed smaller energy window or selected
reaction-family run is possible through existing general JSON controls.

## Restart and acceptance

The user's later instruction selects **saved-results resume in `_r5_dft`**.
Run `examples/anl/tools/resume_on_blodgett.sh` after pulling `composite`.
The helper preflights the current inputs and requires both runs to be inactive.
If any jobs/processes remain, it stops the migration and leaves them untouched.
It archives logs, invalid VRC attempts, database snapshots and original metadata
outside the run, reparses successful SH VPT2 without submitting QC, restages
only failed/unstarted VPT2 with Symmetry=None, and preserves accepted L1/L2/L3
and composite results. It then launches both detached login-host controllers.
Unreviewed failures prevent launch. HSS's accepted reference model is reused
through the general identity check. The entire invalid VRC first-scan trees
are archived; they must be recomputed with the corrected geometry ordering.

The first resume stopped at `gaussian_vpt2: cannot migrate unavailable`.
That was a helper defect: optional native execution failures transition from
`failed` to `unavailable`, but the migration only handled `failed`. The exact
production status is now covered by a regression. Both statuses reach the
same native-error check; only the recognized Gaussian framework-group failure
is restaged. Preflight-unavailable components remain unavailable, and unrelated
native failures still require review. Native-success reparsing likewise handles
an old optional failure only after validating the native files and provenance;
derived skips are released only when their own failed sources have recovered.

Rerunning the same resume command also handles the interrupted migration:
updating both input JSONs during the first invocation no longer prevents the
second invocation from archiving stale derived data for the remaining case.
The interruption regression checks this explicitly and preserves accepted
native records. Another fresh PES restart is not required by this helper error.

Failed TS decisions and compact exit markers remain cached; deleting these
would incorrectly erase genuine reaction exclusions. Their old native attempts
are archived, and the monitor identifies the hash-matched markers as cached
rejections. A new attempt's changed marker is reported normally. Parser-history
records bound by successful-result hashes retain symlinks to their external
archive. Fresh controller logs begin on resume; scientific history is not
rewritten to pretend failed attempts succeeded. Reconstructed reaction decisions
may appear again in the new logs.

The migration updates electronic-level and SO postprocessing data in prepared
input JSONs and refreshes their preparation hashes. See the
[species spectroscopy review](SULFUR_ELECTRONIC_LEVELS_2026-10-17.md).
The approved additive convention now corrects S, SH and CH3S once in ANL
assembly. Current input corrections and MESSPF levels also apply to cached
references, without changing their QC source hashes, geometries, frequencies,
raw electronic energies or ZPE. Old derived profiles/formation records/MESS
outputs are archived and regenerated. Other numerical QC settings remain unchanged.

The helper does not delete shared scratch or unrelated jobs.
Each PES retains its 10-allocation cap and the
eleven configured excluded nodes. All selected active spaces/indices/seeds,
IPEA 0.25, shift 0.1, aDZ/aTZ exponent 3.7, tight convergence, DFT preset,
ordinary VRC and 1% flux-error target remain unchanged.

Local tests check rendered FILE ordering across all three production paths,
rejection of the wrong 30 A geometry, complete/incomplete VPT2 parsing,
framework-error retry archiving, HSS export/import, and rejection of different
reference structures, plus dashboard and scoped restart behavior. These are
not local executions of licensed Gaussian or Molpro. In particular the new
H2S Symmetry=None VPT2 route and corrected seed-to-2.5 A gradient still require
native confirmation after resume; passing Python tests is not an end-to-end guarantee.

The prior SO integration passed **419 regression tests and seven native
RotdPy/MESS/MESSPF checks**, zero failures/skips. After the resume-status
and full CBH electronic-state coverage changes, **428 regression tests passed**
in 24.71 s, with zero failures/skips. One Sella/JAX cache
permission warning occurred in a passing test. Shell syntax and whitespace
checks passed. [Commands, counts, source hashes and limits](../examples/anl/tools/downstream_repair_validation.json)
are retained. These targeted checks cover the changed paths; they are not a
claim that all workflow bugs are eliminated.
