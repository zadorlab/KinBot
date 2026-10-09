# Unreleased

Each KinBot and PES run starts its log with the KinBot and RDKit versions.

## Job results

A finished job's `done` stamp is no longer accepted before its result row is
in `kinbot.db`. The job writes its `.pkl` before the stamp, and over a network
file system the driver can see the stamp first. The last row then still
belonged to the previous step of the same job, and the next step (for example
the final TS optimization of a reaction search) started from that earlier
geometry. KinBot now records how many rows a job has at submission and keeps
waiting until a newer row exists; the existing 60 s grace period applies.

## Stereochemistry and symmetry

**Keep configured stereoisomers separate throughout a calculation.**
KinBot keeps `chemid` as its connectivity identifier. A stereochemical
identifier distinguishes configured tetrahedral centers and double bonds.
Ordinary species keep their existing calculation names. Configured species use
different names for QC jobs, saved results, reused objects, and direct/PES output.
The same comparison checks IRC endpoints and fixed configurations during
conformer searches.

Configured calculation names use `-s` followed by 16 hexadecimal characters.
The full stereoisomer identifier remains in saved records and is checked before
reuse. Names do not depend on `optical_population`: a racemate must still remain
separate from its diastereomers. The shorter names use a new calculation format;
calculations with the earlier 64-character suffix require new directories.
If an input SMILES leaves stereochemistry unspecified, the log reports the
assignment in the generated geometry at INFO level. This note also covers
partly specified SMILES. It is not emitted for supplied coordinates.

RDKit >= 2025.9.3 is required. KinBot selects its stereo-perception settings
explicitly and records these settings and the RDKit version in the log.
SMILES input still requires OpenBabel to generate the starting geometry.
An unsupported initial reactant stops before QC. An unsupported discovered
reaction is omitted with a warning that the network is incomplete. Its
calculation files remain available; unrelated reactions can continue.

The four stereochemical H-transfer example files are valid JSON templates.
Replace their quoted `SET_...` placeholders with local QC and queue settings
before running them. Replace the CPU count and job limit strings with integer
values. The HIR templates use matching placeholders for `method` and
`high_level_method`, and for `basis` and `high_level_basis`; keep these pairs
equal when replacing the placeholders for an L1-only HIR calculation.

Ordinary PAHs and biaryls are no longer rejected by their number of aromatic
rings. Scope checks use an aromatic graph copy so equivalent phenyl arms do
not appear different in one bond drawing. Aromatic radical resonance forms
are not mistaken for cumulenes. Virtual reaction-site labels do not create
a physical biaryl axis. The existing geometric optical calculation
still determines whether a harmonic, HIR or MC model includes its mirror.
No helical stereoisomer identifier is added: selection of one fixed helical
enantiomer and separation of helical diastereomers remain unsupported.

**Start older calculations in new directories.**
This release cannot resume calculations from older KinBot versions. It does
not convert old job names, directories, or result formats. Each new calculation
has a `.kinbot_run.json` format record. A restart requires that record and the
same RDKit version and settings. KinBot checks this before it changes old files.
Save an environment specification, such as a conda export or pinned Python
requirements, with each project so that this runtime can be restored.

### New feature: select one stereoisomer or a racemate

Treating a chiral well as one specified stereoisomer is new. In earlier KinBot
versions, a chiral well identified by the graph-based chirality rule was treated
as a racemate: the code assigned `nopt = 2` and wrote
`SymmetryFactor = sigma_ext / 2` in MESS. There was no option to restrict that
well to one enantiomer. The old rule detected atoms with four different bonded
groups; it did not detect all forms of geometric chirality.

The new input option has two values:

- `"optical_population": "specified"` is the **default**. It includes the
  stereoisomer assigned from the input geometry and excludes its opposite
  enantiomer.
- `"optical_population": "racemic"` includes that stereoisomer and its
  whole-molecule mirror. This restores the earlier counting of both enantiomers
  for those chiral wells. It does not include other diastereomers. For example,
  an RR input includes RR and SS, not RS and SR.

**Existing inputs change behavior when this option is omitted.** KinBot reads
the stereochemistry from the geometry; the input does not need an R/S label or
an explicit stereo option. Every run with the option omitted uses `specified`.
For a chiral well with external rotational symmetry number `sigma_ext = 1`,
represented by one structure, the previous `SymmetryFactor 0.5` becomes
`SymmetryFactor 1.0`. To retain the previous factor for this well, add
`"optical_population": "racemic"`. `SymmetryFactor` is a divisor: 0.5 includes
twice the statistical contribution of 1.0 for the same structure.

**Example: sec-butylperoxy, CH3CH2CH(OO•)CH3.** Its tetrahedral stereocentre
has H, methyl, ethyl, and peroxy groups. Use `charge: 0`, `mult: 2`,
`rotor_scan: 0`, and `multi_conf_tst: 0` to compare one harmonic structure.
For the supplied enantiomer, `sigma_ext = 1`, and MESS receives these factors:

| Treatment | Enantiomers included | `SymmetryFactor` |
| --- | --- | --- |
| Earlier KinBot versions | Both | 0.5 |
| New default, or explicit `specified` | The input enantiomer | 1.0 |
| Explicit `racemic` | Both | 0.5 |

This change does not imply a factor-of-two change in every rate coefficient:
the TS contribution also depends on which reflected reaction its reactants and
products permit. `stereo_reference` preserves the requested stereoisomer in
generated PES inputs.

**Conformational mirrors of achiral molecules are also counted.** For example,
a single gauche ethanol conformer has a distinct mirror although ethanol has
no fixed stereocentre. With `rotor_scan: 0` and `multi_conf_tst: 0`, earlier
versions wrote `SymmetryFactor 1.0`; the new treatment writes `0.5` when the
missing mirror is established. Both conformers belong to the same stereoisomer,
so this also applies in `specified` mode. A HIR model that already includes both
mirrors, or an MC model with both mirrors explicitly present, receives no extra
factor of two. Thus `racemic` restores the earlier choice of enantiomers, but it
does not undo other corrections to optical or rotational counting.

### Reaction pathways and statistical counting

**Find rings without generating every open atom path.**
Ring detection now enumerates closed cycles directly. It keeps the earlier
ordered ring list, including larger perimeters of fused rings. This reduces
the cost for PAHs without changing the resonance-structure search.

**Keep different stereochemical reaction pathways.**
Virtual substitution distinguishes homotopic, enantiotopic, and diastereotopic
sites in the related reaction motifs. It does not change the atoms sent to QC.
The checks include heavy atoms and reactions with several reacting sites.
Forming and breaking bonds distinguish reacting atoms from spectator atoms.
Global atom-equivalence groups remain unchanged.
When equivalent searches merge, KinBot prefers the earlier graph search's atom
selection if it passes the reaction-family checks. For example, propane keeps
`r12_insertion_R_7_2_3` instead of renaming it to `r12_insertion_R_7_2_1`.
Distinct stereochemical selections still receive separate searches.
Direct MESS and PES output keep different stereochemical pathways, also
when `lowestpath` is selected. Repeated observations within one pathway still
use the existing lowest-barrier selection. A MESS Union adds the different
pathway contributions without requiring MC-TST.

**Use the symmetry and properties of each conformer.**
Earlier MC-TST output could write symmetry-equivalent conformer copies as
separate RRHO blocks and overcount them. The old moments-of-inertia filter
was ineffective because `all()` returned a boolean instead of testing each
moment ratio. Duplicate checks now compare geometries under equivalent atom
permutations.
Each selected conformer keeps its geometry, electronic energy, ZPE, frequencies,
Hessian when available, calculation source, and status together.
MC population filtering uses each conformer's rotational and optical weights.
MESS then uses those same weights. Comparisons allow equivalent atom
permutations and proper rotations, so the same conformer is not counted twice.
Duplicate geometry comparisons use an E + ZPE window of 0.5 kcal/mol. This
allows more numerical variation than the earlier 0.2 kcal/mol window while
avoiding expensive comparisons between well-separated energies. Explicit
mirror comparisons retain their separate treatment without this energy gate.
Explicit mirror pairs receive no additional optical multiplier.
Conflicting numerical observations produce warnings. They do not require the
whole calculation to stop. MC-TST continues to disable HIR scans.
An L1 or L2 MC result with missing or non-finite energy or ZPE now fails the
affected optimization with a warning. Other reaction calculations continue.
The code does not substitute zero for these missing conformer properties.
Internal record-association errors still raise an error.

**Keep the graph-based rotational (external) symmetry rules with limited corrections.**
Fixed stereochemistry prevents exchanges between incompatible configurations.
Pyramidal XY3 has a threefold rotational contribution. Planar XY3 has sixfold.
The pyramidal rule and optical factor assume a harmonic single-well inversion
mode; a future explicit umbrella model that covers both minima must account
for their symmetry and mirror contribution without counting them twice.

A departing atom is distinguished from the spectator atoms when rotor (internal)
symmetry is calculated.

**Apply one optical-counting method to harmonic, HIR, and MC models.**
The method compares local molecular parts and accounts for the active rotor
bonds to assign optical symmetry factors. 
The local geometric RMSD tolerance to identify optical isomers is 0.1 angstrom. 
Existing HIR scans are used to establish mirror coverage.
If a HIR model already includes both mirrors, KinBot adds no second factor of
two. If an allowed mirror is clearly absent, its contribution is included.

For unresolved comparisons within the allowed optical population, KinBot uses
the selected calculation's Hessian to estimate the stable-mode harmonic energy
at an aligned mirror midpoint. Each MC conformer uses its own Hessian.
The optical factor is one for
an estimate at or below 4 kcal/mol, and two above that value.
This is an approximation. The estimate is not a calculated inversion barrier.
If the necessary data remain unavailable, KinBot uses factor one with a warning.
`optical_factor_assumptions` permits an explicit factor of one or two for a
named structure, with a necessary written reason.

**Keep usable rotor scans when another rotor fails.**
Shared records contain units, calculation sources, raw and projected
frequencies, rotor axes, symmetry numbers, scan angles, energies, and point
statuses.
Optimization replaces scans whose geometry or rotor definition no longer
matches the selected structure. It retains compatible scans, including usable
partial scans, and keeps the old calculation files. A failed rotor retains its
harmonic motion; the other rotors remain in the model. Projection starts from
the full Hessian of the selected calculation. If that Hessian cannot be used
safely, KinBot restores all harmonic frequencies and omits HIR with a warning.
MESS writing can check saved data and repeat this projection, but cannot submit
QC or wait for scans. Missing scan coordinates alone keep the usable rotor
potential and give a warned optical factor of one. A scan that changes to an
excluded stereoisomer or reaction pathway is unusable; an optical factor of one
cannot correct such a potential.
The existing rotor-zero energy check still disables all rotors if their common
optimized reference is inconsistent.
The initial L1 result is loaded completely when no later optimization replaces
it. Reused product calculations keep the accepted structure and its properties.

**Separate disconnected well networks.**
KinBot follows reactions between bound wells from each requested reactant.
Separated-product channels are exits. A shared product does not connect two
otherwise disconnected well networks. Each requested disconnected network gets
its own MESS input and output. For example, this prevents specified R-butanol from acquiring
an unrelated S-butanol network through a common product.

**New feature: choose whether MESS includes product complexes.**
`me_skip_vdW: 0` is the default. It includes accepted product complexes and
their barrierless exits in both single-well and PES MESS output.
`me_skip_vdW: 1` omits these wells and exits and connects each inner TS to the
separated products. Complex searches and optimizations still run. Inner-TS
Eckart depths use the selected complex energy with either setting. Both
writers use the existing PES approximation of one lowest-energy complex for
each set of product stereoisomers. MC-TST retains its selected-structure
tunneling reference. Omitting a complex removes its stabilization and explicit
capture/redissociation competition, so rates can change.
PES complex selection no longer depends on discovery order. Previously, the
first complex could keep its own energy reference even when a later complex
for the same products had a lower energy. All complex names for that product
set now resolve to the lowest-energy complex.
`correct_submerged` now uses connected bound-well ground energies only. A TS
below separated fragments is not raised to their energy.

**Complete MESS jobs before reporting their result.**
Local, Slurm, and PBS execution waits for the requested calculations to finish.
Earlier versions skipped MESS execution when `queuing: "local"` was selected.
With `run_me: 1`, a local run now starts MESS and waits for the result.
With `run_me: 0`, KinBot writes inputs and scripts without starting MESS.
Local execution, queue scripts, and manual scripts now use `mess_command`;
earlier versions ignored this option and used `mess`. The option specifies an
executable and optional arguments. For example, use
`"mess_command": "/opt/MESS/bin/mess"`. Quote a path or argument that contains
spaces. Shell variables and command substitutions are not expanded. Relative
paths are resolved from `me/`. Input-only writing needs no local MESS executable;
queued jobs find the executable in the compute node's environment. A local run
reports an error if it cannot start the configured command.
KinBot checks the exit result and newly written output. Solver failure is
reported separately from successful reaction generation. Inputs without an
accepted reaction network remain available for inspection, but are not run.
KinBot writes its reaction summary and PESViewer input before starting MESS.
For queued MESS jobs, it allows up to 60 seconds after the job leaves the queue
for the exit marker and new output to appear. A recorded solver failure still
fails immediately. Scheduler queries have a 30-second timeout and allow three
consecutive failures before reporting an error; a query failure does not mean
that the job has finished. Initial input-reference and saved-result conflicts
reported by `StereoRoutingError` are logged and exit with failure, without an
internal traceback. Other internal errors are not suppressed.

**Output formats and reaction-discovery cost.**
Summary files include a `# kinbot_stereopath {json}` comment before each
successful calculated saddle, including ordinary pathways. Homolytic-scission
placeholders can lack this comment. External summary readers must allow comment
lines; the PES reader uses this metadata to keep different pathways separate.
MESS microcanonical files use `<job-stem>.micro`, such as `me/mess_0000.micro`,
instead of the shared `micro.out`. Uncertainty samples and disconnected networks
have separate filenames. `me/mess_networks.json` lists the current MESS jobs.
For ordinary QC calculations with `high_level: 0`, the MESS header now reports
the L1 method and basis instead of unused L2 settings. This corrects the header
description; it does not change the calculation level.
The database includes one `stereochemistry/<job>` input-reference row per
well optimization job checked for reuse. It records the input geometry, full
stereoisomer identity, charge, and spin multiplicity. It is not another QC result,
and job-status checks do not repeatedly add it. Database readers must distinguish these
reference rows from calculation results.
Reaction discovery performs additional stereochemical graph comparisons.
A ten-molecule check with fixed coordinates measured about 2.2 times the
enumeration time of master. This timing excludes imports, characterization,
initial stereoisomer-reference assignment, and QC. It is not a measurement of
the total calculation time.

**Narrow bug fixes and reaction examples.**
H2 elimination now rejects longer ring walks when an equivalent shorter search
was removed as a duplicate. Distinct stereochemical paths remain separate in
PESViewer output, also when their barriers differ by less than 1 kcal/mol.
An exact stereoisomer entry takes priority over a connectivity entry in
`vrc_tst_scan` and `vrc_tst_noscan`, including an empty list. The same rule
applies before VRC calculations and when PES writes rotdPy input.
The initial species reconstruction keeps its checked frequencies. Input-only
MESS writing converts collision-energy units in the same way as normal writing.
Three-fragment product names agree with each other through PES assembly.
Long configured names use bounded file names for MESS intermediates.
Four peroxy H-transfer inputs compare diastereotopic pathways with an achiral
control, using HIR or MC-RRHO. Their collision parameters are illustrative.
Tests include saved methanol, peroxy, pyramidal, and loose-TS structures.

The new stereochemical rate counting targets MESS. MESMER receives the same
species names but does not implement the new pathway and optical sums.

**Wait for delayed calculation results.** After a submitted Slurm or PBS job
leaves the queue, KinBot allows up to 60 seconds for its complete result to
become available. This prevents an early failure or a duplicate submission
when result files appear late. A completed calculation error is still returned
without this extra wait.

**Native Q-Chem constraints are written correctly (#67).** The six Q-Chem job
templates imported ASE's stock `QChem` calculator, which knows nothing about
KinBot's `addsec` keyword and wrote it into `$rem` as `ADDSEC $OPT ...`,
which Q-Chem rejects ("Illegal rem input in read_rem"). KinBot's own
calculator, which writes the constraints as a separate `$opt ... $end`
block, has existed since 2023 and was already used by every Gaussian and
Sella template. The Q-Chem templates now use it too. This affected every
native Q-Chem transition-state search and hindered-rotor scan. Testing against a real
Q-Chem also turned up four further defects, fixed in the same change: bond and
angle constraints were computed but never written (only torsions reached the
`$opt` block); releasing a constraint indexed the geometry with one-based atom
numbers and raised `IndexError`; the `Total energy =` summary line printed by
Q-Chem 6.2 was not recognised, so converged jobs reported no energy; and an
empty `OPT` keyword made Q-Chem frequency jobs run another optimisation.
Validated on Q-Chem 5.4.2 and 6.2.1.

**PES post-processing writes `pesviewer.inp` first and no longer draws its own
graph (#82).** The pyvis-based "interactive graph" duplicated what PESViewer
does from `pesviewer.inp`, and it was the only post-processing step with a
third-party dependency that could fail or hang for reasons unrelated to the
kinetics; a PES whose KinBot runs had all finished could end without a
combined `pesviewer.inp`. The graph is removed and `pyvis` dropped from the
`plot` extra. `pesviewer.inp` is now written immediately after the energies
are assembled, before the rotdPy and MESS inputs, and each post-processing
stage announces itself in `pes.log`.

**Near-linear molecules keep 3N-5 vibrations, and MESS agrees.** `get_frequencies`
decided how many external rotations to project out by an absolute cutoff of 1e-5
on the mass-weighted rotation vectors. Any optimised geometry that is linear only
up to an optimiser residual (a bend of 0.01 degrees is enough) exceeded it, so
three rotations were projected and one component of the degenerate bend was
lost: 3N-6 frequencies instead of 3N-5. With Gaussian or Q-Chem this reached
MESS only through the rotor-projected frequency set; with FairChem or Sella,
where all frequencies come from this routine, it affected every near-linear
species (CO2, HCN, C2H2, ...). Linearity is now decided by three gates: all bond
angles within 2 degrees of 180; the curvature of the Hessian along the rotation
about the molecular axis is a vibration (above 50 cm-1), not a free rotation,
which distinguishes a linear molecule left slightly bent by the optimiser from
a genuinely bent minimum; and that curvature exceeds three times the residual of
the other rigid-body motions. Stored geometries are never modified. A species
judged linear is written to the MESS input as an exactly linear rotor, with a
comment, so that MESS's own moment-of-inertia test (I_min/I_mid < 1e-5) reaches
the same conclusion; otherwise MESS would treat a 179-degree CO2 as a
non-linear top with a tiny third moment, a factor of ~3 in its rotational
partition function. The reverse mismatch is handled too: a genuinely bent minimum so
close to linear that MESS's test would call it linear gets its three rotational
constants written explicitly, so MESS keeps the non-linear rotor that the 3N-6
frequencies assume. Both decisions are logged. The same routine now uses the
symmetric eigensolver, so degenerate modes can no longer come back as complex
numbers.

**Generated job scripts no longer depend on numpy reprs (#84).** Element
symbols and coordinates were pasted into the ASE/Molpro job templates via the
repr of numpy types, which under numpy >= 2 reads ``np.str_('C')`` and
``np.float64(1.06)``. Scripts from the five templates that do not import
numpy (Q-Chem TS search, IRC and HIR; NWChem constrained TS search; Molpro TS
search) failed to run. Values are now converted to builtin Python types
before formatting.

# 2.4.1

Bugfix release. Two groups of changes alter numerical output for inputs that
ran under 2.4.0 and earlier; they are listed first.

## Results change for unchanged inputs

**Hindered-rotor projection (#99).** The internal-rotation vectors that are
projected out of the Hessian were built from mass-weighted atomic
*coordinates* instead of mass-weighting the rotational *displacement*. For
any atom whose mass differs from the rotor-axis atoms this skews the vector,
so the projected frequencies handed to MESS (`reduced_freqs`) were wrong for
every `rotor_scan = 1` run since KinBot 2.0. On butane the corrected
projection reproduces the full spectrum minus the three torsions to within
9 cm-1, where the old one deviated by up to 590 cm-1. Because the error
largely cancels between a transition state and its reactant, the effect on
rate coefficients is typically 10-20% at combustion temperatures; quantities
without a cancelling partner (equilibrium constants, standalone entropies and
heat capacities, master-equation densities of states) see the full
single-species error, up to a factor of 1.5-2 for molecules with several
rotors at high temperature. The regression reference in `tests/frequencies.py`
now comes from a finite rigid rotation of the molecular fragment rather than
from the code itself.

**Multi-conformer TST (#100).**
- Direct (non-PES) runs never resolved the per-member energy offsets, which
  were written only as comments; every Union member sat at the parent energy.
  Offsets are now written as numbers on both routes.
- Member energies are referenced to the accepted parent structure instead of
  to the lowest member, so a conformer that becomes lowest only at L2 is
  placed at its true energy rather than displacing the whole Union.
- Each saddle member gets its own imaginary frequency and shifted Eckart
  depths; members whose L2 optimisation failed are excluded; a single
  surviving member is written as a one-member Union with its own properties.
- Population screening in `find_unique` passes the complete mode list to
  ASE's `IdealGasThermo` (API from ASE 3.23, already implied by the ASE 3.26
  requirement). Previously all imaginary modes were dropped before the call,
  which raises under current ASE for any transition state.
- Conformer ground energies are initialised from E + ZPE at L1, so L1-only
  MC-TST members are no longer equally weighted.

**Free rotors (#100).** The MESS free-rotor block now carries the rotor
`Symmetry` number, as the hindered-rotor block already did. A methyl free
rotor was previously overcounted threefold.

## Fixed

**Selected-structure consistency (#100).**
- With conformer search at L1, the rotor projection used the Hessian of the
  starting structure with the geometry of the selected conformer. The
  selected calculation's geometry, energy, ZPE, frequencies and Hessian are
  now loaded together.
- Every poll of the optimisation re-ran the conformer selection and reset the
  species geometry to the L1 conformer after the L2 geometry had been loaded.
  The selection now runs once.
- When no Hessian is stored for the selected structure (Gaussian conformer
  jobs, native Q-Chem L1), a frequency-only calculation at the fixed geometry
  recovers it. If recovery fails, KinBot warns and continues with harmonic
  frequencies and no hindered rotors.
- The accepted result is published under the conventional job name
  (`<name>_well`, `conf/<name>_low` or `<name>_well_high`), so restarts and
  PES post-processing read the same structure the optimisation accepted.
  Replaced output files are archived, not deleted.

**Hindered-rotor restarts (#100).** With `high_level = 0`, a lower point found
during a scan is now re-optimised at L1 before the scans restart, instead of
being taken as a stationary point directly from the constrained geometry.
Restarts archive the superseded scan results through the database rather
than deleting log files. `qc_opt_ts` honours its `high_level` argument.

**Miscellaneous.** Q-Chem L1 jobs print the Hessian; `formchk` is invoked
through `subprocess` with error checking instead of a silent `os.system`.

## Known issues

- Rotors whose far side has only multiply bonded neighbours (nitro, nitroso,
  isocyanate groups) are not detected, because dihedral detection requires a
  single-bonded terminal atom.
- The semi-empirical conformer pre-search is not yet restricted to wells.

## Contributors

Luka Dockx, Judit Zádor.
