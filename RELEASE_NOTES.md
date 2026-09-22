# Unreleased

## Stereochemistry and symmetry

**Keep configured stereoisomers separate throughout a calculation.**
KinBot keeps `chemid` as its connectivity identifier. A stereochemical
identifier distinguishes configured tetrahedral centers and double bonds.
Ordinary species keep their existing calculation names. Configured species use
different names for QC jobs, saved results, reused objects, and direct/PES output.
The same comparison checks IRC endpoints and fixed configurations during
conformer searches.

Install `kinbot[stereo]` to get RDKit and use the full set of features.
Unsupported assignments use the available legacy treatment with a warning.
This fallback does not establish the missing stereochemistry.

**Select the optical population.**
`optical_population` accepts `specified` (the default) or `racemic`.
The first option keeps the specified configuration. The second includes its
whole-molecule mirror. It does not request all possible diastereomers.
Both reaction endpoints constrain the allowed TS optical population.
`stereo_reference` keeps the requested configuration in generated PES inputs.
`stereo_legacy_inputs` identifies saved inputs that permit reuse of old job names.
KinBot checks the saved structure before it reuses those names.

**Keep different stereochemical reaction pathways.**
Virtual substitution distinguishes homotopic, enantiotopic, and diastereotopic
sites in the related reaction motifs. It does not change the atoms sent to QC.
The checks include heavy atoms and reactions with several reacting sites.
Forming and breaking bonds distinguish reacting atoms from spectator atoms.
Global atom-equivalence groups remain unchanged.
Direct MESS and PES output keep different stereochemical pathways, also
when `lowestpath` is selected. Repeated observations within one pathway still
use the existing lowest-barrier selection. A MESS Union adds the different
pathway contributions without requiring MC-TST.

**Use the symmetry and properties of each conformer.**
Each selected conformer keeps its geometry, electronic energy, ZPE, frequencies,
Hessian when available, calculation source, and status together.
MC population filtering uses each conformer's rotational and optical weights.
MESS then uses those same weights. Comparisons allow equivalent atom
permutations and proper rotations, so the same conformer is not counted twice.
Explicit mirror pairs receive no additional optical multiplier.
Conflicting numerical observations produce warnings. They do not require the
whole calculation to stop. MC-TST continues to disable HIR scans.

**Keep the graph-based rotational (external) symmetry rules with limited corrections.**
Fixed stereochemistry prevents exchanges between incompatible configurations.
Pyramidal XY3 has a threefold rotational contribution. Planar XY3 has sixfold.

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

**Record HIR evidence and recover incompatible scans.**
Shared records contain units, calculation sources, raw and projected
frequencies, rotor axes, symmetry numbers, scan angles, energies, and point
statuses.
If saved scan definitions no longer match the selected structure, KinBot first
tries scan recovery. If recovery fails, it removes HIR and restores the full
harmonic frequencies with a warning. It keeps the old calculation files.
The initial L1 result is loaded completely when no later optimization replaces
it. Reused product calculations keep the accepted structure and its properties.

**Separate disconnected well networks.**
KinBot follows reactions between bound wells from each requested reactant.
Separated-product channels are exits. A shared product does not connect two
otherwise disconnected well networks. Each requested disconnected network gets
its own MESS input and output. For example, this prevents specified R-butanol from acquiring
an unrelated S-butanol network through a common product.

**Complete MESS jobs before reporting their result.**
Local, Slurm, and PBS execution waits for the requested calculations to finish.
KinBot checks the exit result and newly written output. Solver failure is
reported separately from successful reaction generation. Inputs without an
accepted reaction network remain available for inspection, but are not run.

**Narrow bug fixes and reaction examples.**
The initial species reconstruction keeps its checked frequencies. Input-only
MESS writing converts collision-energy units in the same way as normal writing.
Three-fragment product names agree with each other through PES assembly.
Long configured names use bounded file names for MESS intermediates.
Four peroxy H-transfer inputs compare diastereotopic pathways with an achiral
control, using HIR or MC-RRHO. Their collision parameters are illustrative.
Tests include saved methanol, peroxy, pyramidal, and loose-TS structures.

The new stereochemical rate counting targets MESS. MESMER receives the same
species names but does not implement the new pathway and optical sums.

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
