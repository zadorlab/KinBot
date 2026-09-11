# Unreleased

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
