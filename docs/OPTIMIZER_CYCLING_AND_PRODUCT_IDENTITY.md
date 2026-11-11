# Optimizer cycling and product identity

## Coupled external-optimizer evaluations

`run_optimizer` asks the KinBot Gaussian calculator for a coupled energy and gradient evaluation. Sella's energy-first, forces-second API calls then use one native `Force` result at a geometry. Previously those calls could launch separate SCF calculations and supply an energy and gradient from different evaluations. The coupled mode applies during external optimization and restores the calculator's ordinary mode on return or exception.

## Repeated Cartesian minimum steps

The Sella Cartesian minimum observer uses its cached evaluated points. Twelve consecutive steps at the trust-radius floor must show a repeating two-to-four-step gradient pattern, at least three near-opposite gradient reversals, at least three negative model-agreement ratios, and less than ten percent improvement in the minimum force norm. This identifies persistent oscillation while allowing ordinary slow progress.

For that pattern, the floor is reduced by a factor of ten, at most twice. The finite-difference Hessian step remains independently configured. The requested force threshold and subsequent frequency and identity checks continue to determine acceptance. Persistent cycling raises a dedicated diagnostic after the bounded reductions. The `.optimization_budget.json.cycle.json` record retains evaluated energies, positions, forces, model-agreement ratios and decisions across continuation. Saddle, internal-coordinate, cell and explicitly constrained optimizations are outside this recovery's scope.

The captured reverse-complex trace triggers the observer at step 82. The export extends to step 452. This replay establishes that the detector recognizes the observed pattern. Local simulated-calculator tests establish coupled evaluations, ordinary Sella convergence, bounded decisions and continuation behavior. A new native Gaussian optimization has not been performed, so convergence of this particular complex remains unverified. The available native output does not establish the electronic origin of the alternating gradients.

## Colliding legacy product names

A legacy numerical ChemID can represent distinct complete molecular graphs. PES assembly now catalogs the full saved product identities before combining reaction lists. Colliding bare names receive deterministic identity suffixes in the logical network. Original QC filenames and databases remain the source artifacts.

Energy selection checks the observed geometry, charge, spin and complete identity before comparison across workers. A direct reader presented with conflicting full identities raises an error. Cached MESS product and endpoint metadata receive the same contextual aliases as the reaction graph. Aliased L3 results require corresponding qualified filenames, while the existing complete-L2 fallback remains available when L3 is incomplete.

Two captured CSSH radical graphs provide the regression example. H-C-S-S and S-C-S-H previously shared a numerical name, allowing the lower energy of one graph to be combined with the other's geometry. The corrected assembly separates their H2 product channels. This is a graph-identity correction with no molecule-specific branching in installed code.

## Validation and deployment

Focused tests cover optimizer entry points, error propagation, continuation, live-work preservation, controller refresh, Gaussian recovery, product identity, native-file preservation, and cached MESS generation. Both captured PESs were replayed through the normal postprocessing entry point with QC submission disabled. Native MESS validation and the explicitly provisional scouting model are recorded separately in the local pass deliverables.

Updating `composite` installs these procedures for subsequent work. The update itself schedules no calculation. Existing native attempts and provisional-model choices retain their own provenance.
