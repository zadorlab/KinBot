# CH3S: orbital review and ordinary KinBot production

The accepted minimum choice is **C–S (5,4), natural orbitals 11–14**. Native aDZ/aTZ CASPT2-F12
confirmation and the matched IPEA 0.25/SHIFT 0.1 check are complete.
Production uses CASPT2 only; MRCI is historical evidence.

Use the [consolidated production handoff](../caspt2_production/README.md) for
preparation, dependency/native-backend checks, scoped old-run retirement,
two-PES launch, monitoring and acceptance. The [parameter reference](../caspt2_production/PARAMETERS.md)
explains all supplied controls and their master/branch provenance.
The six portable checkpoints are supplied in that package.

The inputs request ordinary single-bottleneck VRC, aDZ-F12 orientation sampling,
aTZ-F12 constrained relaxed geometries, and a high-spin ROHF-UCCSD(T)-F12b/aTZ
anchor plus CASPT2 aDZ/aTZ spin gaps extrapolated with exponent 3.7. The
independent CBH/ANL1-F12 workflow supplies stationary-point thermochemistry.
There is no inner/outer stage or automatically added VRC spin-orbit potential.

The example files are data for the [general orbital-review procedure](../../../docs/ORBITAL_REVIEW_WORKFLOW.md).
The historical `orbital_reconnaissance.json` is an exploratory aTZ recipe, not
the final production input or a new task to rerun. The archived earlier README
is retained in `docs/archive/sulfur_readmes_20261016/`.

Full KinBot/RotdPy/MESS execution remains to be performed and reviewed on
Blodgett. Preparation tests and accepted electronic points do not establish
completed rates or a factor-of-two kinetic uncertainty.
