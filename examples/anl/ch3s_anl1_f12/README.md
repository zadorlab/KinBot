# CH3S: orbital review and ordinary KinBot production

The requested channels are C–S dissociation. These are input-data examples
of the [general orbital-review workflow](../../../docs/ORBITAL_REVIEW_WORKFLOW.md),
not separate molecular implementations. Current native evidence and the status
of each active-space candidate are in the
[local CASPT2 validation record](../../../docs/SULFUR_ORCA_VALIDATION_2026-10-15.md).
The earlier launch history is preserved in
[the archived README](../../../docs/archive/sulfur_readmes_20261016/ch3s_README.md).

`ch3s.json` describes the intended CBH/ANL1-F12 → VRC/RotdPy → MESS run.
Its active-space entry is the local minimum candidate (5,4), pending native
Molpro confirmation. The referenced state/surface-model
files must be populated from reviewed, basis-specific native records before
production; their absence is intentional and must not be bypassed with guessed
orbital indices. Multiple bonds or bottleneck stages can use different models
through the existing surface-model parameters and `install-selection` command.

The requested levels are:

- CASPT2-F12/aug-cc-pVDZ-F12 for orientation sampling.
- CASPT2-F12/aug-cc-pVTZ-F12 for the configured constrained relaxed path.
- RHF/ROHF-UCCSD(T)-F12b/aug-cc-pVTZ-F12 plus the equal mean of CASPT2-F12 and
  MRCI-F12+Q spin gaps at that basis for the electronic correction.
- The independent CBH/ANL1-F12 recipe for stationary-point thermochemistry.

The proposed CASPT2 controls are explicitly IPEA 0.25 and real shift 0.3 Eh
for sampling, the relaxed scan and the CASPT2 spin-gap component. Local ORCA
selection comparisons used IPEA zero consistently. Parameter insensitivity has
not been established; the user chose to assess transfer with the proposed
controls in Molpro rather than repeat a long local sensitivity campaign.

KinBot constructs matching relaxed and frozen points and passes the resulting
correction potential to RotdPy. Geometry optimization follows the configured
scan method; it is not an optimization using a composite spin-splitting gradient.
The fixed-geometry Molpro confirmation only tests whether the local ORCA
orbital/state choices remain usable with native CASPT2-F12. It does not repeat
the production relaxation scans. Conventional ORCA CASPT2 in an F12 orbital
basis is not an explicitly correlated calculation.

`orbital_reconnaissance.json` and the historical finite comparison recipes use
TZ-F12. The latter now record their comparison level explicitly, independently
of the production sampling basis. DZ-F12 sampling needs its own reviewed native
orbital model; a TZ-F12 WFU must not be relabeled or reused across bases.

After native confirmation, install the reviewed selection in the applicable
bond/stage, retain the declared physical states and thermochemistry data, and
run the ordinary `kinbot INPUT.json` workflow. Completion means validated
correction scans, the requested Monte Carlo uncertainty, audited thermochemistry
and an executed MESS calculation. A passing preparation/test suite or a chosen
CAS alone does not establish those results.
