# Reusing accepted TZ/QZ minima for ANL anharmonic corrections

Automatic Gaussian B2PLYP-D3(BJ) L2/VPT2 profiles accept `cc-pVTZ`,
`cc-pVQZ`, `aug-cc-pVTZ`, and `aug-cc-pVQZ`. For an ANL1 or ANL1-F12
minimum, the anharmonic extension uses the matching TZ/QZ family. An
augmented ANL electronic-component setting does not independently change
the L2/VPT2 basis family.

## L1 to L2 handoff

The owning optimizer freezes the selected, normally completed L1 minimum
before L2 refinement. Reuse requires a supported Gaussian B2PLYP-D3(BJ)
profile, verified profile generation, charge/multiplicity, ordered atoms,
default masses, finite properties, and the exact selected geometry. The
actual rendered or reused L2 worker must start from that L1 point. A later
HIR geometry change, an unrelated existing worker, an unsupported physical
model, or an unverified profile does not establish this lineage.

Once the existing optimizer acceptance checks accept L2, both basis points
are retained in `anl_basis_reuse`. Each includes geometry, theory, state,
database/job/row provenance, harmonic results and a content hash. The L1
point also identifies the L2 seed worker. Publication retains this metadata
in the accepted database row. Restart, verified cross-worker atom mapping,
and reference-model export preserve it. No additional geometry, energy or
frequency rejection criterion is imposed on the chemistry search.

## What is reused and what still runs

- A verified aTZ L1 minimum followed by an accepted aQZ L2 refinement supplies
  both basis-specific geometries. ANL stages the two VPT2 calculations at
  these respective minima without repeating their optimizations or adding
  separate harmonic jobs. Gaussian VPT2 still computes the derivatives it
  needs; the stored harmonic spectrum is not a substitute for that work.
- With only an accepted aQZ minimum, only the missing aTZ minimum/harmonic/
  VPT2 branch is added. With only aTZ, the missing branch is aQZ. An
  unaugmented minimum cannot stand in for an augmented calculation.
- `gaussian_vpt2` remains the accepted **L2-basis** task, even when L2 is QZ.
  This preserves a correctly labeled single-basis fallback if the extension
  is explicitly omitted. The legacy live option `anharmonic_qz=false` means
  omit the additional basis branch; for QZ L2, that additional branch is TZ.
- Existing task graphs and submitted work are preserved. New graphs use
  the new provenance; this change does not retrofit unknown old lineage or
  cancel, restart or replace native calculations.
- B2PLYP energies and harmonic derivatives do not replace the separately
  specified coupled-cluster/F12 ANL components or higher-level Hessians.

## Extrapolation and scientific scope

At each basis's own accepted minimum, define the native anharmonic ZPE
difference as `d = ZPE(VPT2) - ZPE(harmonic)`. The existing declared
inverse-power protocol is retained:

```
d(CBS) = d(QZ) + [3^p / (4^p - 3^p)] * [d(QZ) - d(TZ)],  p = 3.7
```

The actual basis pair and exponent are recorded in the assembled component.
This is the existing extrapolation approximation, not a claim that its error
vanishes with augmentation. Only the anharmonic difference is extrapolated;
the assembly uses one separately specified harmonic ZPE.

Thermal fundamentals retain the existing normal-mode treatment: each VPT2
calculation supplies its own native harmonic vectors and shifts at its own
geometry. Verified atom permutations and proper rotations put those vectors
in the accepted species frame; the existing mass-weighted overlap and rotor
projection machinery transfers the shifts. It does not match different
bases merely by sorting frequencies. The upper/lower task roles are explicit
when the L2 reference is QZ. Rotor contributions are included once on the
accepted L2 rotor model. Existing native VPT2 quality and mode-ambiguity
checks remain active.

## Validation

`tests/test_anl_basis_reuse.py` exercises accepted database records through
ANL planning, staging, native-format result parsing and CBS assembly, both
basis families, restart publication, reference export, reordered atoms,
missing/incompatible evidence, no-conformer-search seeds, mismatched
existing workers, repeated staging and a partially completed live graph.
`tests/test_anl_anharmonic_pair.py` also covers a QZ L2 reference with
augmented TZ/QZ mode-resolved shifts.

These are local simulated-site and native-format fixture tests. They do not
constitute a newly completed Gaussian calculation or native production-run
validation.
