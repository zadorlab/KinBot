# Native CAS/MRCI evidence

Generated Molpro CASPT2, CASSCF/recon and MRCI inputs share native printing
controls, including the general VRC renderer, ASE/Sella evaluations and the
profiled RotdPy sampling adapter. The ANL external dispatcher also adds these
to custom MULTI/MRCI decks, respecting explicit custom print thresholds.

Recon already used `GPRINT,ORBITALS,ORBEN,CIVECTOR`. The shared version adds
`REF`, `GTHRESH,THRPRINT=0.05,PRINTCI=0.05` and `ORBPRINT,2` inside MULTI.
These print optimized occupied/active orbitals, two nearby virtuals per irrep,
labeled AO coefficients above magnitude 0.05 and CI/reference coefficients.
They do not change solver accuracy, the Hamiltonian, orbital/root selection,
records or the relaxed-reference Davidson energy reader. `PRINTCI` controls
printing; numerical `CIVEC`/`COEFF` tolerances remain unchanged.

Sources: [Molpro global printing and thresholds](https://www.molpro.net/manual/doku.php?id=program_control),
[MULTI print options](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#print_options)
and [MRCI print options](https://www.molpro.net/manual/doku.php?id=the_mrci_program).
No correlated-density calculation or extra MRCI natural-orbital transformation
is requested. CASSCF natural orbitals describe the reference space; input and
native record-selection lines identify the orbitals passed to correlation,
including canonical orbitals where requested.

In any inline state model, these optional general controls have these defaults:

```json
"orbital_diagnostics": {
  "orbital_coefficient_threshold": 0.05,
  "ci_coefficient_threshold": 0.05,
  "virtual_orbitals": 2
}
```

Thresholds must be finite, positive and less than one; the virtual count is a
nonnegative integer. No separate settings file is required. The examples retain
their existing inline orbital indices, checkpoints and electronic settings.

After each native invocation KinBot writes `<stem>.orbital_audit.json` beside
its `.out`, including failed native exits. Recon uses the executed `.native.inp`;
production uses the executed `.inp`. General external and ASE execution records
hash the report on success. Legacy VRC energy reads refresh the report too.
Missing output produces an incomplete-evidence report for an MR input;
single-reference CC output does not get a CAS report.

Each report includes:

- Exact input/output hashes, paths, byte size and native line references.
- Requested source indices, packing rotations, destination indices, OCC/CLOSED,
  electron/spin/root declarations, weights, records and correlation calls,
  read from executed input rather than assumed from a manifest.
- Every executed MULTI/correlation section, including numerical-gradient calls.
- Native active/closed partitions, spins, root counts, natural occupations,
  printed orbital energies and labeled AO coefficients. Closed and nearby
  virtual orbitals remain available for inspection at the active-space boundary.
- Native geometry/basis evidence, CI vectors, state energies, reference weights,
  selected-energy markers, Davidson variants and printed overlaps/records.
- Missing/unsupported tables, partition/spin/root-count mismatches,
  occupation-sum inconsistencies, native errors and incomplete termination.

This does not replace energy/convergence checks or scientific review.
`evidence_recorded` means extraction found the expected evidence without those
inconsistencies. It does **not** certify no swapping. `review_required` lists
specific missing/inconsistent evidence. Extraction errors are recorded as
`extraction_error`, not disguised as native convergence failures.

AO coefficients are thresholded amplitudes in a nonorthogonal basis, not atom
populations; omitted coefficients are not zero. Signs/components change under
phase, orientation and degenerate-subspace rotations. Occupations and energies
alone cannot establish state identity or detect every collapse. Retain WFUs and
orbital exports for quantitative overlap analysis. Equal energies alone are not
treated as collapse.

To rebuild reports from old files or inspect both PES trees:

```bash
/home/hram/KinBot/.venv/bin/python -m kinbot.molpro_diagnostics \
  /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261016_r3_dft
```

This recursively reads `.out` files, writes adjacent reports and prints a status
line per MR file. It does not submit/cancel/resume jobs or change native files.
Live outputs are snapshots; rerun after completion. Large sampling scans are
separate from the lightweight ten-second monitor; workers already write reports
automatically without a shared concurrent index.

Local validation extracts all 32 retained successful native CASPT2 outputs
(10 CH3S, 22 CH3SS), checking indices, occupations, AO coefficients and state
counts against executed inputs. It also covers repeated MULTI calls, missing
evidence, mismatches, failed exits, dispatcher decks, legacy rendering, RotdPy
and solver-retry interaction. Added printing is documentation-checked and
renderer-tested; its native production execution remains a Blodgett check.
