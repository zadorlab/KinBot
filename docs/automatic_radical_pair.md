# Automatic native verification of simple radical pairs

Set `vrc_tst_active_space` to `auto_radical_pair` to make the owning VRC driver
identify, calculate, verify and install a CAS(2,2) model before production.
This is distinct from the older `radical_pair` shortcut, which only specifies
the electron/orbital count. Automatic mode has no unverified frontier fallback.
An explicit reviewed surface model takes precedence.

## Applicability

Policy version 1 supports accepted homolytic C–C or C–H cleavage of a neutral,
closed-shell singlet, acyclic saturated hydrocarbon into two neutral doublets.
It verifies the actual product graphs against the severed parent graph. The
full-valence comparator is limited to 14 orbitals including all hydrogen 1s
orbitals. Unsupported composition, bonding, charge, spin, explicit state
averages, or incompatible electronic methods produces an ordinary review
request. No species names, campaign paths, or reference energies appear in
the policy.

The production protocol must use CASPT2-F12(2,2) for scanning, sampling and a
single spin-gap level, with a triplet UCCSD(T)-F12b reference. The gap and
reference use the same two-basis CBS pair. Scan and gap shift/IPEA controls
must match. The scan begins at or above 1.8 angstrom; its largest point must
be below the separate asymptote, which must be at least 25 angstrom.

## Native work and scientific decision

1. Independent neutral-doublet ROHF fragment calculations and a triplet dimer
   calculation are run in every sampling, scan and correction basis. Complete
   native MO exports establish that the dimer's two SOMOs span both fragment
   SOMOs and preserve neutral fragment populations.
2. Explicitly selected SOMOs seed separate singlet and triplet single-state
   CASSCF calculations. The production constrained optimizer supplies relaxed
   radial geometries, with all attempts retained. Independent fragment SCFs
   are repeated at the relaxed asymptotic geometries.
3. Compact and all-atom full-valence CASSCF/CASPT2-F12 panels cover the radial
   endpoints, interior points, three angular cuts at the three closest
   separations and their asymptotes, returns at identical geometries, and
   independent starts at each cut's endpoints. Compact panels additionally
   compute a shifted CASPT2-F12 potential, MRCI-F12+Q, and triplet CCSD(T)-F12b
   with native reference-density diagnostics.
4. The assessor reads the executed inputs, native outputs, complete orbital
   XML, saved wavefunctions and reference matrices. It verifies the actual
   MULTI electron/orbital/spin/root definition. The transformed, normalized
   CAS(2,2) CI vector must agree with the exported natural-orbital density.
   At dissociation it must describe a covalent singlet or neutral triplet built
   from the independently identified fragment SOMOs. Orbital rotations are
   compared as subspaces, not by visual orbital numbering. Return-wavefunction
   overlaps include the inactive determinants as well as the active pair.
5. Interaction energies use a common asymptotic zero per model. Both individual
   bases and the configured CBS extrapolation are tested, including the actual
   reconstructed singlet potential and the sampling potential's angular error
   after radial correction. The triplet CC reference must match the MR state
   density, SOMO subspace and independent MR interaction potential.

The fixed, versioned policy records every check and observation. Its principal
bounds are:

| Check | Bound |
| --- | --- |
| Compact versus full-valence, shift, independent-method potential errors | 0.5 kcal/mol |
| Angular residual after radial correction | 0.5 kcal/mol |
| Return / independent-start energy difference | 0.02 kcal/mol |
| Far spin degeneracy, plateau and matched CASSCF/fragment energy | 0.05 kcal/mol |
| Fragment/active and compact/full subspace projection eigenvalues | at least 0.90 |
| Return subspace projection eigenvalues and squared wavefunction overlap | at least 0.99 |
| Asymptotic ionic probability / fragment charge error | at most 0.02 / 0.02 electron |
| Triplet CC T1 diagnostic / reference-density distance | at most 0.03 / 0.20 |
| Full-valence fractional occupation outside the compact subspace | at most 0.10 |

These are explicit operational acceptance bounds, not a claim that finite
probes prove exactness at every geometry or temperature. Standard production
scan, flux convergence and kinetics checks still apply. State averaging and
chemically more complicated pairs require reviewed models.

Molpro's [MULTI documentation](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi)
describes CSF configuration and CI transformation;
[AVAS](https://www.molpro.net/manual/doku.php?id=automated_construction_of_atomic_valence_active_spaces)
provides the independent all-atom valence selection. The scientific decision
also requires the native comparisons above; neither keyword alone is approval.

## Evidence, promotion and restart

Each parent worker stores an immutable campaign under
`vrctst/auto_radical_pair/<contract SHA256>/`. It includes `contract.json`,
`coverage.json`, `relaxations.json`, phase reconnaissance manifests and native
artifacts. `status.json` identifies progress. `verification.json` records
checks, numerical observations and artifact hashes. A passing decision writes
`reference_review.json` and `surface.json` with separate singlet/triplet,
basis-specific checkpoint selections. The physical surface is singlet root 1
with weight 1; the triplet is only its energy-correction reference.

Production validates the bound model, protocol, policy and native files again.
Missing, changed, partial or failing evidence cannot produce approval. A
scientific failure persists as `needs_review` so restarting the driver does
not repeat an exhausted panel. A crash after the scientific decision but
before surface publication completes publication without recalculation.

The workflow runs through the existing reconnaissance and constrained-scan
owners. It joins the PES shared resource pool, preserves live allocations and
completed attempts, and uses at most two allocation retries per probe.
Transient `squeue` failures wait and poll without consuming a chemistry retry.
An ambiguous interrupted submission requires scheduler reconciliation rather
than risking a duplicate job. Controller locks prevent competing writers.
This feature does not promise that a lost node's unsaved in-memory work survives.

To audit a phase without starting new chemistry (this updates its audit fields):

```bash
python -m kinbot.molpro_reconnaissance audit /absolute/path/to/phase/reconnaissance.json
```

Local tests exercise native-deck generation, analytic orbital/CI fixtures,
scientific failure cases, live-work preservation, scheduler outages, shared
capacity, repeat invocation and interrupted publication. They are simulated
site validation. A molecule is native-verified only after its own campaign
passes and writes an accepted `verification.json`.
