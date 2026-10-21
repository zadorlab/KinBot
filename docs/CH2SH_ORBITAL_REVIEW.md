# CH2SH C–S orbital review, 2026-10-21

The r7 request is `471030810040000000002_hom_sci_1_2`, from the
`471030810040000000002` well. The saved request identifies neutral doublet
CH2SH, atom 1 C and atom 2 S, dissociating to triplet CH2 and doublet SH.
This is different from the approved initial CH3S -> CH3 + S channel.

## Candidate spaces and state identities

The first compact candidate is CAS(5,4): two CH2 radical orbitals/electrons
plus both SH pi orbitals with three electrons. SH has a ground 2Pi term;
both spatial components belong in the separated-fragment orbital space.
The C–S sigma/sigma-star pair forms by mixing the fragment orbitals already
counted here. Adding that pair again would double-count the orbital degrees
of freedom. CAS(3,5) does not automatically include the occupied SH pi partner:
adding virtual orbitals alone cannot replace its two active electrons.

The first expansion candidate is CAS(7,6), adding the S–H sigma/sigma-star
pair. These are chemically motivated candidates, not selected production
spaces. Actual orbital coefficients, occupations, state character and
correlated potentials must establish which orbitals are retained and whether
other expansions are needed.

The initial screening recipe runs full-valence CASCI(13,11), including C
2s/2p, S 3s/3p and all three H 1s targets. Its purpose is orbital discovery,
not production sampling. Quartet ROHF supplies a high-spin starting guess;
the production parent's doublet multiplicity is untouched. The CASCI lists
two quartet and four doublet roots. The extra doublets allow inspection of
other nearby states, including possible singlet-CH2-derived character;
their integer energy ranks do not prove their dissociation identities and
do not authorize additional reactive flux surfaces.

There are eight independent calculations: 2.0, 3.5, 6.0 and 30.0 angstrom
at avdz-f12 and avtz-f12. Only the fragment separation changes; these are
fixed-geometry probes, not fragment optimizations or relaxed scans. A first
native job must validate before the remaining jobs can use the three-slot
cap. The gate reuses that job's result: no duplicate pilot calculation.

Molpro recommends an initial AVAS/CASCI occupation inspection before manual
space reduction. Occupations alone are insufficient; retain chemical and
state-character evidence as described in the general review procedure.
References: [Molpro AVAS manual](https://www.molpro.net/manual/doku.php?id=automated_construction_of_atomic_valence_active_spaces),
[SH ground-state terms](https://arxiv.org/abs/1909.02646),
[general review procedure](ORBITAL_REVIEW_WORKFLOW.md).

## Execution and continuation

`examples/anl/tools/requested_orbital_recon.py prepare` takes four arguments:
the live request, its case input, a screening recipe, and a fresh output
directory outside the production pool. It verifies the intended request,
checks all prepared native decks, snapshots the geometry/electronic context,
and inherits the current partition and exclusions. Preparation never submits.
Repeating preparation into the same directory fails without overwriting it.

`run OUTPUT` uses the existing journaled reconnaissance controller, at most
three jobs and no more than the currently configured review reserve. Each
job is exclusive, uses eight Molpro ranks with 1000 MW/rank and requests
104000 MB total, and defaults to the partition maximum walltime. Native
execution uses the existing checked node-local scratch workspace. Only
site setup / executable discovery and scheduler queries run on the login
node. Two allocation-recovery attempts are enabled for recoverable
infrastructure failures. Scientific failures remain failures for review.

`monitor OUTPUT` is read-only and shows the scheduler, native logs and saved
validation results. Both controller entry paths use the existing manifest
lock, so repeating `run` cannot create a second set of allocations.

The prepared recipe and existing native transformations were checked locally;
this new channel has not yet been validated by native Molpro. After the screen:

1. Inspect basis-specific orbitals and assign the proposed compact spaces.
2. Compare optimized compact/expanded CASSCF models, independent starts,
   root character and representative angular configurations.
3. Confirm CASPT2-F12 with the production shift/IPEA, both bases, spin gaps
   and matched interaction-energy references. Expand checks if warranted.
4. Add only the reviewed new-parent entry to `vrc_tst_surface_models` in the
   existing case input, then use `kinbot.run_control resubmit ... --apply`.

This screen does not accept a model, release the review hold, change the
production input, or restart either PES driver. Existing approved VRC paths
and exploration continue throughout. The review reserve is an admission
limit, not a hardware reservation against other cluster users.
