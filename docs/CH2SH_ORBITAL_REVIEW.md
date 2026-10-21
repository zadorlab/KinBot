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

The prepared recipe and existing native transformations were checked locally.
The completed native CASCI screen is recorded below; a production active space
has not yet been validated. After the screen:

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

## Completed CASCI screen and next comparison

The user supplied eight completed native records on 2026-10-21. All passed
the controller's checks without warnings. The state-averaged natural table
and labeled AO coefficients identify the following **source** record indices
in both bases (record 2140.2, C1):

| Candidate | Source indices | Content |
|---|---|---|
| CAS(5,4) | 11–14 | SH pi pair and CH2 frontier pair |
| CAS(7,6) | 10–15 | Same space plus S–H sigma/sigma-star (10 and 15) |

At 30 angstrom the SH pi occupations are approximately 1.4997 each;
the CH2 pair has approximately 1.206/0.794 occupations in this diagnostic
state average. The source orbitals 13 and 14 exchange character between
near and far geometries; the selected *subspace* includes both. No individual
orbital identity follows from occupation order. The lowest two doublets
match their respective quartet partners at separation to printed precision.
The two spatial partners are not exactly degenerate in the unoptimized
CASCI screen; this is a reason to examine optimized orbitals, not to assign
an asymptotic splitting to production.

`coverage_refinement_required` is not native failure. With only 2, 3.5, 6,
and 30 angstrom anchors, the generic three-point plateau check includes
3.5 angstrom. Its suggested 42 angstrom extension does not establish that
30 angstrom is too short. The candidate panel below instead adds 20 angstrom
and the requested 2.75 angstrom midpoint. It remains a finite diagnostic,
not automatic scientific acceptance or exhaustive orientation coverage.

`ch2sh_cs_candidates.json` specifies 22 new fixed-geometry electronic-structure
calculations: CAS(5,4) at 2, 2.75, 3.5, 6, 20 and 30 angstrom in both bases
(12), and CAS(7,6) at matched 2, 3.5 and 30 angstrom in both bases (6).
Both candidates also test a 30-degree rigid rotation of SH about the C–S
axis at 3.5 angstrom in both bases (4). This preserves fragment bond lengths
and C–S separation and checks sensitivity to one changed approach orientation.
Every job starts from a basis-specific saved CASCI record; the new midpoint
and rotated geometry use the 3.5 angstrom seed, and 20 angstrom uses the
30 angstrom seed. Only electronic orbitals are optimized in these probes.
No completed screen is rerun or overwritten. Both candidates use the same
two-quartet/four-doublet state average to isolate the space comparison.

Each job checks CASPT2-F12 using the production 0.1-hartree real shift,
0.25-hartree IPEA shift and state-averaged reference policy. The compact
panel also checks the quartet ROHF-UCCSD(T)-F12b reference. The expanded
panel does not repeat that identical single-reference calculation. Compare
interaction energies against each candidate's own asymptote and spin gaps;
different active spaces' absolute energies are not an accuracy metric.

Prepare using `requested_orbital_recon.py candidates SOURCE RECIPE OUTPUT`.
Preparation checks saved input hashes, request identity, active counts,
WFU/output hashes and every generated deck before publishing. The common
`run OUTPUT` controller manages both candidate manifests under **one** cap
of three allocations, retaining each manifest's submission lock and native
first-job gate. `monitor OUTPUT` and `evidence OUTPUT` are read-only.
There is no automatic model approval, production-input edit, or driver restart.

Once the optimized/correlated comparison and required orbital/state/angle
checks support a choice, the live update will add a new parent entry under
`vrc_tst_surface_models["471030810040000000002"]` with basis-specific
WFU hashes, the chosen CAS, state average and selected reactive root(s).
Use `kinbot.run_control resubmit` on the existing case input; it publishes
an additive review revision to the running PES owner. The native comparison
must finish before this final active-space selection is supplied.

Run the candidate panel on Blodgett (a repeated invocation reuses the
prepared panel; the manifest locks reject a concurrent second controller):

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
KB_PY="$PWD/.venv/bin/python"
KB_TOOL="$PWD/examples/anl/tools/requested_orbital_recon.py"
KB_SOURCE="$HOME/kinbot-runs/r7_ch2sh_cs_orbital_recon"
KB_OUT="$HOME/kinbot-runs/r7_ch2sh_cs_orbital_validation"
if [ ! -f "$KB_OUT/review_launch.json" ]; then
  "$KB_PY" "$KB_TOOL" candidates "$KB_SOURCE" \
    "$PWD/examples/anl/ch3s_anl1_f12/orbital_reviews/ch2sh_cs_candidates.json" "$KB_OUT"
fi
nohup "$KB_PY" -u "$KB_TOOL" run "$KB_OUT" \
  >>"$KB_OUT/controller.log" 2>&1 </dev/null &
printf 'Review controller PID: %s\nLog: %s/controller.log\n' "$!" "$KB_OUT"
BASH
```

Monitor the panel without changing it:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
"$HOME/KinBot/.venv/bin/python" \
"$HOME/KinBot/examples/anl/tools/requested_orbital_recon.py" monitor \
"$HOME/kinbot-runs/r7_ch2sh_cs_orbital_validation"
```

When it finishes, replace `monitor` with `evidence` to print the saved
CASSCF/CASPT2 energies, active occupations and validation flags for review.
