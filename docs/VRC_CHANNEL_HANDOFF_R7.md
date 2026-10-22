# r7 VRC channel and orbital handoff, 2026-10-21

This describes the three initial reviewed channels in the source r7 inputs
and the new CH2SH review recipe. The live refresh adds the latter to the
existing CH3S case input; it does not replace the initial channels. Additional
channels discovered elsewhere in the PES require their own review. This
document does not invent active spaces for them.

## Additional orbital-only batch, October 21 evening

The plan `examples/anl/orbital_review_batches/r7_additional_20261021/batch.json`
stages these **three screens, eight jobs each**, under **one shared cap of
three allocations** (reduced if the current production review reserve is
smaller). Each uses independent ROHF/AVAS CASCI at 2, 3.5, 6 and 30 angstrom
in aDZ-F12 and aTZ-F12. There are no CCSD(T), CASPT2, geometry optimization,
or potential-validation calculations in this batch. Printed roots and
natural-orbital AO coefficients/occupations support subsequent human selection.

| Screen | Source | Diagnostic space | Averaged doublet/quartet roots |
|---|---|---|---|
| C + SH3, C–S cleavage | Actual new review for `471430450000000000002-s06a7b4c221beb748` | CASCI(13,11), C/S valence plus all H 1s targets | 3 / 3 |
| CH2 + HSS, C–S cleavage | Provisional search at `791671580700000000002-s36f1e2b924862f12` | CASCI(19,12), C 2s/2p and both S 3s/3p | 2 / 2 |
| CH2SH + S, S–S cleavage | Same provisional parent, separate bond/record | CASCI(19,12), same target classes, independently generated orbitals | 3 / 3 |

These larger **diagnostic spaces are not production active-space choices**.
The CH3SS spaces retain all valence electrons in heavy-atom AVAS targets;
they do not include independent H-antibonding targets. If the printed
character indicates those are needed, expand the diagnostic space before
selecting production orbitals. The roots are diagnostic averages, not an
instruction to count that many reactive surfaces. C(3P) and S(3P) motivate
three spatial components per coupled spin; CH2 + HSS uses one extra
diagnostic root per spin beyond its nondegenerate ground-fragment channel.

Preparation reads current requests and well inputs and writes only the
separate recon directory. A provisional channel must already be an active
or accepted `hom_sci` search in `kinbot_monitor.out`; rejected searches are
not staged. Known CH2SH orbital work is not in this plan. Channels already
supplied in the case input are skipped, including a published refresh that
has not yet removed its old review marker. Repeating preparation refuses an
existing campaign; repeating `run` resumes its manifests and validated
outputs. Manifest locks prevent a second controller submitting those jobs.
Independent screens continue if one screen fails. Startup/allocation retries,
node exclusions, exclusive jobs, local scratch and partition-maximum walltime
come from the existing recon execution path.

Run from the cluster checkout (no production-controller restart):

```bash
cd "$HOME/KinBot"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
KB_BATCH="$HOME/kinbot-runs/r7_additional_orbital_recon_20261021"
if [ ! -f "$KB_BATCH/review_launch.json" ]; then
  .venv/bin/python examples/anl/tools/requested_orbital_recon.py batch \
    "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7" \
    examples/anl/orbital_review_batches/r7_additional_20261021/batch.json "$KB_BATCH"
fi
.venv/bin/python -u examples/anl/tools/requested_orbital_recon.py run "$KB_BATCH"
```

Use the same helper's `monitor "$KB_BATCH"` and `evidence "$KB_BATCH"`.
When a provisional channel gets its actual review request, bind its completed
evidence before selecting orbitals:

```bash
.venv/bin/python examples/anl/tools/requested_orbital_recon.py bind \
  "$KB_BATCH/ch2shs_cs" \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/ch3ss/791671580700000000002-s36f1e2b924862f12/vrctst/orbital_review/791671580700000000002-s36f1e2b924862f12_hom_sci_1_2/request.json"
```

Binding checks parent, atom ordering/connectivity, bond, charge, spin and
separated fragment identities. A change in optimized coordinates is recorded
and both geometries retained. Binding does not approve an active space or
change production. A reviewed selection recipe still goes through `review
--apply`, which uses the existing live refresh. S–S binding uses the
`ch2shs_ss` child and `hom_sci_2_3` request. Native evidence remains separate
if a tentative reaction is later rejected.

## D1/D2 counting at the MESS boundary

The current single-bottleneck handoff already implements Union's operation:

`N(E) = sum_i weight_i * g_i * N_i(E - offset_i)`.

Each surface contributes once, with `g_i=2` for an unsplit doublet. Explicit
fine-structure levels, if supplied for a surface, replace that default spin
factor. State-averaging weights do not enter. Conserved fragment vibrations
and rotor factors are shared by these surfaces and convolved once after
the sum; MESS receives `ElectronicLevels 0 1` for this combined barrier.
The CH3S well's separate electronic-level partition function is retained.

MESS `Union` likewise sums constituent microcanonical state counts and
Boltzmann-shifted canonical weights: see the
[official manual, section 6.3](https://tcg.cse.anl.gov/papr/codes/mess/messmanual_3_23_16.pdf)
and `UnionSpecies::states` / `UnionSpecies::weight` in the
[official source](https://github.com/Auto-Mech/MESS/blob/master/src/libmess/model.cc).
An additional Union of the two leaves **and** the already summed table
would double count. Replacing this valid sum with a second aggregation is
not required. Tests integrate unequal-threshold doublet fluxes at four
temperatures and verify the analytic Union weights; identical doublets
contribute four spin-resolved states, not eight and not their average.
The manifest reader now also rejects duplicate spin/root entries or a
repeated scalar input within one bottleneck.

This equivalence assumes the configured common conserved-mode model and
the shared sampled energy range. Distinct state-specific fragment vibrational
models would require separate convolutions before summation. It does not
add distance-dependent spin-orbit or nonadiabatic coupling; those remain
limitations of the selected spin-free adiabatic model, not problems solved
by writing `Union`. Existing stereochemical barrier Unions remain separate:
they combine distinct spatial reaction pathways after their own duplicate
classification.

## Ground terms, surfaces and averaging are different quantities

CH3S has a ground electronic term X²E. At C3v and without spin–orbit
coupling, E has two degenerate orbital components; the superscript 2 denotes
spin multiplicity. This is an electronic-structure property, not a selectable
orbital count. Symmetry-lowering Jahn–Teller distortions can split the orbital
components. Spin–orbit coupling changes the physical levels and competes
with the distortion. Two exactly equal energies everywhere along a general
scan are not implied. See [Marenich and Boggs, 2004](https://doi.org/10.1021/jp040443f).

Treating those components as two independent, uncoupled spin-free adiabatic
VRC surfaces is a **model choice**. The production input selects doublet roots
1 and 2 for CH3S. It does not implement coupled spin-vibronic dynamics or a
distance-dependent spin–orbit Hamiltonian. The chosen roots must retain the
intended character; root numbers alone do not prove that they do.

CH3SS's ground-parent term is X²A″ in Cs, which is orbitally nondegenerate;
see the [photoelectron assignment](https://doi.org/10.1063/1.1343081).
Consequently its two chemical dissociations each select one ground-parent
doublet surface. The S–S products nevertheless include CH3S(²E) and S(³P),
so that channel needs a much larger *state average*. Product-manifold
degeneracy does not make every asymptotic component a separate reaction into
the ground-parent well. Selecting only the ground-parent adiabatic pathway
is not proof that excited-surface or nonadiabatic contributions are zero.

| Chemical channel | Reactive scalar relaxed paths | Active space | Averaged roots (doublet / quartet) | Relative input weight per root (doublet / quartet) |
|---|---|---|---|---|
| CH3S -> CH3 + S, C–S | D1 and D2 | CAS(5,4) | 3 / 3 | 1/3 / 1/3 |
| CH3SS -> CH3 + S2, C–S | D1 | CAS(7,5) | 1 / 1 | 1 / 1 |
| CH3SS -> CH3S + S, S–S | D1 | CAS(7,5) | 6 / 6 | 1/6 / 1/6 |
| CH2SH -> CH2 + SH, C–S, new review | D1 | CAS(5,4) | 4 / 2 | 1/4 / 1/2 |

D1/D2 mean roots within spin multiplicity 2, not the two spin projections of
one state. Every reactive leaf has spatial weight 1 and a common spin-free
ground-fragment asymptotic reference. Averaging weights normalize to a total
of 1/2 per spin block: each CH3S root has total weight 1/6, each CH3SS C–S root
1/2, each CH3SS S–S root 1/12, and each new CH2SH doublet/quartet 1/8 and 1/4.

State averaging optimizes one compromise orbital set across the specified
states. It reduces bias toward a single component and helps avoid orbital
optimization root flipping; it is not a guarantee against crossings or an
occupation-based proof of state identity. [Molpro's state-averaging guidance](https://www.molpro.com/manual/doku.php?id=quickstart#state-averaged_mcscf)
describes this distinction. These weights are not thermal populations,
spin-statistical weights or an instruction to sum all averaged roots as flux.
The VRC assembly separately combines the declared reactive leaves, includes
their spin factors, and writes its combined MESS barrier with electronic
degeneracy 1 to avoid multiplying that factor again.

The occupations below are eigenvalues of the state-averaged one-particle
density, not populations of individual CASPT2 states. A value such as 1.5
can arise from averaging complementary occupations. AO coefficients locate
an orbital on particular atoms and show its angular/bonding character, but
their squares are not population percentages in a nonorthogonal AO basis.
The overall orbital sign is arbitrary; relative phases and the complete
subspace matter. Occupations alone are therefore insufficient for selection.
Each index range below refers to its own reviewed WFU file and record 2140.2,
not to a universal ordering shared by the molecules or channels.

## CH3S C–S: five electrons in four orbitals

CH3 contributes its singly occupied carbon frontier orbital (one electron),
and atomic sulfur contributes its complete 3p manifold (four electrons in
three orbitals). All three sulfur p directions are retained because angular
sampling must not privilege one atomic component. S(³P) has three spatial
components; coupling its triplet spin with doublet methyl yields three
doublets and three quartets in the spin-free ground-fragment manifold.
The parent X²E components motivate two reactive doublet paths; the third
averaged doublet is not automatically another ground-parent flux path.

The reviewed DZ/TZ natural records use indices **11–14**. Larger-parent
coefficients place 11–13 on sulfur p and 14 on methyl carbon. Their confirmed
far occupations are `(1.33333, 1.33333, 1.33333, 1.00000)` in both bases. At
the inner TZ point they become approximately
`(1.49748, 1.49748, 1.47637, 0.52867)`, showing bonding redistribution inside
the retained space. C–H bonds and sulfur 3s are inactive.

Adding sulfur 3s in CAS(7,5) changed the tested doublet angular residual by
0.00039 kcal/mol and the doublet–quartet angular-gap residual by 0.01640.
The occupied-valence (13,8) comparison gave 0.10803 and 0.03367. These support
the minimum choice on the tested cuts, not a universal accuracy bound.
[CH3S evidence and native inventory](../examples/anl/ch3s_anl1_f12/orbital_reviews/ACTIVE_SPACE.md).

## CH3SS C–S: seven electrons in five orbitals

The fragments are methyl doublet and ground triplet S2. Keep methyl's
one-electron frontier orbital plus the complete S2 pi/pi-star manifold:
six electrons in four orbitals. S2's intact sigma bond is outside this
minimum CAS. Unlike atomic S(³P), the ground S2 term is an orbitally
nondegenerate Sigma term; there is one spatial ground-fragment state for
each of the doublet and quartet spin couplings, hence one root per spin.

The larger parent record requires original indices **17,18,20,21,22**.
Preparation swapped 17 and 19 so the saved production records contain the
chosen set contiguously at **18–22**. Do not apply that swap again. Bonding pi
coefficients have corresponding same-phase sulfur contributions; the pi-star
partners have opposite phases. The methyl orbital is carbon-dominated.
Confirmed far TZ occupations are
`(1.95828, 1.95828, 1.04172, 1.04172, 1.00000)`.

The omitted sigma/sigma-star pair is not exactly uncorrelated: its larger
parent occupations are about 1.94564/0.06127. The explicit CAS(9,7) expansion
changed the tested doublet angular residual by 0.00239 kcal/mol and the
angular spin-gap residual by 0.00095, supporting its exclusion for this
channel. [CH3SS evidence](../examples/anl/ch3ss_anl1_f12/orbital_reviews/ACTIVE_SPACE.md).

## CH3SS S–S: a different seven-electron, five-orbital space

Keep the CH3S sulfur-centered frontier pair (three electrons) and all three
departing atomic sulfur p orbitals (four electrons). The two CH3S orbitals
support complementary (2,1)/(1,2) occupations. Freezing one as a permanent
lone pair would remove the other hole configuration. This is the reason for
CAS(7,5), rather than atomic S plus only one CH3S radical orbital.

The two CH3S spatial components times the three atomic S components give
**six spatial roots per coupled spin**. That motivates six doublets and six
quartets in the average. They need not be all degenerate at a distorted
fragment geometry. Only the selected ground-parent D1 path contributes
reactive flux in this model.

The channel-specific natural seeds also use indices **18–22**, but their
content differs from the C–S channel. Indices 18/19 belong mainly to bound
sulfur, and 20–22 to departing sulfur. Both bases give far occupations
`(1.50000, 1.50000, 1.33333, 1.33333, 1.33333)`. At the inner TZ point they
become `(1.71816, 1.63491, 1.49536, 1.36269, 0.78888)`. Canonical orbitals can
reorder these fragment groups inside the same CAS without losing the span.

Adding the intact C–S sigma/sigma-star pair in CAS(9,7) changed tested angular
doublet/gap residuals by 0.04494/0.05703 kcal/mol. At 2.1 Å, corresponding
interaction/gap differences reached 0.114/0.316 kcal/mol. There is no claim
of uniformly sub-0.1 kcal/mol agreement. The tested active-span overlaps and
return calculations support continuity, not universal tracking at every
future geometry. [Detailed native evidence](../examples/anl/ch3ss_anl1_f12/orbital_reviews/ACTIVE_SPACE.md).

## New CH2SH C–S review: five electrons in four orbitals

Triplet CH2 has two frontier electrons in two orbitals; SH(²Pi) contributes
three electrons in its pi pair. This gives CAS(5,4), source indices **11–14**
in each basis's original CASCI record. SH contributes two spatial ground
components, yielding two ground-fragment doublets and two quartets. The
four-doublet/two-quartet average deliberately includes two extra diagnostic
doublets; their energies/occupations alone do not identify their asymptotes.
The selected reaction into the nondegenerate parent uses D1 only. The new
recipe retains explicit Q1 for the MR reference; the earlier three channels
use their reviewed lowest-correlated-quartet selector because their native
quartet root ordering was observed to change.

At 2 Å, the completed aDZ CASSCF(5,4) occupations are
`(1.70808, 1.48606, 1.02341, 0.78245)`, with the intended SH/CH2 AO character.
The CAS(7,6) checks at 2, 3.5 and 30 Å leave the added S–H sigma pair near
1.98/0.02 and preserve the frontier pattern. Both original CASCI bases give
consistent near/far character. This supports orbital selection without
completing the optional 22-job potential panel. Full two-basis optimized-CAS
and angular comparisons are not claimed. [Evidence, recipe and refresh command](CH2SH_ORBITAL_REVIEW.md).

## Shared production sequence and controls

Each reviewed reactive leaf uses its own chosen CAS, state average and
basis-specific hashed natural-orbital seed. The source inputs prescribe:

1. A constrained separated-supermolecule optimization at the 30 Å reference,
   using CASPT2-F12/aDZ-F12; extract the fragment geometries. This reuses the
   full channel CAS and requires no separate fragment active spaces.
2. Fixed-geometry L2 Gaussian B2PLYP-D3BJ/cc-pVTZ fragment orbital jobs for
   HOMO pivot directions. These do not replace MR fragment geometries or
   enter the MR interaction energy. Atomic fragments need no such job.
3. CASPT2-F12/aDZ-F12 relaxed paths and matched frozen-fragment calculations,
   retaining the existing distance constraints, without imposing new C3v/Cs
   geometric symmetry. The 47 distances are 2–8 Å by 0.2, 8.5–12 by 0.5,
   and 13–20 by 1. The separate asymptote remains 30 Å.
4. High-level reference curves from ROHF-UCCSD(T)-F12b on the quartet plus
   CASPT2-F12 doublet-minus-quartet gaps. Both the CC reference and MR gaps
   use aDZ/aTZ-F12 extrapolation with power 3.7. In Molpro, `rhf` is the
   restricted-open-shell HF reference for these open-shell calculations.
   For one gap method the construction is
   `E_target = E_CC,quartet,CBS + E_MR,target,CBS - E_MR,quartet,CBS`.
5. Assemble the relaxation/high-level correction and perform CASPT2-F12/aDZ
   orientation sampling with RotdPy at the same 47 dividing distances.
   The relative flux sampling target is 5%; configured sample bounds still
   apply. Several pivot alternatives may exist at a distance, so 47 is not
   necessarily the number of RotdPy Surface objects.

CASPT2 uses IPEA 0.25 Eh, real shift 0.1 Eh and the reviewed same-spin
state-average reference policy. The earlier three channels have separate
DZ/TZ WFU files already in their case directories. The new recipe points to
the completed reconnaissance's two hashed seeds. The source configuration
therefore has four reviewed chemical channels and five reactive leaves after
this addition; it does not describe every channel future exploration may find.

The live command preserves the actual input's grids, levels, numerical and
resource controls rather than reinstalling defaults. Existing scratch-first,
exclusive-node, partition-walltime and shared-admission policies are untouched.
Other reviewed paths continue; a new unspecified channel remains held for
review. ANL still waits for outer-PES L2 selection. Neither successful orbital
review nor a passing parser certifies eventual convergence of every native
point or the complete kinetics calculation.
