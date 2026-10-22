# CH + H2S orbital selection, 2026-10-22

This is the orbital-model selection for
`471230610040000000002_hom_sci_1_2` in the CH3S r7 PES. It does not approve a
finished potential, scan, or rate. The user supplied eight complete CASCI
records plus native AO coefficients at 2 and 30 angstrom in both bases.
No further reconnaissance QC is required to install this initial selection.

## Geometry and evidence

The live request supplied neutral doublet parent geometry with atom order
C(1), S(2), H(3), H(4), H(5), C-S cleavage [1,2], and products CH doublet plus
H2S singlet. H(5) belongs to CH and H(3)/H(4) to H2S. The source campaign is
`$HOME/kinbot-runs/r7_ch_h2s_orbital_recon_20261022`:

* `request.json`: source atom order, fragment states and geometry.
* `native/reconnaissance.json`: all eight CASCI job records, occupations,
  energies, wavefunction and native-output hashes.
* `native/{adz,atz}__cs_r_{2,3p5,6,30}/`: native inputs, outputs, WFU records,
  `execution.json` and referenced orbital-audit reports.

These were independent ROHF/AVAS CASCI(13,11) screens, using four equally
weighted doublet roots, at 2, 3.5, 6 and 30 angstrom in avdz-f12 and avtz-f12.
The source request was used without substituting another well's geometry.
The live publisher checks the saved/request geometry and exact basis-specific
WFU/output hashes before publishing the selection. Production uses its current
geometry and performs its own optimization.

## Selection and chemical interpretation

Use **CAS(5,4)** with source indices **8,12,13,14** from each basis's **2 A**
natural-orbital record **2140.2**. This is two sulfur donor electrons plus
three electrons in the CH frontier subspace. CAS(3,3) using only 12-14 would
omit the explicitly active sulfur donor needed to describe C-S formation.

| Character | 2 A source index | 30 A source index | DZ occupations, near / far | TZ occupations, near / far |
|---|---|---|---|---|
| S p-type lone-pair donor, mixing with C at short range | 8 | 7 | 1.97918 / 1.99881 | 1.97918 / 1.99881 |
| CH sigma-type frontier orbital | 12 | 12 | 1.48472 / 1.47041 | 1.48447 / 1.47035 |
| First CH p frontier component | 13 | 13 | 0.99467 / 0.76629 | 0.99490 / 0.76637 |
| Second CH p frontier component, strongly mixed at short range | 14 | 14 | 0.50273 / 0.76595 | 0.50297 / 0.76594 |

The donor's ordering changes: DZ source 8 at 2 A has leading S px/py/pz
coefficients 0.64401/-0.31925/0.16919 and C s 0.36763. At 30 A it is source
7, with S px/py/pz 0.76264/-0.58958/-0.26605. The S s-dominated near orbital
7 becomes far orbital 8. The corresponding TZ coefficients show the same
exchange of character. Applying index 8 to the far record would therefore
choose the wrong sulfur orbital. Only the 2 A records are production seeds.

At 30 A, source 12 is C/H(5) sigma-like, while 13/14 are C p combinations;
their leading coefficients contain no appreciable sulfur term in the supplied
tables. At 2 A, 14 has substantial S/H(3) mixing. Source 9 is the occupied
C-H bond combination; 10/11 and 15/16 belong mainly to H2S bond/antibond
combinations. Whole-orbital signs are arbitrary and AO coefficient squares
are not atom-population percentages in a nonorthogonal basis.

Molpro repacks source 8 into active slot 11, leaving 12-14 in place; ten
orbitals are inactive (`closed,10`, `occ,14`). The selected span, not its
occupation order at a later geometry, defines this seed. Subsequent accepted
checkpoints carry the active block forward. All other channels' seeds remain
unchanged.

## State average and correction

Retain **four doublet roots at weight 1/4 each**, as in the screen. Select
**D1 only** as this ground-parent reactive path. D2-D1 at 30 A is 2.7566 cm-1
(DZ) or 3.5138 cm-1 (TZ), compared with 0.09034533/0.09020927 Eh at 2 A.
The near merging is consistent with the CH ground orbital pair, but finite
distance CASCI energies do not prove exact degeneracy or path connectivity.
CH's ground X2Pi assignment is listed in the
[NIST microwave spectral database](https://physics.nist.gov/cgi-bin/MolSpec/hydrosearch.pl?fixedpoint=1&lowerfreq=&molecule=CH&units=MHz&upperfreq=&x13CHxXx2Pi1sh2=1&xCHxXx2Pi1sh2=1&xCHxXx2Pi3sh2=1).
Additional averaged roots are orbital support, not additional reactive flux.
Selecting D1 does not establish that excited-surface contributions are zero.

The specified fragment spins (1/2 and 0) allow only total spin 1/2. Therefore
this channel uses the existing **direct CASPT2-F12 correction** branch; a
quartet CC anchor and doublet-quartet gap would not describe these ground
fragment terms. The selected recipe makes `vrc_tst_correction_model=direct`
explicit. Scan, sampling and trusted levels use CASPT2-F12(5,4), retaining
the live input's bases, DZ/TZ extrapolation, grids, shifts and convergence
controls. Existing channels keep their own spin-gap treatment.

## Limits that remain visible

Source 15 is mainly S-H antibonding, with occupation 0.08577 (DZ) / 0.08553
(TZ) at 2 A, falling to 0.03627 / 0.03628 at 30 A. It is not the missing
sulfur donor. Leaving it and its bond partner inactive/virtual is a compact
model choice, **not** a demonstrated active-space convergence result. Inspect
S-H behavior, active-orbital character and root continuity in production;
retain all connectivity and whole-scan review gates. A larger S-H-containing
space remains a possible scientific revision if the production evidence
warrants it. Coefficient/occupation agreement between bases does not validate
angular accuracy or every production geometry.

`coverage_refinement_required` is retained in the source manifest. Its
suggestions (2.75 and 42 A) arise from potential/occupation coverage checks;
they do not turn successful fixed-geometry orbital probes into failed jobs,
nor certify an asymptotic plateau. No compact optimized-CASSCF or correlated
validation panel has been run for this selection. This follows the manual's
[AVAS screening workflow](https://www.molpro.net/manual/doku.php?id=automated_construction_of_atomic_valence_active_spaces)
while distinguishing orbital choice from potential validation.

## Live handoff

After pulling the selection recipe, publish it through the existing owner:

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
.venv/bin/python examples/anl/tools/requested_orbital_recon.py review \
  "$HOME/kinbot-runs/r7_ch_h2s_orbital_recon_20261022" \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/ch3s/input.json" \
  examples/anl/ch3s_anl1_f12/orbital_reviews/ch_h2s_cs_selected.json --apply
.venv/bin/python -m kinbot.run_control resubmit \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/ch3s/input.json" --status
BASH
```

This adds the new parent's reviewed channel through input refresh. It starts
no replacement driver and submits no further reconnaissance jobs. The owner
may initially acknowledge `queued` and applies the addition at its refresh
boundary. `applied` and the channel's model-acceptance receipt establish the
handoff; neither means native VRC work is complete. The separate C + SH3
whole-scan review remains unresolved.

Local handoff validation: **115 passed** across the orbital-review helper,
multistate VRC, correction construction and live input-refresh tests. The
new selection test exercises basis-specific source hashes, donor repacking,
doublet-only direct DZ/TZ components, unchanged prior channels and repeated
publication without starting a driver or submitting reconnaissance QC. These
are software checks, not additional native electronic-structure validation.
