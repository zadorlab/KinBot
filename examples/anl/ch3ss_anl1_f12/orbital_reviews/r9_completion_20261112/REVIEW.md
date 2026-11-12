# Lean R9 orbital review, 2026-11-12

The supplied `orbital-evidence_nov10_1735hrs.tar.gz` has SHA-256
`bdb511e0193693d4ac614520ffdce758e9c066bf4c8173298009669e41c14e11`.
All 5,901 members, totaling 8,601,273,471 bytes, matched the embedded manifest.
The collector reported zero issues. The selected native review covers 184
correlated points, including sixteen earlier points reused by the campaign.
All 184 have completed native records and pass the input/output identity checks.

## Decisions for the retained channels

| Channel | Decision | Additional work in this campaign |
|---|---|---|
| CH3S to CH3 + S | Retain CAS(5,4) | None |
| CH3SS to CH3 + S2 | Retain CAS(7,5) | None |
| HSCH2S to CH2S + SH | Retain the chemically motivated CAS(7,6); correct angular state following | 16 fixed-point checks |
| CH2S/SH complex to CH2S + SH | Retain CAS(5,4) | 8 checks binding the selected model to the R9 atom order and correction coordinate |
| HCS/H2S complex to HCS + H2S | Retain CAS(3,3) for the outer capture channel | None |
| CH2SH to CH2S + H | Exact model selection pending | New odd-electron candidate campaign |
| CH3SS to CH2SS + H | Exact model selection pending | New odd-electron candidate campaign |

The remaining six unfinished relaxation diagnostics receive zero optimization
retries. Their completion is unnecessary for the current active-space decisions.
Production geometry optimization remains part of the later VRC workflow.

### Previously accepted C-S models

The original CH3S and CH3SS native evidence archives were reverified against
their inventories, covering 50 and 110 members respectively. Thirty documented
occupation vectors were independently matched to their native printed tables.
`accepted_review_check.json` and `native_occupation_crosscheck.json` record the
hashes and native line numbers. The respective `ACTIVE_SPACE.md` documents and
joint minimum-space decisions were reviewed alongside these records.

- CH3 + S requires the complete atomic sulfur 3p shell and methyl radical
  orbital. The five electrons occupy four orbitals, with three sulfur
  occupations of 1.33333 and the methyl occupation of 1.00000 at separation.
  The retained natural indices are 11–14. Existing S-3s and larger-space
  comparisons give doublet angular residuals of 0.00039 and 0.10803 kcal/mol.
  The three-doublet/three-quartet orbital support and the original two reactive
  doublet components remain the selected definitions.
- CH3 + S2 uses six electrons in the complete S2 pi/pi-star manifold and the
  methyl radical electron. The seven electrons occupy five orbitals, natural
  indices 18–22 in the already reordered checkpoints. The aTZ separated
  occupations are 1.95828, 1.95828, 1.04172, 1.04172 and 1.00000. The existing
  sigma/sigma-star enlargement changes the tested doublet angular residual by
  0.00239 kcal/mol. One doublet and one quartet support the orbitals, with D1
  supplying the retained reactive contribution.

The historical comparisons include explicitly documented shift-0.3 work and
matched shift-0.1 checks. The R9 production definition uses shift 0.1 and IPEA
0.25 consistently. The old comparisons support these chemical space choices
over their sampled geometries. They supply no uniform production-rate bound.

### Case 7, HSCH2S to CH2S + SH

The extra pair in CAS(7,6) remains necessary under the prior criterion of
relative interaction-potential accuracy. Removing it changes the near potential
by 1.0982/1.1140 kcal/mol in aDZ/aTZ. At the completed relaxed near geometry the
changes reach 1.1394/1.1574 kcal/mol. The aTZ retained radial return residual at
6 angstrom is 0.0000439 kcal/mol, and the near return residual is zero at the
printed precision. The 20-to-30 angstrom change is 0.000822 kcal/mol.

The angular tests expose an electronic-branch problem. The unrotated aTZ
asymptote has a D2–D1 gap of 0.000232 kcal/mol. The two angular asymptotes have
gaps of 26.088487 and 26.088261 kcal/mol. The discrepancy is also present in the
native CASSCF energies before perturbative correlation. These solutions fail
the intended SH Pi-component correlation.

Four corrected chains use the verified radial 30-angstrom checkpoint, covering
both orientations and both bases. Each starts at the rotated asymptote, moves
inward, returns to that asymptote and repeats the inward point. Their sixteen
calculations test branch reproducibility and the separated-fragment degeneracy.
The established radial calculations and the rejected smaller CAS remain source
evidence. The CAS(7,6) choice remains conditional on this angular-state check.

### Cases 8 and 13, CH2S/SH outer approach

Case 13 closely matches the retained R9 complex, with a proper-map RMSD of
0.0005295 angstrom after exchanging the two equivalent carbon hydrogens.
Its CAS(5,4) contains CH2S pi/pi-star and both SH frontier components. The aTZ
near occupations are 1.890571, 1.500000, 1.499949 and 0.109480.

The candidate has stable radial returns and the expected nearly degenerate
SH components at all three tested asymptotic orientations. The occupied-union
comparison changes the own-asymptote potential by at most 0.005083/0.004763
kcal/mol in aDZ/aTZ, including the completed relaxed-midpoint check. The aTZ
20-to-30 angstrom change is 0.002303 kcal/mol. Its 6-angstrom return difference
is 0.0000565 kcal/mol.

The opened source orbital 18 has only approximately 2–3e-5 projection into the
optimized enlarged active subspace. The original four frontier directions stay
in that subspace. Thus the enlarged solution supplies a nearby occupied-space
sensitivity check, while its energy difference does not measure an enforced
source-18 donor contribution. The earlier native chemical-character and fixed
initialization evidence, together with the current stable frontier span and
matched potentials, supports retaining CAS(5,4) for this outer channel.

The auxiliary case 8 uses a different approach geometry. Its first angular
branch has substantial electronic-branch anomalies, including an approximately
150 kcal/mol aTZ comparison difference. Those data are explicitly retained in
`numerical_review.json`. They contribute no positive validation of the R9
candidate. Case 13 supplies the directly relevant evidence.

The eight new binding checks use the actual R8 L2 complex in its R9 atom order
and the C(1)–H(4) correction anchor. They cover the starting endpoint,
intermediate distance, separated limit and return in both bases. These checks
provide the exact new checkpoint records needed by the selected reaction.

### Case 25, HCS/H2S outer approach

Retain the HCS occupied pi/SOMO/pi-star CAS(3,3), natural indices 20–22, with
one doublet D1. The aTZ near occupations are 1.899430, 1.000000 and 0.100570.
The H2S fragment remains intact along this outer capture coordinate. Its internal
H-transfer reaction is represented by the separate inner abstraction saddle.

The candidate passes the tested radial returns, angular comparisons and
separated-limit checks. The aTZ maximum return residual is 0.0002071 kcal/mol,
and its 20-to-30 angstrom change is 0.002215 kcal/mol. The optimized enlarged
trial differs by at most 0.01056 kcal/mol in the matched interaction potentials.

The optimized enlarged trial again expels the intended H2S donor, with source-19
projection below 8e-7 in both bases. It therefore provides an alternative
optimized HCS correlation comparison. The earlier fixed-orbital control that
actually retained the H2S donor changed the near energy by 1.7e-7/1.5e-7 hartree.
These observations, the chemical identity of the intact H2S partner, and the
retained HCS active character support the compact outer-channel choice.

Both failed case-25 optimizations stopped at a bookkeeping defect. Sella added
one coordinate dummy for the nearly linear C...H-S approach coordinate, and
the checkpoint writer
mistook that appended X for a changed physical molecule. Both captured frames
contain six physical atoms, finite energies and finite physical forces. The
shared checkpoint and continuation readers now remove only verified appended
zero-force X entries while retaining physical atom order and evaluated data.
The native frames and a real Sella linear-molecule execution exercise the fix.
This repair supports future calculations. Those two diagnostics receive no retry.

## H-loss campaign and final decision

The two parents have 25 and 41 electrons and doublet spin. With doubly occupied
inactive orbitals, `N_active = N_total - 2 N_closed` is odd. Their conventional
models therefore test CAS(1,1), CAS(3,3), and chemically relevant unions of the
valence pairs. The separate two-radical singlet/triplet CAS(2,2) workflow has
different electron parity.

The new campaign uses native R8 L2 parent and closed-fragment geometries,
explicit atom mappings, an independent atomic-H orbital, and counted
heavy-atom-valence-plus-departing-H references CAS(13,9) and CAS(19,13).
It checks isolated fragment singlet/triplet energies and quartet asymptotes.
The H SOMO is identified by an actual same-basis orbital projection. When
singly occupied natural orbitals are degenerate, the test retains their full
odd-electron span and verifies its hydrogen projection and electron count.

Every occupied-like/virtual-like fragment valence-pair combination in the reference is
screened with the H SOMO at near, intermediate and far geometries in both bases.
Own-asymptote interaction differences allocate the later diagnostic work.
Equivalent optimized screen spaces are recognized using same-geometry MO
subspace overlaps and energies. Up to three distinct leading pair trials,
their distinct unions, and the one-orbital limit receive radial/return/angular
and bound-parent checks. State-average sensitivity is tested for the leading
pair and largest retained union.

The initial 92 points include the 24 mandatory existing-channel checks. There
are at most 246 compact screen points and 264 further path/state checks, giving
at most 602 distinct native points. Reused screen points and equivalent spaces
reduce that count. The controller executes these dependencies in one submitted
campaign under the existing 20-job pool. There are zero optimization retries.

The ranking allocates diagnostics only. Final selection still uses native
orbital/CI character, proper fragment-state correlation, branch continuity,
matched potentials and the controls' actual optimized subspaces. The text
report records the evidence and hashes required for that decision. A finite
campaign can expose a new physical issue, so the final CAS choices will follow
the completed native evidence.
