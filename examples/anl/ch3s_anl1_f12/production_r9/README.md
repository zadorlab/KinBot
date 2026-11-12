# Lean R9 review stage

Prepared on 2026-11-10 from the selected lean R8 models, including CH3 + S.
The native-review update of 2026-11-12 is recorded in
`../../ch3ss_anl1_f12/orbital_reviews/r9_completion_20261112/REVIEW.md`.
The current package is a review draft. Production submission remains pending the
specific scientific and result-import items below.

## Retained scope

| PES | Wells | Product sets | Connections | VRC connections |
|---|---:|---:|---:|---:|
| CH3S | 2 | 2 | 4 | 2 |
| CH3SS | 5 | 4 | 9 | 5 |

CH3S retains CH3S and CH2SH, the CH2S + H product set, and CH3 + S. Its tight
connections are the direct CH3S H-loss route and CH3S/CH2SH isomerization.

CH3SS retains CH3SS, CH2SSH, HSCH2S, the CH2S/SH encounter complex, and the
HCS/H2S encounter complex. It retains CH3 + triplet S2, CH2SS + H, CH2S + SH,
and HCS + H2S. Two isomerization saddles, the inner route into the CH2S/SH
complex, and the verified reverse abstraction connecting the two encounter
complexes complete the nine connections.

`scope.json` in each case records the exact native identities and endpoints.
`lean_scout.inp` preserves the tested L1/L2/PST baseline. Its PST blocks are the
connections selected for VRC refinement. The optional HSCH2S to HC(=S)SH + H
extension is outside this chosen lean scope.

## Every retained VRC model

Indices refer to the natural orbitals in record `2140.2` of the exact reviewed
WFU files. State-average roots support the orbitals. The reactive-root column
specifies the flux contributions.

| PES and connection | CAS and orbital content | Natural indices | State average / reactive roots | Current status |
|---|---|---|---|---|
| CH3S `nobar_1`, CH3S to CH3 + S | (5,4), methyl radical and all three sulfur 3p orbitals | 11–14 | Three doublets and three quartets / D1 and D2, original weight 1 each | Existing native-reviewed model ready |
| CH3S `nobar_4`, CH2SH to CH2S + H | Exact H-loss model remains to be selected | Pending | Pending | New exact review required |
| CH3SS `nobar_1`, CH3SS to CH3 + S2 | (7,5), methyl radical and full S2 pi/pi-star manifold | 18–22 | One doublet and one quartet / D1 | Existing native-reviewed model ready |
| CH3SS `nobar_2`, CH3SS to CH2SS + H | Exact H-loss model remains to be selected | Pending | Pending | New exact review required |
| CH3SS `nobar_14`, HSCH2S to CH2S + SH | (7,6), cleavage/reorganization pairs and both SH frontier components | 18–23 | Two equally weighted doublets / D1 | Angular checks verified, corrected seed installation pending |
| CH3SS `nobar_39`, CH2S/SH complex to CH2S + SH | (5,4), CH2S pi/pi-star and both SH frontier components | 19–22 | Two equally weighted doublets / D1 | Exact endpoint checks verified, production installation pending |
| CH3SS `nobar_41`, HCS/H2S complex to HCS + H2S | (3,3), HCS occupied pi, SOMO and pi-star | 20–22 | One doublet / D1 | Selection supported; production installation pending |

The first and third models retain their basis-specific records and settings from
the existing CH3S and CH3SS example inputs. Their native evidence is described in
the respective `orbital_reviews/ACTIVE_SPACE.md` files. Both bases' hashes are
included in `input.draft.json`.

The installed HSCH2S (7,6) selection is preserved from the R7 review. Its compact
native comparison found an own-asymptote potential difference of 1.0982 and
1.1140 kcal/mol at 2 angstrom when the extra pair was removed, in aDZ-F12 and
aTZ-F12 respectively. The 30-to-42 angstrom changes were below 0.002 kcal/mol.
Case 7's new radial and completed near-relaxation evidence supports retaining
the extra pair. The sixteen corrected angular checks now recover SH doublet
degeneracy to 0.000270 kcal/mol, with return differences below 0.000572 kcal/mol.
The new separated-start seeds must accompany the production model. The saved
earlier selection is in the CH3SS package's `prior_reviews/7.json`.

For the CH2S/SH outer connection, existing case 13 has the closest endpoint.
Its proper-map RMSD to the retained R8 L2 endpoint is 0.0005295 angstrom. Matching
the atom order exchanges the last two hydrogens. The geometric mapping is saved
in `endpoint_mapping.json`. All eight exact endpoint checks have now passed,
including return energy differences below 0.000032 kcal/mol and native active
subspace return overlaps above 0.999999999965. Production seed installation is pending.
Case 13 now supplies completed radial, angular, return and relaxed-midpoint
evidence. Its tested occupied-union potential differences are at most 0.0051
kcal/mol. The auxiliary case 8 angular results have electronic-branch anomalies
and are recorded separately in the new review.

For the HCS/H2S outer connection, case 25 supplies the selected (3,3) model.
The old fixed-orbital H2S-donor control gave energy changes of approximately
1.7e-7 and 1.5e-7 hartree. The submitted optimized donor comparison is needed to
test that proposed enlargement. The new optimized union expels the H2S donor
again, while the actual HCS subspace remains stable. The reviewed evidence
supports retaining (3,3) for the outer HCS/H2S channel. Forward and reverse discovery names share this endpoint
and review. The current L2 endpoint has the same atom order and a 0.03327
angstrom proper-map RMSD to the original review geometry.

The two exact H-loss requests are saved in `new_h_loss_reviews.json`, with native
L1/L2 parent coordinates, product multiplicities, bond indices and source hashes.
Their active spaces remain unassigned. Existing C-S cleavage reviews cover
different orbital reorganization.

## Reuse and submission requirements

The local prepared copy is under
`local_runs/r9_lean_stage_20261110/staged`. It contains six SHA-verified installed
WFU seeds, the prior candidate review records, and frozen R8 source snapshots
for seven workers. Those snapshots contain 4,261 normally completed database
rows and 115 terminal rejected-search records. These counts cover the source
workers, including their unselected reactions. `reuse_inventory.json` and
`source_snapshot_manifest.json` record their provenance.

Remaining work before production submission comprises the following items.

- Install the reviewed case-7 separated-start and case-13 exact-endpoint seeds,
  together with the supported case-25 model. The sixteen angular and eight
  endpoint checks are complete and reviewed in `initial_native_review.json`.
- Complete the two exact H-loss orbital reviews. Their initial native reference
  panel is complete. Compact-space screening resumes from those saved results.
- Implement and validate independent-run result adoption. The frozen copies
  currently serve as evidence. Runtime adoption still needs path rebinding,
  geometry/theory/state validation, old scheduler-ID reconciliation, and
  preservation of terminal search outcomes. Existing cross-worker reuse and
  same-run refresh routines do not complete this independent-run import.
- Integrate the retained reverse saddle and the HCS/H2S encounter complex into
  the lean production ANL/VRC network. Their source covalent worker is outside
  the selected exploration scope. This handoff must retain the selected endpoint
  and saddle while respecting the chosen five-population network.

Previously completed R7 ANL/VRC calculations can be adopted once their complete
result contracts match the selected geometry, theory and surface model. The WFU
seeds have already been verified and copied. Adoption of other native results
remains pending the independent-run importer.

## Settings and local checks

The drafts inherit the requested ANL1-F12 augmented thermochemistry, direct
CASPT2-F12 aDZ/aTZ CBS correction, reviewed spin/root definitions, electronic
levels and spin-orbit centroid declarations. They retain the R8 150 kcal/mol L1
ceiling, matched-ZPE 0.5 kcal/mol vdW test, existing L1/L2 theory profiles,
conformers and hindered rotors. The shared-pool planning settings remain 20
jobs, eight processors each, exclusive nodes, and the latest R8 node exclusions.
Temperatures span 300–2000 K with the existing six-pressure list. MESS writing
is enabled and automatic MESS execution is disabled for the eventual review.

The following local checks completed with external calculations disabled.

- Full parameter resolution for both drafts and default-removal equivalence.
- Normal reaction discovery for all five retained covalent wells. Every selected
  motif was found. Five automatic dissociation-saddle companion motifs were
  also generated by the existing discovery logic.
- Production Molpro rendering for all three installed channel models, giving
  twelve sample/correction/root/basis components with verified WFU hashes.
- Hash-verified source copying, preserved rejection histories, and repeated
  stage construction with matching source bytes.

`validation.json` contains the detailed discovery and rendering records.
These checks establish input and evidence consistency. A complete independent
R9 import, controller restart and native production execution remain unverified.

The repository drafts retain original seed locations. The local prepared copies
use verified relative `orbitals/` paths. Use the completed import/staging handoff
when it is available to produce the final submission input.
