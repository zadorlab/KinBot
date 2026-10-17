# Electronic levels and spin–orbit corrections for the sulfur PES runs

The 17 October user-supplied `network_request.json` inventory covers both PES
and their CBH references. All entries are neutral. Each main example JSON
contains the data, reference convention and sources under
`electronic_state_models`; no second configuration file is required.

## Approved additive model

The user approved the following practical approximation on 17 October:
retain the scalar ANL recipe, including its scalar-relativistic component and
harmonic ZPE plus anharmonic correction; add one experimental ground-level SO
correction. There is no separate Jahn–Teller energy correction. This supersedes
the earlier levels-only molecular treatment in this document and the inputs.

| Species | Excitation levels (cm-1, electronic degeneracy) | ANL SO correction (cm-1) | kcal/mol |
| --- | --- | ---: | ---: |
| S, 3P | (0,5), (396.055,3), (573.640,1) | -195.756111 | -0.559695 |
| SH, 2Pi | (0,2), (376.835,2) | -188.4175 | -0.538713 |
| CH3S, effective 2E manifold | (0,2), (255,2) | -127.5 | -0.364541 |

All three inputs explicitly request `spin_free_reference: "term_centroid"`.
The general implementation computes

```
delta_SO = -sum(g_i * excitation_i) / sum(g_i)
E0 = E_ANL_scalar + delta_SO + ZPE_harmonic + delta_ZPE_anharmonic
```

The same additive correction is included once in any selected L2/L3 fallback
expression. Its component, value, source and reference are retained in each
composite JSON (`components.spin_orbit`, `electronic_terms`, and
`spin_orbit_model`). The raw QC energies, native output hashes and ZPE
components are not rewritten. Corrected zero-K energies enter the CBH reaction
balance; resulting heats of formation then receive the existing MESSPF thermal
increments. Experimental ATcT formation enthalpies are used as supplied, not
SO-corrected a second time.

`levels` separately supplies excitations relative to the lowest component for
MESS/MESSPF `ElectronicLevels[1/cm]`. They change electronic populations and
thermal increments; they do not shift the zero-K baseline again. Molecular
rotational degeneracy belongs to the rotor, not to these electronic weights.
The [MESS manual](https://tcg.cse.anl.gov/papr/codes/mess/messmanual_3_23_16.pdf)
defines the level energies separately from `ElectronicEnergy`/`ZeroEnergy`.

This is the existing **composite-branch** ANL `spin_orbit` component and JSON
hook (introduced in commit `4793a0f6`, 14 October), not an origin/master feature.
The ANL `scalar_relativistic` DKH difference is a separate spin-free component;
it does not split the sulfur 3P term. Separate scalar-relativistic and SO terms
are used in the [ANL0-F12 application](https://www.osti.gov/servlets/purl/1878167)
and [ANL0F' supporting information](https://www.osti.gov/servlets/purl/2311339).

The interface is general: matching uses SMILES/declared aliases and charge,
with multiplicity validation. A complete supplied term can use the centroid
option for an atom or a molecule. A user can instead supply
`spin_orbit_correction_hartree` for a different reference treatment; both
options together are rejected. A levels-only model requests no energy shift.
There is no sulfur-specific correction formula in KinBot code. An OH 2Pi model
uses the same interface when its measured gap is supplied; OH is absent from
these H/C/S networks and is not added to their inputs.

## Experimental sources and approximation limits

- S: [NIST S I levels](https://physics.nist.gov/PhysRefData/Handbook/Tables/sulfurtable5_a.htm).
- CH3S: [Lee, Han and Fu, gas-phase IR/VUV photoionization (2011)](https://www.sciencedirect.com/science/article/pii/S0009261411008153),
  255 +/- 4 cm-1 between the two observed components.
- SH: ground-state A0=-376.835 cm-1 tabulated by
  [Milan et al., CPL 239 (1995) 326–331](https://doi.org/10.1016/0009-2614(95)00488-2),
  [author-repository text](https://pure.uva.nl/ws/files/2012241/28655_ChPhL_1995_239_326.pdf).
  Their photodissociation population ratio is not an equilibrium degeneracy.
  The separate SO/scalar-relativistic thermochemical treatment of SH is also
  discussed by [Csaszar, Leininger and Burcat (2003)](https://doi.org/10.1021/jp026605u).

SH uses the standard first-order equal-component centroid approximation.
For CH3S, applying half the *observed effective* gap to the scalar optimized
minimum plus scalar ZPE is an explicitly accepted approximation, not an exact
reference conversion derived from spectroscopy. The gap alone does not specify
an absolute energy shift from that scalar minimum. Its 4 cm-1 experimental
uncertainty is not the uncertainty of the whole thermochemical model; no
sub-kcal/mol accuracy guarantee follows solely from this prescription.

[Marenich and Boggs (2004)](https://doi.org/10.1021/jp040443f) distinguish the
spin-free JT stabilization (93 cm-1), pseudorotation barrier (15 cm-1), and
calculated electronic SO splitting (358 cm-1). Their
[coupled spin-vibronic study (2005)](https://doi.org/10.1021/ct0501452) solves a
Hamiltonian containing JT coupling, SO and vibrational motion. The experimental
255 cm-1 gap reflects their combined effect; it is neither the bare 358 cm-1
constant nor an instruction to subtract the JT energy from it arithmetically.

Ordinary Gaussian VPT2 treats local anharmonic vibrations on one electronic
surface. An optimized distorted geometry already includes static JT lowering
on that scalar surface, but VPT2 does not solve coupled JT/SO dynamics. We keep
the accepted scalar geometry and its validated anharmonic correction, then
apply the chosen effective SO shift. We add neither an extra -93 cm-1 JT term
nor a second correction using 358 cm-1. This is an RRHO/experimental-level
model, not a full spin-vibronic partition function.

JT coupling is a plausible contributor to the earlier CH3S negative mode,
not a proven attribution. Independent finite displacements gave negative
curvature from both energies and forces (about -1.58 and -1.54 eV/A2 at
+/-0.01 A), so that stationary point required reoptimization regardless.
Existing minimum recovery and fresh Hessian acceptance remain enabled.

## Complete supplied inventory and present decisions

| Inventory identity | Treatment at this time |
| --- | --- |
| `[S]`, `[SH]`, `C[S]` = `[CH3][S]` | The three corrections and level sets above. |
| `CS[S]`, CH3SS | Ground doublet only. The observed neutral ground state is nondegenerate X 2A-double-prime; do not assign an atomic 3P or CH3S 2E pair. |
| `[S]S` / `S=[SH]`, HSS | Ground doublet only; these labels have the same explicit-H connectivity. Its A-X electronic band is not a ground-term SO pair. |
| `[S][S]`, S2 | Ground X 3Sigma-g-minus, g=3. Rotational spin fine structure is not an independent atomic-J excitation ladder. No first-order orbital SO centroid correction is assigned. |
| `C`, `S`, `CS`, `C=S`, `[H][H]` | CH4, H2S, CH3SH, CH2S, H2: ground singlets. Plain C and S here are hydrogen-saturated SMILES, not atoms. |
| `[CH3]`, `[H]` | Ground doublets with no additional orbital SO pair. |
| `C=[SH]` | CH2SH connectivity: nondegenerate ground doublet, not the CH3S 2E term. |
| `C=S=[SH]`, `C=[SH][S]`, `[S]S=C` | Retain the discovered doublets; no established experimental ground-term SO multiplet correction for these specific isomers was identified. No correction is transferred by elemental composition. |

Thus **no additional experimental first-order ground-term SO correction is
added for the current inventory**. This does not assert that second-order SO,
spin-rotation or spin-spin effects are identically zero. New species outside
this inventory would require a new electronic-state review.

Supporting sources: [CH3SS photoelectron study](https://doi.org/10.1063/1.1343081),
[HSS rotational spectroscopy](https://pubmed.ncbi.nlm.nih.gov/10712872/),
[CH2SH photoelectron study](https://pubs.rsc.org/en/content/getauthorversionpdf/D4CP03906H),
and [S2 electronic-state calculations](https://doi.org/10.1021/acsearthspacechem.3c00116).
S2 has higher singlet terms near 0.545 and 0.989 eV; they are distinct electronic
excitations, not a missing first-order SO correction. They remain outside this
ground-term model. See [Sun et al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6699204/).

## Restart and VRC behavior

`resume_on_blodgett.sh` requires both saved PES runs inactive. It archives old
inputs, logs, failed attempts and derived thermochemistry/MESS artifacts in a
sibling directory, updates the reviewed electronic data, then resumes both
controllers. Completed QC, frozen network/QC source hashes and rejected
reaction decisions remain intact. Only newly requested SO postprocessing is
allowed to differ during frozen-plan reuse; geometry, raw energies, methods,
frequencies and spin changes still invalidate that reuse.

Automatic assembly consults current input SO data even for previously cached
CBH reference models. Stable-species corrections are not inherited by TSs
merely because a TS record uses reactant SMILES as a routing label. MESSPF
similarly reads current level data when rendering cached references. This
makes a clean run and a saved-results resume use the same electronic model.

`ground_multiplicity` and `term_orbital_degeneracy` retain their separate roles
in QC spin and spatial-root counting. The reviewed CAS manifolds remain 3/3,
1/1 and 6/6 doublet/quartet roots. None of these experimental corrections is
inserted into CASPT2 sampling, the spin-splitting correction potential, or an
invented distance-dependent SOC surface. MESS uses the corrected asymptotic
thermochemistry and the supplied level populations within the existing model.
