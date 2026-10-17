# Electronic levels for the saved sulfur PES runs

The 17 October user-supplied `network_request.json` inventory covers both PES
and their CBH references. All entries are neutral. The main example JSONs
contain the level data and sources under `electronic_state_models`; no second
configuration file is required.

## Where the data enter

`levels` supplies energies above the lowest level in cm-1 and total electronic
degeneracies. KinBot writes these into the existing MESS/MESSPF
`ElectronicLevels[1/cm]` grammar. Thus they affect electronic partition
functions, thermal thermochemistry and statistical kinetics. MESSPF can run
before the final network MESS calculation; this is not exclusively a final
output decoration. Molecular rotational degeneracy belongs to the rotor and
is not multiplied into these electronic degeneracies again.

Other fields in the same model have separate meanings: `ground_multiplicity`
sets the QC spin, and `term_orbital_degeneracy` can determine spin-free spatial
root counts when they are not explicitly supplied. The three production CAS
manifolds remain 3/3, 1/1 and 6/6 doublet/quartet roots. Experimental SO
splittings do not become CAS roots or modify the scalar-relativistic VRC PES.

Atomic S already explicitly requests `spin_free_reference: term_centroid`.
That separate option adds the known centroid-to-lowest-J correction to ANL
electronic energies. Its value and convention are unchanged. Adding molecular
levels alone requests **no** such correction. In particular, half the CH3S
spectroscopic splitting cannot automatically be subtracted from an optimized,
Jahn-Teller-distorted spin-free minimum: the reference energies differ.

For S the correction is
`-(5*0 + 3*396.055 + 1*573.640)/(5+3+1) = -195.7561 cm-1`,
approximately -0.560 kcal/mol. MESS does not infer this subtraction from
`ElectronicLevels`: its [manual](https://tcg.cse.anl.gov/papr/codes/mess/messmanual_3_23_16.pdf)
defines those as excitation energies relative to the ground level, separately
from the species `ElectronicEnergy`/`ZeroEnergy`. KinBot's ANL component assembly
adds an explicit SO correction only if it is absent from the assembled component
list. A user supplying already SO-corrected electronic energies must not also
request this centroid correction.

This input hook is a **composite-branch addition**, present in commit
`4793a0f6` (14 October 2026), not an origin/master feature and not newly
introduced by this restart. Master writes the historical single ground
electronic level with spin multiplicity. The branch's ANL recipe declares
`scalar_relativistic` and `spin_orbit` as separate components; the former is
the spin-free DKH-minus-nonrelativistic difference. DKH2 in these decks does
not split a 3P term into its J levels. Separate terms follow the literature:
[ANL0-F12 application with explicit DKH and experimental SO terms](https://www.osti.gov/servlets/purl/1878167)
and [ANL0F' supporting information with DKH and atomic Cl SO terms](https://www.osti.gov/servlets/purl/2311339).
The new molecular level data in this restart do not introduce another ANL SO
correction. Their energy-zero treatment therefore retains the existing
spin-free molecular approximation; level populations alone cannot correct
that approximation.

## Species review and input decisions

| Inventory identity | Treatment |
| --- | --- |
| `[S]`, triplet atomic sulfur | Existing levels (cm-1, g): (0,5), (396.055,3), (573.640,1); existing centroid correction retained. |
| `C[S]` = `[CH3][S]`, CH3S | Add X 2E3/2 (0,2) and X 2E1/2 (255,2). Gas-phase IR/VUV photoionization reports 255 +/- 4 cm-1. |
| `[SH]`, SH | Add X 2Pi3/2 (0,2) and X 2Pi1/2 (376.835,2), from the experimental ground-state spin-orbit constant A0=-376.835 cm-1. |
| `CS[S]`, CH3SS | Retain ground doublet g=2. Its nondegenerate ground electronic state is not a sulfur-atom 3P manifold. No experimental neutral fine-structure doublet is assigned here. |
| `[S]S` / `S=[SH]`, HSS | Retain ground doublet g=2. These two labels have the same explicit-H connectivity. The known A-X electronic band is not an SO pair to insert as atomic fine structure. |
| `[S][S]`, S2 | Retain X 3Sigma-g ground term g=3. Its rotational spin fine structure is not represented as independent atomic 2J+1 levels. Higher singlet electronic terms are discussed below. |
| `C`, `S`, `CS`, `C=S`, `[H][H]` | CH4, H2S, CH3SH, CH2S and H2: retain ground singlet g=1. Plain `C` and `S` here are hydrogen-saturated SMILES, not atoms. |
| `[CH3]`, `[H]` | Retain ground doublet g=2; no additional low-lying SO components supplied. |
| `C=[SH]`, `C=S=[SH]`, `C=[SH][S]`, `[S]S=C` | Retain the discovered doublet g=2. No experimentally sourced low-lying level set for these particular isomers was established in this review. Do not transfer CH3S or atomic-S levels merely by elemental composition. |

Sources for the additions:

- [Lee, Han and Fu, IR spectrum of mass-selected CH3S (2011)](https://www.sciencedirect.com/science/article/pii/S0009261411008153).
- [Milan et al., Chemical Physics Letters 239, 326-331 (1995)](https://doi.org/10.1016/0009-2614(95)00488-2),
  [author-repository full text](https://pure.uva.nl/ws/files/2012241/28655_ChPhL_1995_239_326.pdf).
  A photodissociation population ratio reported in this paper is not an
  equilibrium degeneracy and is not used as one.

S2 also has measured excited singlet terms near 0.545 and 0.989 eV; these are
electronic excitations, not splitting of its ground triplet. They are not added
as an unqualified SO repair. Using identical rotational/vibrational partition
functions, their electronic contribution relative to the ground term would
be about 2.8% and 0.1%, respectively, at 2000 K, and much smaller at 300 K.
This is a ground-term thermochemistry approximation, not a claim of an exact
S2 partition function. See [Sun et al., JPCA 123 (2019), 6886-6896](https://pmc.ncbi.nlm.nih.gov/articles/PMC6699204/),
which cites the measurements. General JSON controls already allow these
levels; a state-specific rovibrational treatment is a separate refinement.

## CH3S negative curvature and the limits of this model

CH3S's X 2E state has Jahn-Teller coupling; published spectroscopy describes
it as relatively weak. It can make a symmetric spin-free stationary point
unstable to distortion. This is a plausible contributor to the earlier
negative mode, not a demonstrated diagnosis of that particular B2PLYP result.
The user's independent finite displacements measured negative curvature from
both energies and forces (approximately -1.58 and -1.54 eV/A2 at +/-0.01 A).
That minimum therefore had to be rejected and reoptimized, irrespective of
its attribution. The existing general displacement/reoptimization and fresh
Hessian checks remain enabled. No switch to restricted-open-shell B2PLYP or
relaxation of frequency acceptance is made.

[Bouallagui et al., PCCP 19 (2017), 31245-31254](https://pmc.ncbi.nlm.nih.gov/articles/PMC5995305/)
discuss the ground term and JT/SO coupling. The two measured CH3S origins used
here are an RRHO electronic-level approximation. They do not implement a full
coupled spin-vibronic partition function or cure a bad geometry. Rate accuracy
and thermochemistry still need the completed workflow's scientific review.

## Restart behavior

`resume_on_blodgett.sh` copies the reviewed level specifications into the
saved main/worker/reference inputs after verifying that existing spins,
spatial term degeneracies and energy-zero conventions agree. Everything else
in those JSONs is retained. Original inputs and preparation hashes are saved
outside the run. MESSPF rendering applies current input levels to cached
reference models without changing accepted QC data. Tests verify this behavior
and the unchanged three CAS root manifolds. New clean runs use the same main
example JSONs directly.
