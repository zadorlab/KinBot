# VRC orbital selections and publication, 24 October 2026

The five new channels below are selected together for one additive live revision
of the CH3SS case input. They are ready for production optimization, subject to
the saved native files passing the publication checks. This is completion of
orbital selection, not certification of an entire relaxed potential or rate.
The user has requested retirement of the separate relaxation campaigns. Their unfinished,
timed-out and failed structures remain diagnostic evidence; none is promoted
as an accepted geometry. Production relaxed scans retain their connectivity,
root/orbital and correction-readiness review gates.

All active counts below mean **CAS(electrons, spatial orbitals)**. Source indices
refer to each named channel's own 2 angstrom natural-orbital record 2140.2, with
separate aDZ-F12 and aTZ-F12 checkpoints. Equal index ranges across channels do
not mean the same orbitals. The two original production spaces use their already
reordered/confirmed checkpoints. No source ordinal is a universal orbital label
at a later geometry or another basis.

D and Q denote doublet and quartet roots included in orbital optimization;
weights are equal within each listed multiplicity. The reactive-flux column is
separate. Averaging additional roots does not add their flux. All current and
new entries use direct CASPT2-F12 rather than a CC energy anchor or a spin-gap
correction. New entries retain the existing shift 0.1 Eh, IPEA 0.25, aDZ/aTZ
extrapolation exponent 3.7 and production geometry/scan settings.

| Status and channel | Space / source indices | Chemical contents | Averaging; reactive flux |
|---|---|---|---|
| Installed: initial CH3S C–S → CH3 + S | (5,4), 11–14 | Complete S 3p manifold: 4e/3o; methyl radical: 1e/1o | 3D + 3Q; D1 and D2, weight 1 each |
| Installed: CH2SH C–S → CH2 + SH | (5,4), 11–14 | SH pi frontier pair: 3e/2o; CH2 frontier pair: 2e/2o | 4D + 2Q; D1 |
| Installed: CH + H2S channel | (5,4), 8,12,13,14 | Sulfur p donor: 2e/1o; CH sigma/frontier subspace: 3e/3o | 4D; D1 |
| Previously installed, reaction rejected: C + SH3 | (3,4), 12–15 | SH3 frontier orbital: 1e/1o; complete C 2p set: 2e/3o | 3D + 3Q; formerly D1, now excluded |
| Installed: initial CH3SS C–S → CH3 + S2 | (7,5), ordered 18–22 | S2 pi/pi-star manifold: 6e/4o; methyl radical: 1e/1o | 1D + 1Q; D1 |
| Installed: initial CH3SS S–S → CH3S + S | (7,5), 18–22 in a distinct seed | CH3S sulfur frontier pair: 3e/2o; atomic S 3p manifold: 4e/3o | 6D + 6Q; D1 |
| New: CH2SHS C–S → CH2 + HSS | (5,4), 19–22 | HSS frontier pair: 3e/2o; CH2 frontier pair: 2e/2o | 2D + 2Q; D1 |
| New: CH2SHS S–S → CH2SH + S | (5,4), 19–22 in a distinct seed | CH2SH frontier orbital: 1e/1o; atomic S 3p manifold: 4e/3o | 3D + 3Q; D1 |
| New: CH2(SH)S C–S → CH2S + SH, homolysis approach | (7,6), 18–23 | Reorganizing cleavage/CH2S pi pair: 2e/2o; SH frontier pair: 3e/2o; retained C–S sigma pair: 2e/2o | 2D; D1 |
| New: HCS2 + H2 product exit | (7,6), **15,19,20,21,22,25** | HCS2 pi/pi-star pair: 2e/2o; sulfur-rich frontier pair: 3e/2o; H2 sigma/sigma-star: 2e/2o | 2D; D1 |
| New: CH2S + SH product exit | (5,4), 19–22 | CH2S pi/pi-star: 2e/2o; SH pi frontier pair: 3e/2o | 2D; D1 |

**Previously installed selections.** The original CH3S space closes C–H bonds
and sulfur 3s, retaining all three S p directions for orientation sampling.
Confirmed far occupations are 1.33333 three times and 1.00000. Adding S 3s
changed the tested doublet angular residual by 0.00039 kcal/mol; the larger
(13,8) comparison gave 0.10803. These are sampled relative-potential checks,
not universal error bounds. See the retained
[CH3S review](../../ch3s_anl1_f12/orbital_reviews/ACTIVE_SPACE.md).

For original CH3SS C–S, the larger-space selection was **17,18,20,21,22**;
17↔19 was already performed in the published checkpoint to produce 18–22.
Do not repeat that swap. The excluded S2 sigma/sigma-star pair had occupations
about 1.946/0.061; its (9,7) expansion changed the tested doublet angular
residual by 0.00239 kcal/mol. Original CH3SS S–S needs a different five-orbital
space because the bound CH3S sulfur has two frontier components, whereas the
new CH2SH + S channel retains only one bound-fragment frontier orbital.
The original S–S C–S sigma-pair expansion produced a 0.04494 kcal/mol angular
residual, with a larger 0.114 value at the tested inner point. See the
[CH3SS review](../../ch3ss_anl1_f12/orbital_reviews/ACTIVE_SPACE.md).

CH2SH uses SH's three-electron pi pair and CH2's two frontier electrons.
The larger (7,6) comparison left its additional **S–H** sigma pair near
1.98/0.02, supporting exclusion from the selected space. That pair is not
the C–S sigma pair added for the new CH2(SH)S channel. The earlier CH2SH review
did not complete a matched two-basis correlated-potential/angular validation.

For CH + H2S, the sulfur donor is source 8 at 2 angstrom but source 7 in the
far diagnostic. Selection therefore uses its checkpoint identity, not a common
energy-index window. Source 15 is S–H antibonding (~0.086 near), not the missing
donor. Its exclusion is a compact-model choice, without demonstrated correlated
potential convergence. Only four doublets remain in this model. The failed
auxiliary quartet/reference checks do not establish a spin-gap correction.

C + SH3 kept carbon 2s inactive and all three carbon p directions active.
Its subsequent relaxed scan changed S–H connectivity, so the reaction remains
rejected. The model's presence in input is historical provenance, not permission
to include that channel in rates. This publication does not reopen it.

**The two earlier provisional CH2SHS screens.** These belong to parent
`791671580700000000002-s36f1e2b924862f12`, not the initial CH3SS well.
C–S retains HSS frontier occupations near 1.5095/1.5050 and CH2's 1/1 pair;
the S–S sigma pair (~1.958/0.049) remains inactive. S–S retains the three
atomic S p occupations near 1.3333 and one CH2SH frontier orbital near 1;
the retained C–S sigma pair (~1.977/0.027) remains inactive. The small change
from the provisional geometry supports reuse as an initialization, but it does
not prove ordinal or state identity after relaxation. Binding checks the actual
request's atom order, products/spins and connectivity and records the geometry
difference. No relaxed-probe completion is claimed for either model.

**New CH2(SH)S homolysis comparison.** All 24 fixed-point comparisons completed.
The following is V(7,6) − V(5,4), where each D1 potential uses its own 30 angstrom
energy as zero; units are kcal/mol. Values are calculated from the user-reported
native energies in [compact_comparison_evidence.json](compact_comparison_evidence.json).

| C–S distance / angstrom | aDZ-F12 | aTZ-F12 |
|---|---:|---:|
| 2.0 | +1.0982 | +1.1140 |
| 2.75 | +0.4923 | +0.4886 |
| 3.5 | +0.0473 | +0.0538 |
| 6.0 | −0.0033 | −0.0019 |

Although the extra 18/23 C(1)–S(3) sigma pair stays nearly filled/empty, this
is a reproducible change in the inner potential shape in both bases. It
motivates choosing (7,6), rather than treating the absolute energy lowering as
validation. Each compact model's 30→42 angstrom change is below 0.002 kcal/mol.
The 19/22 and 20/21 character reorganizes between 2 and 3.5 angstrom in the
broad screens; smooth occupations and near-degenerate far roots alone do not
prove root or active-span continuity. Production must retain that review.

**Product-exit selections.** Both screens use the saved L2 product complexes,
not a stretched H–H distance inside the initial well. For HCS2 + H2, [6,4]
separates fragments [1,2,3,6] and [4,5]. Source 15/25 is H2 sigma/sigma-star
(~1.976/0.024); 19/22 is the delocalized HCS2 pi pair (~1.917/0.139); 20/21 is
the sulfur-rich three-electron pair (~1.542/1.450). Including H2's complete
pair permits its internal bond to respond in the active space. This is a
chemical inclusion, not a measured energetic advantage over (5,4). Other
sigma pairs remain inactive. The same assignments are present in both bases
at all four separations. No compact CASSCF/CASPT2 or angular comparison was
performed for this exit.

For CH2S + SH, [5,4] separates [1,3,5,6] from [2,4]. Source 19/22 is CH2S
pi/pi-star (~1.896/0.110 near; ~1.893/0.113 far), and 20/21 is the SH pair
(~1.500 each). The SH orbitals rotate/exchange order at large separation;
retain the pair as a span. C–S and S–H sigma pairs remain inactive (~1.97/0.03),
though their occupied and virtual members mix at short range. This is a
separate approach from homolysis, with its own seed and no matched compact
correlated-potential expansion test. Its (5,4) choice is therefore provisional
with respect to energetic convergence, not validated by the homolysis (7,6)
results. A common-product rate comparison will require a consistent energy
reference and scrutiny of the difference between these channel models.

AO coefficients used in these assignments are thresholded amplitudes, not
atomic populations; signs alone have no invariant meaning. The active-electron
counts follow the retained chemical manifolds. Occupation sums from a larger
space need not equal those integer counts exactly.

**Handoff.** First retire only the standalone review owners and jobs with
`queue_progress_review.py PRODUCTION --stop`. It waits for submission
bookkeeping, verifies owner identities, cancels the user's jobs whose WorkDir is
inside the review root, and verifies cancellation. It preserves all native
files, never cancels production jobs, restores the shared production cap to
23 and reserve to 0, and reconciles holds made by KinBot's pool manager.
A scheduler error remains an error; it is not reported as an empty queue.

Then run the generic `kinbot.requested_orbital_recon review-batch` with the
`all_selections.json` plan, campaign root `~/kinbot-runs`, and `--require-all`.
Preview before adding `--apply`. The plan uses five independent, basis-specific
screens and hash-checks the 24-job compact comparison as supporting evidence.
Missing/changed evidence blocks publication. One additive input revision
contains all five models; no review job or production controller is launched.
The existing PES owner consumes the revision at its safe gate. A queued receipt
means publication is pending, not that scans or rates have completed.

For this code update, the example `publish_progress_review.py` wraps that generic
publication with the authorized retained-driver refresh. It preflights the entire
batch before replacement, refreshes both production driver trees through
`refresh_retained_run.sh`, then publishes and exports evidence. Existing native
allocations and independent RotdPy owners are retained. Future sampling workers
load the new immutable code snapshot through their shared setup; already running
workers finish on their original code. Same-code repeated refreshes use the saved
journal. See [README.md](README.md) for the single handoff command.

The generic engines contain no molecular IDs or example-specific recovery
branches. This directory contains the scientific choices; the legacy example
launcher only identifies its historical review owners and selects the requested
capacity. No solver thresholds, scientific retry settings or accepted production
results are changed by this handoff.

**Correction scans versus sampling surfaces.** The two product exits explicitly
use `coordinate: com` in their surface models. RotdPy varies the separation of
the two fragment centers of mass. Their relaxed/frozen one-dimensional correction
still uses KinBot's mapped interfragment atom anchors, starting from the saved
product complex. It never stretches an H–H bond within the original well. The
asymptotic-fragment optimization and subsequent constrained relaxation follow
the normal production driver; no unfinished diagnostic geometry is imported.

Correction distances come from `vrc_tst_scan_points`, with the separately
configured asymptote (30 angstrom here). Sampling distances come from
`rotdpy_dist`. The publication preserves both lists and archives their effective
values, anchor atoms, sampling coordinate and source-geometry proof in
`orbital_reviews/live_refresh/*/review_summary.json`. Equal numerical grid values
do not make anchor and COM distances interchangeable. RotdPy evaluates the 1D
correction at the recorded atom-pair distance in each sampled configuration.
The small COM surfaces may reject overlapping fragments; the final variational
minimum must be checked to lie inside the usable surface grid.

For later review, preserve each reaction's `vrctst/scan_progress` and
`vrctst/scan_review` records, the native point graphs, relaxed/frozen geometries,
`corr_*.json` and scan-provenance sidecars, and the generated RotdPy input,
validation metadata and surface results. Check the actual atom maps and distances,
intrafragment connectivity, active-orbital span/root character in both bases,
asymptotic behavior, correction interpolation coverage and sampling convergence.
The scan tests exercise product-side geometry, a nondefault distance grid, the
separate asymptote and evidence records; MESS tests retain the inner saddle and
the outer sampled barrier, including required-flux failure behavior. These are
local integration tests, not completed native production scans.

**ANL and VPT2 for product complexes.** Complexes retain the requested minimum
ANL calculations, including inherited VPT2 jobs and their validation. There is
no covalent CBH fragmentation of a disconnected complex. Its formation energy is
constructed as

```text
Hf(complex, 0 K) = sum Hf(products, 0 K)
                + E0_harmonic(complex) - sum E0_harmonic(products).
```

The binding difference uses the best common accepted electronic expression and
harmonic ZPE on both sides. Product formation energies can retain their validated
VPT2 contributions. VPT2 files/profiles are preserved for the complex as well,
but their anharmonic corrections are not used in this binding subtraction. The
handoff explicitly references the harmonic complex profile, so MESS receives
that choice consistently. Native normal termination alone does not establish
that VPT2 treats floppy intermolecular modes accurately. Existing checks reject
incomplete or nonphysical results; retained results remain available for a
mode-by-mode review. VPT2 can be useful for some molecular dimers, but large
amplitude motions need particular scrutiny; see the primary
[V30 dimer study](https://arxiv.org/abs/2209.04392). No specific inherited native
VPT2 result is newly certified by this code change.

The update also fixes profile assembly for disconnected complexes: electron
counts come from the native atoms and charge, not the CBH connectivity routine
that intentionally rejects disconnected fragmentation targets. Neither fragment
CBH thermochemistry nor their VPT2 task definitions are disabled.

**MESS uses accepted L3 thermal data.** The handoff selects the highest completed,
verified L3 geometry and harmonic task independently. The main RRHO geometry,
rotational constants and external symmetry use that L3 geometry. The vibrational
spectrum uses the accepted L3 harmonic result. Original L2 geometry, Hessian,
conformer and hindered-rotor records remain intact. Only the selected conformer
that supplied the ANL graph inherits its L3 properties; uncalculated conformers
retain their own data. Atomic species have no vibrational modes.

For a stable species with an accepted VPT2 handoff, the frequencies are

```text
nu_i(MESS) = nu_i(L3, retained after internal-rotation projection)
           + [nu_j(L2, VPT2 fundamental) - nu_j(L2, harmonic)].
```

The mapping is a global maximum-overlap assignment of mass-weighted L3 and L2
normal vectors after atom mapping and proper rotational alignment. It is not
frequency sorting. The default requires squared overlap at least 0.5 and rejects
competing assignments within 0.1 overlap when their shifts differ by more than
1 cm-1. A reviewed one-based mode map can resolve mixing, but still must pass the
harmonic mismatch bound (100 cm-1 by default) and positive corrected-frequency
checks. Missing vectors, changed hashes, wrong geometries and ambiguous mapping
block the thermal handoff; they do not silently substitute unverified L3 data.

Molpro's full Cartesian normal-mode table supplies a vibrational spectral
Hessian to its printed precision. KinBot projects the same validated L2 rotor
identities from this L3 Hessian and retains mode vectors for matching. Invalid
or demoted rotors return their vibrational mode to the spectrum. Each MESS rotor
block contains the L2 geometry associated with its L2 scan, so its inertia does
not accidentally use the main L3 geometry. This is a complete eigenbasis
reconstruction, not an unprojected frequency list with the lowest modes deleted.
See the native [Molpro frequency documentation](https://www.molpro.net/manual/doku.php?id=harmonic_vibrational_frequencies_frequencies)
and the rotor geometry keywords in the
[MESS manual](https://tcg.cse.anl.gov/papr/codes/mess/messmanual_3_23_16.pdf).

Without accepted L3 thermal data, the existing L2 model remains the fallback;
any L2-only VPT2 correspondence follows the legacy L2 matching policy. A partial
ANL graph can upgrade geometry without inventing an unavailable Hessian. TS
entries use L3 real modes and the L3 imaginary frequency for tunneling, without
automatic VPT2 shifts. Ordinary products retain their own accepted VPT2 shifts.

Complex wells retain their L3 intermolecular harmonic modes and receive no
automatic complex VPT2 shifts. Their outer RotdPy barrier contains the products'
internal modes and sampled relative-motion flux, not the complex's soft harmonic
modes a second time. This is a harmonic model of the bound complex, not proof
that a floppy well is quantitatively converged or a newly fitted large-amplitude
intermolecular potential. Final complex retention is reassessed using accepted
ANL/formation energies and its L3 harmonic stability when available.

The complete ANL energy expression separately uses its accepted harmonic ZPE
and L2 VPT2 ZPE correction; the L3 profile uses L3 harmonic ZPE, whereas fallback
profiles can use L2 harmonic ZPE. Frequency changes never add a second ZPE to
that energy. External CBH MESSPF references use the same verified L3/shift policy
while retaining their L2 rotor data. The handoff and thermal-model provenance
record native hashes, atom maps, mode assignments, overlaps, geometry and rotor
sources. Local tests exercise native Molpro/Gaussian vector parsing, signed TS
modes, atom permutations/rotations, quality/hash gates, rotor demotion, and the
final MESS inner/complex/outer assembly. Synthetic VPT2 shifts used in those
tests are explicitly labeled; neither a native VPT2 calculation nor a production
MESS solver run is claimed by these tests.

**Why the non-homolysis requests exist.** Automatic VRC selection in
`kb._accepted_homolytic_channels` includes accepted two-product reactions with
`do_vdW=True`, in addition to homolytic scissions. The function name is narrower
than its actual scope. The supplied accepted L2 records for `h2_elim_4_5` and
`r12_insertion_R_2_1_5` both contain `product_complex`; `pes_snapshot.snapshot`
writes that field only when this flag is true. Thus these are product-complex
exits, not an attempt to replace either tight saddle with VRC.

With `me_skip_vdW=0` (the default), the MESS connections are:

```text
parent well -- inner tight TS -- product complex -- outer VRC -- fragments
```

`MESS.write_input` retains the ordinary saddle and product complex, and
`pes.create_mess_input` connects the sampled outer barrier to that complex.
For multiple routes to the same products, the assembler selects a shared
complex and avoids duplicate outer barriers; it does not simply add every
request as independent flux. The inner saddle's tunneling reference uses the
selected complex. Required/production ROTD mode refuses missing or invalid
sampled flux rather than silently substituting phase-space theory.
`me_skip_vdW=1` explicitly omits the complex and outer exit; that setting is
not changed by this publication. A final generated MESS file is still needed
to verify the particular live network after its outstanding QC completes.

The supplied live inputs have `vdW_detection=0.0`. Positive binding can therefore
qualify, but the last log extraction did not return numerical depths. KinBot
also reassesses the optimized complex's binding and frequencies before retaining
it. The saved `product_complex` entries establish the accepted classification;
this report does not invent a numerical depth or certify that every weak well
will remain in the final composite network.

**Reaction identities and coordinates.** Atom numbers below are one-based in
the retained request order, not numbers in an independently reordered fragment.

| Parent routing name | Reaction suffix | Breaking bond or outer coordinate |
|---|---|---|
| `470831050000000000002` | `hom_sci_1_2` | C1–S2, CH3 + S |
| `471030810040000000002` | `hom_sci_1_2` | C1–S2, CH2 + SH |
| `471230610040000000002` | `hom_sci_1_2` | C1–S2, CH + H2S |
| `471430450000000000002-s06a7b4c221beb748` | `hom_sci_1_2` | C1–S2, C + SH3; rejected |
| `791471490990000000002` | `hom_sci_1_2` | C1–S2, CH3 + S2 |
| `791471490990000000002` | `hom_sci_2_3` | S2–S3, CH3S + S |
| `791671580700000000002-s36f1e2b924862f12` | `hom_sci_1_2` | C1–S2, CH2 + HSS |
| `791671580700000000002-s36f1e2b924862f12` | `hom_sci_2_3` | S2–S3, CH2SH + S |
| `791472110370000000002` | `hom_sci_1_2` | C1–S2, CH2S + SH |
| `791472110370000000002` | `h2_elim_4_5` | Outer H6···H4 separation between HCS2 and H2; no covalent bond is broken by that coordinate |
| `791472110370000000002` | `r12_insertion_R_2_1_5` | Outer H5···H4 separation between CH2S and SH; no covalent bond is broken by that coordinate |

For the elimination's inner path, the endpoint geometries show H4 leaving S2
and H5 leaving C1 to form H2 (H4–H5); H6 stays on HCS2. For the insertion-family
path, the saved endpoints give net C1–S2 separation into CH2S and SH; the
family label alone does not establish the detailed sequence along the IRC.
The two outer H···H anchor distances are interfragment distances, not H2 bond
stretches. Their scan geometry is the saved product-side complex.

**Monitor reconciliation, snapshot 2026-10-24 14:14:48 −07:00.** The configured
production limit was 22 with a three-slot orbital-review reserve, hence the
effective limit was 19. Retirement sets 23 and reserve 0; existing production
node exclusions remain unchanged. This is a simultaneous-allocation cap, not
a guarantee that Slurm has 23 usable idle nodes.

The previously failed L3 CC task `ts_0450141d9a891e9c/l3_geometry` was running
as 52928790. Its bounded iteration recovery was implemented, but this snapshot
does not establish native convergence. Three CH3S channels had complete
48-point correction scans and four active RotdPy owners; each still reported
0/116 final surface files, so neither sampling nor rates were complete. The
two initial CH3SS scans were progressing (27 and 15 accepted points; running
at 7.4 and 5.0 angstrom), and its PES exploration remained active. There were
18 production-root allocations, including 13 submitted ANL tasks. Existing jobs
are retained. The five review holds are the five models in this publication;
the old branch reports are explicitly marked as belonging to earlier controllers.

The −3.62 kcal/mol L1 path `Intra_R_Add_Endocyclic_F_3_1` was rejected because
neither IRC returned to the initial well. Separate TS searches with no qualifying
imaginary mode are scientific validation rejections. Neither should be converted
into an iteration-limit retry. Old optimizer tracebacks and an isolated warning
without a terminal failure do not establish a new live failure class. The old
failed production diagnostic probe was inactive. The separate review root is
outside this monitor's scope and is inventoried explicitly during retirement.

**Evidence transport.** Use [EXPORT.md](EXPORT.md) and [export.json](export.json)
to preserve native outputs, WFU files, inputs, selection receipts and attempt
history. The user runs remote commands and SSH transfers. Local byte verification
after transfer is separate from the chemical review; new remote orbital files
are not described as locally archived until that transfer has succeeded.
