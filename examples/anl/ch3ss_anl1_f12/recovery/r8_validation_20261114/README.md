# R8 connection and literature validation — 14 November 2026

**16 November update:** use [the reconciliation handoff](RECONCILIATION.md) for
the existing panel. Native evidence corrects the earlier TS5 endpoint inference
and shows that several apparently completed IRCs stopped after one point.
The updated owner holds the unsent independent IRC grid while a finite diagnostic
wave runs. The original launcher below is retained as historical setup guidance.
The second-wave handoff repairs the rejected DVV time steps and separates the
TS2 CC and TS3 ROHF failures. Use `resume_probe_repairs.sh` for the existing
failed probes; native convergence remains to be verified.

This independent diagnostic campaign leaves the R8 PES controllers, accepted
results, rejection history and MESS files unchanged. It does not launch MESS
or promote a channel. It shares R8's existing 20-node pool, using at most
eight cores on an exclusive node for each calculation.

## Why this work, rather than repeating all of L1

The retained scouting model already uses L2 geometries/frequencies. Its
unverified part is the use of some **L1 IRC connections for L2 saddles**.
The direct-L2 export demonstrates that a saddle can move or connect different
species at L2. That requires connection checks, not indiscriminate repetition
of every unsuccessful L1 guess.

The matched database audit contains 27 CH3S and 211 CH3SS L1/L2 observations
(including repeated searches). Relative energy changes reach about 3.94 and
11.86 kcal/mol, respectively. These include geometry changes and are **not**
a method-error bound. Thus this panel also checks all 24 explicit L1 energy
exclusions recorded in the supplied R8 snapshot at L2. It does not claim that
all possible paths, or every failed L1 guess, have been exhausted.

## Remaining paper–R8 questions and the staged resolution

Paper labels refer to DOI [10.1039/d4cp04873c](https://doi.org/10.1039/d4cp04873c).
W4 = CH3SS; W2 = CH2SSH; W3 = HSCH2S; W1 = CH2S···SH.
Labels are proposed correspondences, never acceptance criteria.

| Question | Supplied direct-L2 observation | New checks |
|---|---|---|
| TS1, W1 to HCS + H2S | Paper-seeded saddle gives the abstraction; a separate native reverse seed has one unconverged IRC direction | Independent native Gaussian IRCs from both saddles, both directions and two step sizes; matched T/Q energies |
| TS2, W2 to W1 | Paper seed instead connects W2 conformers; native `rxn_54_1` connects W2 to complex A | Independent IRCs for both; SI-seeded coupled-cluster optimization, frequencies and both QSD paths; T/Q comparison |
| TS3, W3 to W1 | L2 paper seed connects W2 to complex B, not W3 | Independent L2 integration; SI-seeded coupled-cluster saddle/frequencies and both QSD paths |
| TS4, W4 to W2 | L2 endpoints match; native `rxn_2_1` duplicates this saddle | Independent L2 integration and T/Q energies, one representative; separately check `rxn_2_2`, which instead leads from W2 to CH2···HSS |
| TS5, W3 to HC(=S)SH + H | Geometry stays close to SI; the final optimized endpoint is CH(SH)2. The later native audit shows H separation on the raw IRC followed by recapture during optimization, so that endpoint does not establish a path disagreement | Independent native L2 integration and SI-seeded coupled-cluster saddle/frequencies and both QSD paths |
| TS6, W2 to W3 | Paper-like saddle and compact `rxn_12_1` both connect W2/W3; the old lower `rxn_12_2` instead connects W2/complex B at L2 | Independent L2 integration of all three; coupled-cluster refinement of the two W2/W3 saddles, T/Q energies and mapped geometry/mode comparison |
| W1 identity and two vdW arrangements | Complex A and B differ by about 0.90 Å mapped RMSD; both are observed minima | SI W1 coupled-cluster optimization/frequencies and T/Q energies at both observed complexes; retain separate identities unless equivalence is demonstrated |
| CH2SH H loss | Scale-1.6 seed gives a first-order saddle connecting CH2SH to CH2S···H; scale-1.9 duplicates it | Independent native L2 IRCs and matched T/Q saddle/endpoint energies, one representative |
| Direct CH3SS H loss | Previous guesses give internal rearrangements or a product-side minimum, not the proposed direct saddle | Two different product-side H-addition seeds, C–H 1.5/1.8 Å, followed by independent IRCs only if a first-order saddle is obtained |

`evidence.json` contains the frozen observations, mapped geometry/mode audit,
SI coordinate provenance, exclusion ledger and exact seeds. Its source text is
`r8-direct-l2-results.DsXyP5iu.txt`, SHA-256
`2af5e4050709a976e21d6b483dcd621b65de482fcd67629f00e16301c02548b1`.
The SI geometry source has SHA-256
`2c21c3f6905d07ad652a6b56be8f91284ae44154fb890b60d237326d7cbf0863`.

## Cost and reuse

| Work | Scope |
|---|---|
| Remaining retained connections | 44 saved L2 saddle observations, reused with exact native input/output hashes; two L2 IRCs and endpoint frequencies per observation |
| Independent native Gaussian paths | 15 selected observations, one Hessian followed by both directions at 0.10 and 0.05 sqrt(amu) bohr steps |
| New H-loss guesses | Two product-side saddle searches, with the same path checks when a saddle validates |
| L1 exclusion audit | 17 saddle refinements with two IRCs each, plus seven distinct fragment calculations, covering all 24 recorded decisions |
| Coupled-cluster geometry work | Five saddle and four minimum refinements/frequencies; both QSD directions plus endpoint optimization/frequencies for TS2, TS3 and TS5 only |
| Electronic benchmarks | T/Q single points at the nine coupled-cluster-refined structures and 19 selected existing L2 structures |

The native dispatcher contains **69 groups / 220 tasks**. The 44 reused-saddle
checks require up to 176 additional IRC/endpoint-frequency jobs. These are task
counts, not equal-cost calculations or 396 new transition-state searches.
The expensive numerical coupled-cluster paths are restricted to three specific
connection disagreements. Failed prerequisites block only their own descendants;
relaunching the controller does not reset an exhausted attempt.

Gaussian uses B2PLYP-D3(BJ)/cc-pVTZ, an ultrafine grid, tight SCF, a native
numerical Hessian and mass-weighted EulerPC IRCs with `GradientOnly` and `RCFC`.
The two step sizes have matched maximum path lengths. The final optimization
and frequency calculation is a separate Link1 step: its normal termination
does not establish that the preceding IRC converged. Review both sections.
See the [Gaussian IRC documentation](https://www.conflex.co.jp/gaussian_support/irc.php).

The independent coupled-cluster convention is explicitly
ROHF-UCCSD(T)-F12b with `scale_trip=1`, i.e. **(T*)**, frozen core,
cc-pVDZ-F12 geometry/frequencies and cc-pVTZ/QZ-F12 energies. It is not an exact
replication of an unspecified F12a/b/triples convention in the article/SI.
The native output retains the triples contributions and diagnostics. QSD
optimization/path history and gradients are printed; see the
[Molpro OPTG manual](https://www.molpro.net/manual/doku.php?id=geometry_optimization_optg)
and [F12 conventions](https://www.molpro.net/manual/doku.php?id=explicitly_correlated_methods).

## Scientific decision requirements

- Match charge, spin, theory, atom mapping and energy reference before comparing
  barriers. Do not compare an energy relative to W1 with one relative to W4.
- Require a converged first-order saddle and inspect its reactive displacement.
  Native IRC convergence, the actual path and optimized endpoint connectivity
  must agree; a final minimum alone cannot repair an incomplete or diverted IRC.
- For disputed L2 paths, compare the two integration step sizes. For the three
  coupled-cluster paths, inspect the native descent and stopping reason before
  assigning endpoints. A path that reaches its point limit is not approved by
  a subsequent unconstrained optimization.
- Retain distinct parallel saddles. Remove duplicates only using the existing
  atom-mapped geometry and reactive-mode evidence, including an explicit mirror
  mapping where applicable. Similar endpoints or a lower energy alone are
  insufficient.
- Examine T1, convergence, reference/state continuity and basis dependence before
  using coupled-cluster benchmarks. The paper itself reports difficulty near
  the vdW region; a single-reference result is not an authority there merely
  because it is higher level. Existing orbital reviews remain evidence for their
  sampled chemistry; they do not automatically certify a newly diverted path.
- Re-evaluate recorded 150-kcal/mol exclusions using L2 E+ZPE and the existing
  matched L2 initial-well reference. Atomic ZPE is zero; fragment multiplicities
  are preserved. Do not combine L1 reference energies with new L2 products.
- An unsuccessful new saddle search is a bounded search outcome, not proof that
  a physical channel cannot exist. This panel tests the known disputes and
  retained model; it cannot prove global PES completeness or exact experimental
  rates. PST, finite basis, electronic-method and statistical-mechanics
  approximations remain explicit scouting assumptions.

No result is automatically written into MESS. Connections and exclusions must
first be reconciled against the above evidence; MESS remains paused as requested.

## Run on Blodgett

From the updated `composite` checkout:

```bash
bash examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/run_blodgett.sh
```

The launcher runs detached. It preserves existing R8 allocations/controllers
and performs no controller refresh. It requires Gaussian and Molpro through
the existing site setup. Re-running the launcher adopts this campaign's saved
work and its controller lock prevents duplicate ownership.

Read-only status:

```bash
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_CHECK="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114"
"$PWD/.venv/bin/python" -B -m kinbot.saddle_recheck status "$KB_CHECK/remaining_l2"
"$PWD/.venv/bin/python" -B examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/campaign.py status "$KB_CHECK/native"
```

When finished, collect full native outputs, inputs, paths, geometry files,
results and consistent database snapshots. This writes only a new export
directory; it does not walk checkpoint stores or rerun calculations:

```bash
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
"$PWD/.venv/bin/python" -B -u \
  examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/collect.py collect \
  "$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114"
```

It prints the verified `validation-evidence.tar.gz` path and a compact diagnostic
`review.txt`. The archive is needed for complete path/mode review; text summaries
alone do not contain all native vectors and path geometries.

## Local validation

Tests exercise the real dispatcher with a simulated scheduler/execution layer:
live-job preservation, delayed execution publication, interrupted preparation,
dependency release, terminal failure, repeated invocation, changed definitions,
native Hessian validation and bounded evidence collection. The 44 retained
native source pairs were also staged successfully against the local R8 export.
Gaussian/Molpro calculations have not been run locally for this new panel;
neither unit tests nor staging constitute successful native chemistry.
