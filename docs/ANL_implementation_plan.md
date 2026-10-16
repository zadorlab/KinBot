# KinBot profiled theory, CBH-ANL, ROTD_py, and MESS plan

This document records the current design and remaining work. Historical
external-site transcripts, one-off recovery commands, and completed HPC smoke
harnesses were removed; they remain available in Git history. Reusable
validation belongs in ordinary inputs and general orchestration code.

## Scope

The extension is opt in. Inputs that do not request profiled theory,
CBH-ANL, or ROTD_py keep the established KinBot behavior and input format.
The added path provides:

1. FairChem UMA or another configured calculator for L1 search work.
2. B3LYP/cc-pVTZ or B2PLYP-D3(BJ)/cc-pVTZ for L2 geometry, frequencies,
   VPT2, and hindered rotors.
3. Declarative ANL component calculations through Gaussian, Molpro, CFOUR,
   and MRCC.
4. Connectivity-derived CBH reactions and pinned ATcT CBH-0 references.
5. Recursive CBH laddering for formation enthalpies at 0 K and finite
   temperature.
6. Strict transfer of composite energies, formation enthalpies, frequencies,
   and ROTD_py number-of-states data into MESS.
7. Restartable Slurm execution with native-output hashes and parser
   provenance.

## Compatibility rules

- Existing KinBot parameters retain their meaning.
- Profile selection occurs per job; the package-wide legacy backend remains
  the fallback only for jobs without a profile.
- Geometry, energy, frequencies, Hessian, and ZPE from an accepted result stay
  coherent. Data from unrelated calculations are never silently combined.
- If the optional L3 optimization is unavailable, only its L3 harmonic
  calculation is unavailable. Independent L3 and ANL single points run on the
  accepted L2 geometry and record that fallback in their recipe provenance.
- If VPT2 is unavailable, thermochemistry uses the accepted harmonic ZPE and
  labels the selected ladder entry as harmonic rather than anharmonic.
- Site paths and executable locations are discovered from the environment,
  modules, and user configuration. Repository code contains no cluster path.
- Licensed calculations are not repeated when a hash-verified native result
  can be reparsed.
- A named ANL recipe is emitted only when every required component is present
  at the declared method and state.

## Theory hierarchy

The built-in higher profile is:

```text
L1  FairChem UMA
L2  B2PLYP-D3(BJ)/cc-pVTZ
L3  requested ANL recipe
```

The automatic default for every ANL head is B2PLYP-D3(BJ)/cc-pVTZ at L2.
The explicit lower profile substitutes B3LYP/cc-pVTZ for reproducing the
published anharmonic correction. Gaussian VPT2 is a frequency-only
calculation on the accepted L2 geometry. Hindered-rotor scans also use the
corresponding L2 surface.

All conventional and F12 coupled-cluster components use an RHF determinant
for closed shells or an ROHF determinant for open shells, followed by the
unrestricted coupled-cluster method requested by the recipe. Molpro inputs
therefore use `rhf` followed by `uccsd(t)` or its F12 form. Higher-order MRCC
and CFOUR inputs must state and verify the supported determinant and
unrestricted correlation variant in their native output.

Electron-count reductions are exact:

- one electron: one HF calculation at each required basis;
- two electrons: one CCSD or CCSD-F12 calculation at each required basis;
- triples and higher increments that are algebraically zero are not staged.

Derivative jobs remain separate when geometry gradients or Hessians are
required, even when their electronic rank reduces.

## ANL energy model

A composite record stores each native component once and gives explicit
signed expressions for the electronic energy and zero-point energy.  Define

```text
CBS_p(E_(n-1), E_n) = E_n + alpha_n (E_n - E_(n-1))
alpha_n              = (n - 1)^p / (n^p - (n - 1)^p)
p                    = 3.7

delta_anh(M/B) = ZPE_anh(M/B) - ZPE_har(M/B)
delta_CV       = CBS_3.7[E_CCSD(T,all-electron)(cTZ,cQZ)]
               - CBS_3.7[E_CCSD(T,frozen-core)(cTZ,cQZ)]
delta_REL      = E_CCSD(T,all-electron,DKH2)(acTZ-DK)
               - E_CCSD(T,all-electron,nonrelativistic)(acTZ-DK)
```

Every CBS extrapolation in the four built-in profiles uses the same inverse
angular-momentum power of 3.7.  Ram et al. reported different fitted alpha
values for their calculations; those fitted values are not used by these
profiles.  Each result records the exponent, cardinal numbers, basis labels,
geometry hashes, and native output hashes.

The implemented equations are:

```text
ANL0 = CBS_3.7[CCSD(T)/aug-QZ,aug-5Z] // CCSD(T)/TZ
     + ZPE_har[CCSD(T)/TZ]
     + delta_anh[L2/TZ]
     + {CCSDT(Q)-CCSD(T)}/DZ
     + delta_CV
     + delta_REL
     + DBOC[HF/TZ]
     + SO

ANL0-F12 = 0.5 CBS_3.7[CCSD(T)-F12b/TZ-F12,QZ-F12]
           + 0.5 CBS_3.7[CCSD(T)/aug-QZ,aug-5Z]
           all on the CCSD(T)/TZ geometry
           + ZPE_har[CCSD(T)/TZ]
           + delta_anh[L2/TZ]
           + {CCSDT(Q)-CCSD(T)}/DZ
           + delta_CV
           + delta_REL
           + DBOC[HF/TZ]
           + SO

ANL1 = CBS_3.7[CCSD(T)/aug-5Z,aug-6Z] // CCSD(T)/QZ
     + CBS_3.7[ZPE_har(CCSD(T)/TZ), ZPE_har(CCSD(T)/QZ)]
     + delta_anh[L2/TZ]
     + {CCSDT(Q)-CCSD(T)}/TZ
     + {CCSDTQ(P)-CCSDT(Q)}/DZ
     + delta_CV
     + delta_REL
     + DBOC[HF/TZ]
     + SO

ANL1-F12 = 0.5 CBS_3.7[CCSD(T)-F12b/QZ-F12,5Z-F12]
           + 0.5 CBS_3.7[CCSD(T)/aug-5Z,aug-6Z]
           all on the CCSD(T)/QZ geometry
           + CBS_3.7[ZPE_har(CCSD(T)/TZ), ZPE_har(CCSD(T)/QZ)]
           + delta_anh[L2/TZ]
           + {CCSDT(Q)-CCSD(T)}/TZ
           + {CCSDTQ(P)-CCSDT(Q)}/DZ
           + delta_CV
           + delta_REL
           + DBOC[HF/TZ]
           + SO
```

`aug-nZ` denotes aug-cc-pVnZ on every atom, including hydrogen; `cVnZ` is
cc-pCVnZ; and `acVTZ-DK` is aug-cc-pCVTZ-DK.  The ANL0, ANL0-F12, and ANL1
names retain their correction hierarchy, but these production profiles use
full hydrogen augmentation.  ANL0-F12 averages its TZ/QZ-F12 branch with the
fully augmented conventional QZ/5Z limit; ANL1-F12 averages QZ/5Z-F12 with
the fully augmented conventional 5Z/6Z limit.  Setting
`anl_f12_augmented` to `1` replaces the F12 branch with the corresponding
fully augmented `aug-cc-pVnZ-F12` pair and records that choice in the recipe
label and workflow intent.  Every CBS term uses the 3.7 power above.

The published anharmonic term is B3LYP/cc-pVTZ.  Automatic production for
all four ANL heads defaults to B2PLYP-D3(BJ)/cc-pVTZ and records that choice
in the recipe label.  A harmonic-only fallback is also labeled explicitly.
`SCALE_TRIP=1` is pinned for Molpro F12b calculations and is likewise
included in the profiled recipe label.

Closed shells use an RHF determinant and open shells use a semicanonical
ROHF determinant.  All post-HF coupled-cluster correlation is unrestricted:
RHF-UCC for closed shells and ROHF-UCC for open shells.  Direct MRCC
CCSDT(Q) and CCSDTQ(P) parsing selects the native `/B` perturbative variant.
The ordinary coupled-cluster terms are frozen core; the all-electron side of
`delta_CV` and both sides of `delta_REL` are explicit exceptions.  Exact
one-electron and two-electron rank reductions are applied before submission.

ANL1-only components remain excluded from an ANL0 result. A fallback result
keeps the actual recipe label; it is never renamed to the requested ladder
head.

## CBH reaction generation

The generator works directly from state-explicit molecular connectivity:

1. Parse and canonicalize the supplied SMILES.
2. Verify charge, electron parity, multiplicity, and the KinBot atom
   inventory.
3. Canonicalize the graph and convert aromatic systems to a deterministic
   alternating single/double Kekule representation.
4. Retain atom-local formal charge, radical electrons, isotope, and stereo
   metadata while building the heavy-atom graph.
5. Generate hydrogen-capped connected fragments for the requested CBH rung.
6. Expand ring overlaps by inclusion-exclusion, split disconnected overlap
   terms into connected molecular fragments, and collect identical fragments.
7. Add the hydrogen-balancing species required by the radical rule.
8. Verify elemental and formal-charge balance before accepting the equation.

The implementation is independent and does not copy a public CBH library.
The same graph rules cover acyclic, cyclic, fused aromatic, heteroaromatic,
ionic, zwitterionic, radical, multiple-bond, and hypervalent covalent SMILES.
Equivalent SMILES traversals produce the same canonical reaction.  Malformed,
disconnected, charge-inconsistent, or electron-parity-inconsistent identities
are rejected before any calculation is planned.

Signed coefficients are negative for reactants and positive for products.
The target coefficient is always `-1`.

## Formation enthalpy at 0 K

For a reaction with coefficients `nu_i`, the calculated 0 K reaction energy is

```text
Delta_r E(0) = sum_i nu_i E_i(0).
```

The target formation enthalpy is

```text
Hf_target(0) = sum_(j != target) nu_j Hf_j(0) - Delta_r E(0).
```

CBH-0 is the only rung that reads direct pinned ATcT formation enthalpies.
Every CBH-1, CBH-2, or CBH-3 immediate reference carries a computed
`FormationEnthalpy` from a lower-rung reaction. The solver descends until all
leaves are CBH-0 ATcT references. Each reaction node may select its own
highest complete method tier, which implements ANL laddering without forcing
a large target and all smaller references onto one method.

The per-state method orders are:

```text
ANL1-F12 -> ANL1 -> configured L3 -> L3 single point on L2 -> L2
ANL1     -> configured L3 -> L3 single point on L2 -> L2
ANL0-F12 -> ANL0 -> configured L3 -> L3 single point on L2 -> L2
ANL0     -> configured L3 -> L3 single point on L2 -> L2
```

ANL1 does not silently fall through to ANL0 within one state because their
optimized geometries, harmonic terms, and higher-order increments differ.
The CBH size policy may assign different requested recipe heads to different
states, and different CBH nodes may therefore use different exact levels.

A node uses one exact common method for every energy in that reaction. Missing
or inconsistent data move that node to the next tier or stop the solve.
Different nodes do not need the same method or component inventory.  In
particular, a target reaction may use ANL1-F12 while one of its reference
formation enthalpies is supplied by an ANL0 or L3 reaction below it.  Only
the species participating in the same individual reaction node must share
one exact level label.

Native runtime, memory, scheduler, or optimizer failures in an optional ANL
component mark that component unavailable.  The dispatcher does not submit
the explicitly related tasks that remain unsent, but it continues independent
components.  Parser, provenance, method, and quality mismatches remain hard
failures because treating those as missing calculations would conceal a code
or chemistry error.  A named ANL recipe is still assembled only when every
term in that exact recipe is present; otherwise the CBH solver selects the
next complete composite profile (ANL0, configured L3, or L2) supplied for that
species.  It never relabels a partial result as the failed ANL tier.

## Formation enthalpy at finite temperature

MESSPF supplies each species' ideal-gas heat content:

```text
H_i(T) - H_i(0) = R T^2 d ln(Q_i)/dT + R T.
```

The same balanced CBH reaction is then evaluated at temperature `T`:

```text
Delta_r H(T) = Delta_r E(0)
             + sum_i nu_i [H_i(T) - H_i(0)]

Hf_target(T) = sum_(j != target) nu_j Hf_j(T) - Delta_r H(T).
```

The reference rule is identical to the 0 K rule:

- CBH-0 uses pinned ATcT values at 0 and 298.15 K.
- Higher rungs recursively use lower-rung computed formation enthalpies at
  both temperatures.
- Direct ATcT use above CBH-0 is invalid.
- Elemental reference-state heat increments are not part of this path.

MESSPF inputs are therefore generated for every species in the recursive CBH
reaction tree. Species in the live KinBot network retain their accepted L2
geometry, accepted VPT2 modes when available, and accepted L2 hindered rotors.
External CBH references use hash-traceable thermochemistry models from
completed KinBot/composite runs. A failed VPT2 calculation falls back to
harmonic modes; a failed rotor scan remains explicitly identified as a
harmonic-oscillator fallback in the MESSPF provenance.

NASA-7 fits use the CBH-derived Hf(298.15 K) as their enthalpy anchor and the
native MESSPF entropy as their entropy anchor. The audit reparses every native
table and independently reconstructs the recursive CBH temperature identity.

## KinBot to MESS handoff

A strict handoff maps exact KinBot routing identifiers to:

- one accepted composite 0 K energy record;
- one recursive CBH formation record;
- one accepted frequency source;
- one accepted L2 geometry;
- the existing L2 hindered-rotor representation.

Stable-well and product relative energies use differences in ANL 0 K
formation enthalpies. Stationary barriers use consistent composite energy
differences. Barrierless channels consume the selected, hash-verified ROTD_py
number-of-states surface through MESS `Core Rotd`. The traditional KinBot
energy path remains active when CBH-ANL is not requested.

KinBot preserves the hash-verified native ROTD_py number-of-states table and
stages a separate MESS copy with the conventional zero-flux origin.  When the
native zero-energy row has positive flux, the staged table follows the
published MESS convention: `0 0`, then that positive-flux row at 0.01 cm-1.
Both hashes and the transformation are recorded in `rotdpy_sources.json`.

Thermochemistry fitting uses the same CBH reaction network at 298.15 K. It
does not switch to an atomization or elemental-reference equation.

## ROTD_py production path

The completed production methyl recombination validation used:

```text
sampling:    CASPT2(2e,2o)/cc-pVDZ
correction:  MRCI+Q(2e,2o)/cc-pVTZ
surfaces:    24, from 2.3 to 4.6 angstrom
```

A production result must contain the multipoint correction potential, pivot
and dividing-surface definitions, converged Monte Carlo sampling, nonempty
numeric surface data, number-of-states output, and hashes tying outputs to
the generated input. Reduced interface smoke data cannot satisfy the
production gate.

The sampling basis may be one basis or an adjacent two-basis list.  A list is
extrapolated point by point with the same inverse-power expression used by the
ANL recipes; `vrc_tst_sample_extrapolation_power` defaults to 3.7.  Thus a
high-accuracy input can request, for example,

```json
"vrc_tst_sample_method": "caspt2-f12(9,7)",
"vrc_tst_sample_basis": ["avtz-f12", "avqz-f12"]
```

Spin splitting is state neutral.  `target_multiplicity`,
`reference_multiplicity`, and the two `fragment_multiplicities` are checked
against angular-momentum coupling before any task is submitted.  This covers
singlet/triplet surfaces from two doublets and doublet/quartet surfaces from a
doublet plus a triplet. Explicit CASPT2 and MRCI active electron/orbital
counts are preserved, while automatic policies replace template counts; both
paths are validated for electron parity at every requested spin.
Automatic reference selection uses the maximum allowed spin, including when
fragment spins are declared in the input. An already-maximum-spin target
requires direct correction or an explicit expert reference; it never silently
uses a lower-spin CC anchor. Explicit settings cannot change the accepted
parent or product spins.
Unqualified F12 multireference methods select variant B; MRCI+Q selects the
relaxed-reference Davidson energy (`ENERGD`).

`vrc_tst_active_space` has four policies. `explicit` preserves the CAS in
each method string. `full_valence` derives a whole-molecule main-group AVAS
space from atom-labelled valence shells in C1 and replaces the template CAS
in every scan, sampling, and spin-gap method. `radical_pair` is an opt-in
CAS(2,2) shortcut accepted only for two doublet fragments on coupled singlet
and triplet surfaces, such as CH3 + CH3. It is rejected for doublet + triplet
and other fragment-spin combinations. The `anl1-f12` VRC accuracy preset
requires an explicit active-space policy or orbital/surface model. It specifies
energy levels without choosing `full_valence` or `radical_pair` implicitly.

The fourth policy, `selected`, supplies counted active electrons/orbitals
and requires an explicit orbital model. The sulfur examples use provisional
selected spaces with basis-specific external model files; their counts are
not a certification of orbital identity. `electronic_state_models` supplies
ground spins, optional resolved MESS levels and spin-free term degeneracies.
A term-only record controls spatial root counts without inventing measured
electronic level spacings. In particular, retaining CH3S X(2E) and S(3P)
gives six spatial product components per coupled spin manifold.

The [14 October audit](ANL_workflow_audit.md) describes the general adaptive
reconnaissance, reviewed-record gate, software corrections, remaining
chemistry acceptance checks, and the current default-F12 basis instruction.

An expert orbital model can use a prior Molpro wavefunction instead of AVAS:

```json
"vrc_tst_state_model": {
  "orbital_selection": {
    "mode": "molpro_record",
    "path": "/absolute/path/selected.wfu",
    "record": "2140.2",
    "sha256": "<64 lowercase hexadecimal characters>"
  },
  "state_average": [
    {"multiplicity": 1, "roots": 1},
    {"multiplicity": 3, "roots": 1}
  ]
}
```

KinBot validates the record and hash before submission, copies it atomically
into each calculation directory, rechecks the hash, and emits Molpro `FILE`
and `START` directives. This makes a manually selected orbital ordering
reproducible without embedding molecule-specific orbital numbers in KinBot.

Automatic ANL calculations use `anl_scf_convergence=1e-10` and
`anl_cc_convergence=1e-8` by default.  The values are translated into native
Gaussian, Molpro, CFOUR, and MRCC controls.  VRC calculations independently
use `vrc_tst_scf_convergence=1e-10` and
`vrc_tst_correlation_convergence=1e-8`, allowing expert inputs to change VRC
thresholds without weakening the thermochemistry graph.

ROTD_py owns and records every submitted sample job. Resume logic counts
completed samples, avoids duplicate submissions, stops orphan controllers,
and never resubmits after every surface has converged.

## Scheduler and runtime behavior

The generic dispatcher:

- requests user-specified nodes and exclusive allocation;
- discovers node CPU and memory capacity from Slurm;
- caps ranks by the method's useful scaling range and minimum memory per rank;
- fails before submission if the minimum cannot be met;
- isolates Molpro scratch and repository files per invocation;
- detects Gaussian scratch and library requirements on compute nodes;
- discovers Gaussian, Molpro, CFOUR, MRCC, modules, basis files, and runtime
  libraries without cluster-specific constants or prior module loads; it
  searches explicit configuration, `PATH`, conventional versioned install
  prefixes, and available Environment Modules/Lmod entries, then pins the
  resolved setup in each worker script;
- reconciles jobs that have left `squeue` without treating an invalid job ID
  as a dispatcher failure;
- stages all dependency-ready post-geometry tasks concurrently under the
  global job limit;
- preserves and reparses native outputs before considering a rerun.

MRCC restart support retains a nonempty amplitude checkpoint and adds the
native restart option. Failed attempts are archived with their resource and
input provenance.

## Validation contract

Unit tests establish input generation, parsing, algebra, restart state, and
handoff behavior. Chemistry validation uses a small number of published
comparisons rather than rerunning every expensive component:

- source-geometry CH4 and CH3 higher-order increments;
- ethane CCSD(T) and higher-order differences;
- ethane and methyl CBH-ANL Hf(0 K);
- native VPT2 and harmonic ZPE parsing;
- one complete L2 hindered rotor;
- one production ROTD_py calculation;
- MESSPF partition functions and recursive CBH Hf(298.15 K);
- CH3 + CH3 high-pressure kinetics from MESS.

The final report contains native hashes, selected recipe labels, component
sums, CBH reaction trees, MESSPF temperature traces, NASA continuity checks,
ROTD_py source hashes, MESS output hashes, and independent literature
comparisons.

## Supported production boundary

The ordinary KinBot input drives ANL0, ANL0-F12, ANL1, or ANL1-F12. It
freezes the accepted reaction network, generates every CBH reference graph,
runs missing reference conformer/L2/frequency/hindered-rotor work, schedules
the size-selected composite DAGs, selects the best complete common tier at
each reaction node, solves 0 and 298.15 K formation enthalpies, applies the
handoff, and runs MESSPF/ROTD_py/MESS when requested.

The exact per-state tier orders are listed in the formation-enthalpy section
above.  Size-infeasible states fall to L3 or L2 under the current three-band
size policy.  A future independent ANL1-to-ANL0 size boundary would require a
new user-visible policy key; it must not be inferred by relabeling ANL1 work.

The graph generator covers CBH-0 through CBH-3 for covalent SMILES carrying
formal charges and atom-local radical electrons.  External CBH references are
ordinary restartable KinBot calculations and retain real L2 rotor models.
The general reaction selector uses family names plus atom motifs, and the
product-search stop is a normal input option rather than a ChemID-specific
test hook.

The internal NASA-7 output is identified as KinBot's fitted representation;
it is not labeled as PAC99 output.  Native ATcT, electronic-structure,
partition-function, ROTD_py, and MESS hashes remain in the final report so an
independent uncertainty or alternate fitting workflow can consume the same
verified data without changing the computed chemistry.

## Repository organization

```text
kinbot/anl/atct.py             pinned ATcT API/cache reader
kinbot/anl/cbh.py              CBH graph reactions and recursive Hf solver
kinbot/anl/tasks.py            declarative QC task construction
kinbot/anl/dispatch.py         restartable scheduler execution
kinbot/anl/results.py          native result parsers
kinbot/anl/workflow.py         verified component extraction and CBS algebra
kinbot/anl/recipes.py          named composite equations
kinbot/anl/validation.py       composite assembly and scientific audits
kinbot/anl/handoff.py          strict network attachment
kinbot/anl/thermochemistry.py  MESSPF parsing, CBH Hf(T), NASA fitting
kinbot/anl/end_to_end.py       final reusable verification orchestration
kinbot/rotdpy.py               ROTD_py execution and provenance
kinbot/mess.py                 MESS/MESSPF serialization and execution
```

The next compact external-site acceptance network and its scientific checks
are defined in `docs/ANL_small_system_validation.md`; implementation history
belongs in Git, not in this plan.
