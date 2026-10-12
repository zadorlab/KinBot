# KinBot profiled theory, CBH-ANL, ROTD_py, and MESS plan

This document records the current design and remaining work. Historical
external-site transcripts and one-off recovery commands were removed; they
remain available in Git history. The executable validation procedure is in
[`ANL_full_HPC_validation.md`](ANL_full_HPC_validation.md).

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

The lower profile substitutes B3LYP/cc-pVTZ at L2. Gaussian VPT2 is a
frequency-only calculation on the accepted L2 geometry. Hindered-rotor scans
also use the corresponding L2 surface.

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

A composite record stores each component once and gives explicit signed
expressions for electronic energy and ZPE. The common ANL0-F12 profile uses

```text
E0 = E_ref,CBS(F12)
   + Delta_higher_order
   + Delta_core_valence,CBS
   + Delta_scalar_relativistic
   + DBOC
   + spin_orbit
   + ZPE_harmonic
   + Delta_ZPE,anharmonic
```

Every difference is evaluated on the same accepted highest-level geometry
available for that recipe. The selected program or reference determinant may
differ between components when the underlying method is the same and the
native program supports the requested calculation.

Two-point CBS terms use the recipe's pinned coefficients. A generic pair is
stored as

```text
E_CBS = c_low E_low + c_high E_high
```

with the coefficients, basis labels, source hashes, and extrapolated quantity
in the result. Harmonic ZPE and the anharmonic correction may each be
extrapolated from TZ/QZ when both calculations are available. `SCALE_TRIP=1`
is retained for Molpro F12 calculations where the method requires it.

ANL1-only components remain excluded from an ANL0 result. A fallback result
keeps the actual recipe label; it is never renamed to the requested ladder
head.

## CBH reaction generation

The generator works directly from explicit-hydrogen molecular connectivity:

1. Parse and canonicalize the supplied SMILES.
2. Verify charge, multiplicity, supported valence, and the KinBot atom
   inventory.
3. Build the heavy-atom graph with integral bond orders.
4. Generate capped fragments for the requested CBH rung.
5. Canonicalize identical fragments and collect integer coefficients.
6. Add the hydrogen-balancing species required by the radical rule.
7. Reject a self-referential or elementally unbalanced equation.

The implementation is independent and does not copy a public CBH library.
Unsupported aromatic, ionic, or unusual-valence cases fail closed until a
state-aware fragmentation rule is implemented and tested.

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

The default method order is:

```text
ANL1-F12 -> ANL1 -> ANL0-F12 -> ANL0 -> configured L3 -> L2
```

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
geometry, VPT2 modes, and L2 hindered rotors. External CBH references use
hash-traceable thermochemistry models from completed KinBot/composite runs.
A rotatable external reference requires a real L2 rotor model; it cannot be
silently represented as RRHO.

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

Thermochemistry fitting uses the same CBH reaction network at 298.15 K. It
does not switch to an atomization or elemental-reference equation.

## ROTD_py production path

The production methyl recombination validation uses:

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
  libraries without cluster-specific constants;
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

## Current acceptance boundary

The completed implementation supports ANL0 and ANL0-F12 component assembly,
recursive CBH reaction solves, ATcT CBH-0 leaves, reviewed VPT2 corrections,
L2 hindered rotors, ROTD_py production surfaces, and MESS/MESSPF handoff.

The following remain explicit future work:

1. Finish and validate every pinned ANL1 and any separately defined ANL1-F12
   component equation.
2. Add state-aware aromatic, ionic, and unusual-valence CBH fragmentation.
3. Automate generation of complete external CBH reference KinBot workflows,
   including their conformers and hindered rotors.
4. Propagate ATcT, electronic-structure, CBH ladder, and partition-function
   uncertainties into NASA and kinetics products.
5. Validate a conventional stationary transition state, IRC, and restart path
   in addition to the barrierless methyl recombination case.
6. Compare the internal NASA fitting result with the chosen production PAC99
   workflow before claiming interchangeable output.

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

The external-site commands and current acceptance checklist are maintained in
`docs/ANL_full_HPC_validation.md`; implementation history belongs in Git, not
in this plan.
