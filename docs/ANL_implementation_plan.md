# KinBot Profiled ASE/Sella + ANL Composite Methods
## Final revised implementation plan

**Repository reviewed:** `zadorlab/KinBot`, current `master` as of 2026-09-16
**Observed upstream state during review:** KinBot 2.4.1; current master head reported by GitHub activity as `924faba`; no open pull requests at the time of review.
**Primary objective:** implement a compact, opt-in, backward-compatible hierarchy

```text
UMA / FairChem L1
        ↓
B3LYP/cc-pVTZ L2 (lower ANL tier)
or B2PLYP-D3(BJ)/cc-pVTZ L2 (higher ANL tier)
        ↓
ANL0 / ANL0-F12 / ANL1 / future pinned ANL-family recipe
```

while retaining full user control over L1, L2, and every constituent L3/composite calculation.
For each tier, VPT2 and hindered-rotor scans use that tier's L2 surface.
The VPT2 frequency job consumes the accepted L2 geometry directly; Gaussian
does not optimize that geometry unless `Opt` is explicitly requested.

The implementation should:

1. use ASE calculators as the program interface;
2. use Sella for geometry optimization and constrained hindered-rotor optimization;
3. use the most appropriate Hessian/frequency path for each backend;
4. route single points to Molpro, CFOUR, MRCC, or Gaussian according to a declarative recipe;
5. preserve all existing KinBot behavior when the new feature flag is not active;
6. make only narrow, flag-guarded edits to existing hot files;
7. keep method definitions, software details, algebra, and result provenance separate;
8. avoid redundant high-level calculations;
9. never silently alter a named ANL method;
10. remain restartable, auditable, and practical on a cluster.

---

# 1. Core conclusions from the repository review

The current KinBot architecture is already close to what is needed, but three global assumptions must be isolated behind a new opt-in path.

## 1.1 KinBot already vendors custom ASE calculators

Current source tree includes:

```text
kinbot/ase_modules/calculators/
    gaussian.py
    qchem.py
    orca.py
```

The current Gaussian calculator is a KinBot-owned `FileIOCalculator` implementing:

```text
energy
forces
dipole
```

The current ORCA wrapper already follows ASE 3.26's newer:

```text
BaseProfile
CalculatorTemplate
GenericFileIOCalculator
```

pattern.

This establishes the correct precedent.

**Recommendation:** new Molpro, CFOUR, and MRCC interfaces should live in the same package and should preferentially use `GenericFileIOCalculator`, following the current ORCA architecture.

Do not create ANL-specific shell scripts that bypass ASE for ordinary energy/gradient calls unless a program capability genuinely requires a separate auxiliary calculation.

## 1.2 Current Sella templates already instantiate KinBot calculators directly

For example, the current Sella optimization path effectively does:

```python
from kinbot.ase_modules.calculators.<code> import <Code>
mol.calc = <Code>(**kwargs)
opt = Sella(mol, ...)
```

and current hindered-rotor Sella logic similarly places a calculator on an ASE `Atoms` object and applies a `sella.Constraints` dihedral constraint.

Therefore the desired architecture is evolutionary, not a rewrite:

```text
current custom ASE calculator pattern
        +
new calculator factory/capabilities
        +
per-calculation theory profile
```

## 1.3 The present blocker is global backend selection

`QuantumChemistry.__init__` currently stores one:

```python
self.qc = par["qc"]
self.qc_command = par["qc_command"]
```

and `get_qc_arguments()` begins by branching on `self.qc`.

For a Gaussian calculation, `high_level=1` changes:

```python
method
basis
```

but not the QC backend.

Likewise `qc_hir()` asks for:

```python
high_level=1
hir=1
```

but then decides which calculator/template to use from global `self.qc`.

Consequently this does **not** work correctly today:

```text
L1 = FairChem/UMA
L2 = Gaussian B2PLYP-D3(BJ)/cc-pVTZ
```

without a small per-job theory/profile router.

## 1.4 Current result readers also assume one backend

Examples include:

```text
check_qc(job)
read_qc_hess(job, natom)
hessian_is_massweighted()
```

where output type or Hessian interpretation is selected from global `self.qc`.

That assumption becomes invalid as soon as one run contains:

```text
FairChem L1
Gaussian L2
Molpro L3 tasks
CFOUR L3 tasks
MRCC L3 tasks
```

The new path therefore needs persistent job metadata and a single `backend_for_job(job)` resolver.

## 1.5 KinBot 2.4.1 now preserves a coherent accepted result

`calculation.py` deliberately loads:

```text
geometry
energy
ZPE
frequencies
Hessian
```

from one selected result row.

This behavior must be preserved.

The new profile system must never synthesize a species object from, for example:

```text
UMA geometry
Gaussian energy
UMA frequencies
Gaussian Hessian
```

When L2 is accepted, the L2 record must remain coherent.

---

# 2. Scientific defaults and user control

## 2.1 Default hierarchy when the new preset is requested

The proposed built-in preset is:

```text
L1:
    FairChem / UMA

L2:
    Gaussian
    B2PLYP-D3(BJ)/cc-pVTZ

L3/composite:
    whichever ANL method the user requests
```

Example:

```json
{
  "theory_preset": "uma-b2plyp-anl",
  "composite_method": "ANL0-F12",
  "fc_model_path": "/path/to/uma.pt"
}
```

If the user requests the preset but does not request a composite method:

```text
UMA → B2PLYP-D3(BJ)/cc-pVTZ
```

is still valid.

If the user requests an ANL method with no explicit L1/L2 profiles, the preset can provide the recommended L1/L2 defaults automatically.

## 2.2 This must remain opt-in

Do not change the package-wide legacy default to FairChem.

Reasons:

- FairChem is an optional dependency.
- UMA requires a model checkpoint.
- many existing KinBot installations do not have FairChem;
- existing inputs must remain byte-for-byte semantically equivalent when no new flags are used.

Recommended activation logic:

```text
profiled_theory = 0 by default

automatically activate profiled_theory if any of:
    theory_preset is set
    l1_profile is set
    l2_profile is set
    composite_method is set
```

The user may also explicitly set:

```json
"profiled_theory": 1
```

## 2.3 User control over L1

Support either compact aliases:

```json
"l1": "uma"
```

or a full profile:

```json
"l1_profile": {
  "calculator": "fairchem",
  "model_path": "/path/to/model.pt",
  "task_name": "omol",
  "device": "cuda:0"
}
```

A non-UMA example:

```json
"l1_profile": {
  "calculator": "gaussian",
  "method": "wB97XD",
  "basis": "6-31+G(d,p)"
}
```

## 2.4 User control over L2

Default preset:

```json
"l2_profile": {
  "calculator": "gaussian",
  "method": "B2PLYP",
  "basis": "cc-pVTZ",
  "calculator_kwargs": {
    "EmpiricalDispersion": "GD3BJ"
  },
  "optimizer": "sella",
  "frequency_mode": "native_hessian"
}
```

The user can replace any of those values.

## 2.5 User control over L3 subcalculations

Every ANL constituent must have a stable task ID.

Example:

```json
"l3_overrides": {
  "reference_f12_qz": {
    "basis": "cc-pVQZ-F12"
  },
  "dboc": {
    "calculator": "cfour",
    "memory_mb": 32000
  },
  "higher_order_q": {
    "nproc": 16
  }
}
```

Separate **theory-changing overrides** from **resource-only overrides** if possible:

```json
"l3_resource_overrides": {
  "ccsdtqp_dz": {
    "nproc": 32,
    "memory_mb": 120000
  }
}
```

Resource changes must not make a canonical method “custom.”

Theory changes do.

## 2.6 Parameter precedence

Use:

```text
legacy KinBot defaults
        <
built-in theory preset
        <
explicit L1/L2 user profiles
        <
canonical composite recipe
        <
explicit L3 theory overrides
        <
runtime resource overrides
```

If an L3 theory override changes a canonical recipe term, store:

```text
protocol_exact = false
base_recipe = "ANL0-F12"
```

and use a label such as:

```text
custom:ANL0-F12
```

Do not continue to call it exact `ANL0-F12`.

---

# 3. Literature-grounded ANL design

## 3.1 Original ANL definitions

The original 2017 ANL paper defines:

```text
ANL0
ANL0-F12
ANL1
```

as complete 0 K composite methods.

The equations include:

- a reference CCSD(T)/CBS term;
- harmonic ZPE;
- anharmonic ZPE correction;
- higher-order excitation correction;
- core-valence correction;
- scalar relativistic correction;
- DBOC;
- spin-orbit correction.

The original ANL paper also states that the general RUCCSD(T) calculations use restricted-spin HF wavefunctions within the unrestricted coupled-cluster formulation in Molpro, while UUCCSDT(Q) calculations use UHF as required by the MRCC implementation. Vibrational frequencies were obtained by Molpro numerical differentiation in that implementation; DBOC by CFOUR; DFT calculations by Gaussian.

That is a critical reference-policy precedent.

## 3.2 Original ANL0 structure

Represent the exact published formula as data:

```text
E_ANL0 =
  E_CCSD(T)/CBS(a'QZ,a'5Z) // CCSD(T)/TZ
+ ZPE_harm[CCSD(T)/TZ]
+ (ZPE_anh[B3LYP/TZ] - ZPE_harm[B3LYP/TZ])
+ (E_CCSDT(Q)/DZ - E_CCSD(T)/DZ)
+ (E_CCSD(T,full)/CBS(cTZ,cQZ) - E_CCSD(T,frozen-core)/CBS(cTZ,cQZ))
+ REL_DKH
+ DBOC_HF/TZ
+ SO
```

The implementation should not encode this by procedural if/else chains.

It should encode absolute calculation nodes plus a linear expression.

## 3.3 Original ANL1 structure

```text
E_ANL1 =
  E_CCSD(T)/CBS(a'5Z,a'6Z) // CCSD(T)/QZ
+ ZPE_harm[CCSD(T)/CBS(TZ,QZ)]
+ (ZPE_anh[B3LYP/TZ] - ZPE_harm[B3LYP/TZ])
+ (E_CCSDT(Q)/TZ - E_CCSD(T)/TZ)
+ (E_CCSDTQ(P)/DZ - E_CCSDT(Q)/DZ)
+ (E_CCSD(T,full)/CBS(cTZ,cQZ) - E_CCSD(T,frozen-core)/CBS(cTZ,cQZ))
+ REL_DKH
+ DBOC_HF/TZ
+ SO
```

The published ANL1 expression above uses a B3LYP/TZ anharmonic correction.
The user-selected higher-tier profile instead uses B2PLYP-D3(BJ)/cc-pVTZ
for L2, hindered rotors, and VPT2. Its anharmonic correction must therefore
be recorded as a **profiled ANL1 variant**, with method provenance explicit;
it must not be reported as reproducing the published ANL1 formula. The lower
tier uses B3LYP/cc-pVTZ for L2, hindered rotors, and VPT2. Both tiers run
VPT2 frequency analysis at their accepted L2 geometries without a Gaussian
`Opt` route. A user-requested in-job optimization would be a separate,
explicit task policy.

## 3.4 Original ANL0-F12 structure

The published definition replaces the main conventional ANL0 reference term with:

```text
CCSD(T)-F12b/CBS(TF,QF)
```

while retaining the ANL0 geometry/frequency framework unless a distinct later variant is explicitly selected.

Canonicalize this internally to:

```text
CCSD(T)-F12b
```

never an ambiguous generic `-F12`.

## 3.5 Later ANL-family practice

Later ANL-family kinetic studies commonly use:

- Molpro for most CCSD(T), F12, core-valence, and relativistic pieces;
- CFOUR for DBOC and selected closed-shell/singlet CCSDT(Q);
- MRCC for open-shell CCSDT(Q) and CCSDTQ(P);
- Gaussian for DFT/VPT2 anharmonicity.

Later papers also use improved variants such as:

```text
ANL0F
ANL0-F12'
ANL1-QZF
```

with higher-level F12 geometries, altered ZPE choices, or other refinements.

These are **not synonyms** for the original ANL methods.

## 3.6 Ram et al. perfluoroalkyl work

“Association Kinetics for Perfluorinated n-Alkyl Radicals” is especially useful as a model for a flexible automated implementation because it uses ANL1/ANL0-style reference energies in a laddered thermochemistry workflow for fluorinated systems.

The important software lesson is not “add one hard-coded ANL1-F12 name.”

It is:

> represent high-level thermochemistry as a reusable graph of absolute calculations and algebraic corrections so that a recipe can scale to the highest feasible basis/correlation level on a species-by-species basis.

The paper’s published abstract describes ANL1/ANL0-style reference energies rather than defining a universally named canonical `ANL1-F12`.

Therefore the initial implementation should ship:

```text
ANL0
ANL0-F12
ANL1
```

and optionally a separately documented later recipe such as:

```text
ANL1-QZF
```

only after its exact equation has been pinned.

A Ram-style species-dependent ladder should be represented as:

```text
ANL-LADDER
```

or:

```text
ANL-EXTENDED
```

rather than falsely labeling it canonical `ANL1`.

---

# 4. F12 policy

This is an area where the implementation should be stricter than a handwritten workflow.

## 4.1 Always resolve generic F12 energy methods to F12b

For every **energy** node:

```text
CCSD(T)-F12
UCCSD(T)-F12
```

should normalize to:

```text
CCSD(T)-F12b
UCCSD(T)-F12b
```

unless the recipe explicitly requests F12a/F12c.

Store both:

```text
requested_method
resolved_method
```

in provenance.

## 4.2 Default `SCALE_TRIP=1` for energy-only F12b calculations

For energy nodes, use:

```text
SCALE_TRIP=1
```

by default.

This produces the familiar scaled triples correction based on the MP2-F12/MP2 correlation ratio.

Represent configuration as:

```python
scale_trip = "auto"
```

not as a global boolean.

Then:

```text
energy-only F12b task:
    auto -> SCALE_TRIP=1
```

## 4.3 Do not blindly use `SCALE_TRIP=1` for analytic F12 gradients

Molpro documents analytic gradients for DF-MP2-F12, DF-CCSD-F12, and DF-CCSD(T)-F12.

However:

> the scaled version of (T) is not available in the analytic-gradient implementation.

Therefore:

```text
F12 derivative task + scale_trip="auto":
    analytic gradient requested
    -> use scale_trip=0
```

unless the recipe explicitly says the derivative surface itself must use scaled triples.

## 4.4 If the user explicitly requests scaled F12 derivatives

Allow:

```json
{
  "scale_trip": 1,
  "derivative_policy": "numerical"
}
```

Then the Molpro calculator may request `FORCE`; Molpro can fall back to numerical gradients when analytic gradients are unavailable.

But the code must warn that:

- Sella will receive numerical forces;
- each Sella step is much more expensive;
- a finite-difference Hessian of those numerical forces would be a nested numerical differentiation and may be prohibitively expensive/noisy.

Require an explicit flag for this expensive combination:

```json
"allow_nested_numerical_derivatives": true
```

when applicable.

## 4.5 F12 geometry recipes must pin their own derivative policy

Original `ANL0-F12` does not require an F12 geometry.

A later recipe such as `ANL1-QZF` may.

For such a recipe, do **not** automatically inherit the energy-node `SCALE_TRIP=1`.

The recipe should contain something like:

```python
geometry = TaskSpec(
    method="CCSD(T)-F12b",
    derivative_method="DF-CCSD(T)-F12b",
    scale_trip=0,
    derivative_policy="analytic",
)
```

if that matches the published definition.

## 4.6 CABS and analytic F12 derivatives

Molpro's analytic F12 derivative implementation has specific CABS/RI restrictions.

Do not silently assume an energy input can be reused unchanged as a gradient input.

The Molpro renderer must have a distinct derivative mode with explicit:

```text
DF method
ansatz
RI basis
CABS settings
CABS singles behavior
```

and it must record whether a CABS singles correction was included numerically.

---

# 5. Reference and correlation policy

Keep these distinct:

```text
SCF reference
CC correlation form
```

They are not the same thing.

## 5.1 Ordinary closed-shell calculation

Default:

```text
SCF: RHF
CC: unrestricted-form CC where required by the pinned recipe
```

For a closed-shell case where the actual program's restricted and unrestricted equations are numerically equivalent, record what was requested and what implementation was used.

## 5.2 Ordinary open-shell calculation

Default:

```text
SCF: ROHF / restricted-spin open-shell HF
CC: unrestricted CC
```

This matches the established RUCCSD(T)-style convention.

## 5.3 Higher-order exceptions

Use explicit, recipe-pinned exceptions:

```text
open-shell CCSDT(Q):
    UHF + unrestricted CC
    MRCC

CCSDTQ(P):
    reference pinned by recipe/backend
    MRCC
```

Do not use:

```text
try ROHF
if crash:
    use UHF
```

A reference change is a method change.

## 5.4 Cross-program differences are valid

Do not impose a same-program denominator rule.

For example, established ANL-family work may form:

```text
UUCCSDT(Q)/DZ
-
RUCCSD(T)/DZ
```

where the two calculations are produced through different software paths.

The implementation should reproduce the recipe exactly and reuse an existing denominator when possible.

---

# 6. Custom ASE calculator architecture

Add:

```text
kinbot/ase_modules/calculators/
    factory.py
    gaussian.py      # existing
    qchem.py         # existing
    orca.py          # existing
    molpro.py        # new
    cfour.py         # new
    mrcc.py          # new
```

Do not move or rename existing calculators in the first implementation.

That minimizes merge conflicts.

## 6.1 Use the modern ASE GenericFileIO pattern for new calculators

Pattern new wrappers after current `orca.py`:

```python
class ProgramProfile(BaseProfile):
    ...

class ProgramTemplate(CalculatorTemplate):
    ...

class Program(GenericFileIOCalculator):
    ...
```

Benefits:

- isolated working directory;
- explicit input/output filenames;
- normal ASE `Atoms.calc` semantics;
- clean `energy` / `forces` capability exposure;
- ASE profile command handling;
- easier unit testing.

## 6.2 Calculator capability descriptor

Add a small immutable capability object.

Example:

```python
@dataclass(frozen=True)
class CalculatorCapabilities:
    energy: bool = True
    forces: bool = False
    native_hessian: bool = False
    dipole: bool = False
    rohf: bool = False
    uhf: bool = False
    f12: bool = False
    numerical_forces: bool = False
```

The factory returns both the calculator class and capabilities.

Do not advertise a property merely because a program can theoretically calculate it.

Advertise only what the KinBot wrapper correctly renders and parses.

## 6.3 Factory

Example:

```python
CALCULATORS = {
    "gaussian": GaussianSpec(...),
    "molpro": MolproSpec(...),
    "cfour": CFourSpec(...),
    "mrcc": MRCCSpec(...),
    "orca": OrcaSpec(...),
    "qchem": QChemSpec(...),
    "fairchem": FairChemSpec(...),
}
```

Interface:

```python
build_calculator(profile, directory, task) -> Calculator
capabilities(name) -> CalculatorCapabilities
```

This removes repeated backend mapping from `qc.py`.

---

# 7. Gaussian calculator

## 7.1 Keep current wrapper

Do not rewrite existing `gaussian.py` unless required.

It already supports:

```text
energy
forces
dipole
```

and works with Sella.

Only add:
- capability metadata if necessary;
- native-Hessian helper integration outside the class;
- perhaps a stable executable/version accessor.

Avoid converting it to `GenericFileIOCalculator` in the same PR.

That would create unnecessary conflict risk.

## 7.2 L2 geometry optimization

For the new profiled path:

```text
Gaussian is used only as energy/gradient calculator
Sella performs the geometry optimization
```

Do not use Gaussian `Opt` for the default profiled L2 path.

Default route should correspond to:

```text
B2PLYP/cc-pVTZ
EmpiricalDispersion=GD3BJ
Symmetry=None
SCF=XQC
Force
```

as required by the current calculator writer/property request.

## 7.3 TS optimization

Use:

```python
Sella(atoms, order=1, ...)
```

with the Gaussian calculator supplying forces.

## 7.4 HIR

Use the same L2 Gaussian calculator with:

```python
Constraints.fix_dihedral(...)
```

and Sella constrained optimization.

This preserves a single B2PLYP-D3(BJ)/cc-pVTZ surface for:
- stationary-point L2;
- HIR scans.

---

# 8. Molpro ASE calculator

This is the most important new backend.

## 8.1 Intended capabilities

Initial supported properties:

```text
energy
forces
```

Optional helper capability:

```text
native_hessian
```

through a separate frequency driver.

The calculator should support:

- RHF/ROHF/UHF reference rendering;
- CCSD(T);
- UCCSD(T);
- CCSD(T)-F12b/UCCSD(T)-F12b;
- frozen-core/all-electron modes;
- DKH settings;
- arbitrary Dunning basis names and explicit basis blocks;
- charge and spin;
- F12 auxiliary/RI settings;
- `SCALE_TRIP`;
- deterministic tagged output.

## 8.2 Stable output parsing

Do not grep “the last CCSD(T) energy.”

The input renderer should explicitly print a KinBot-tagged scalar after the requested method.

Conceptually:

```text
KB_FINAL_ENERGY = <resolved program variable>
```

and then print/table that value with a stable label.

The parser should:
1. require normal Molpro termination;
2. locate the expected KinBot marker;
3. parse exactly one final energy;
4. validate method metadata from the job spec.

## 8.3 Forces

If ASE requests `forces`, render:

```text
energy method
FORCE
```

Molpro will use analytic gradients where available and numerical gradients otherwise.

The parser converts Molpro gradients to ASE forces:

```text
force = -gradient
```

with correct units.

Unit tests must verify:
- sign;
- Hartree/Bohr -> eV/Å conversion;
- atom ordering.

## 8.4 Reference rendering

Normal open-shell tasks:

```text
restricted-spin HF/ROHF orbitals
unrestricted CC
```

Higher-order MRCC tasks do not go through this calculator.

## 8.5 F12 rendering

Energy node default:

```text
UCCSD(T)-F12b
SCALE_TRIP=1
```

Derivative node auto behavior:

```text
DF-UCCSD(T)-F12b
SCALE_TRIP=0
```

when analytic gradients are requested and supported.

Do not render generic `-F12` and later guess which energy was intended.

---

# 9. CFOUR ASE calculator

## 9.1 Initial scope

Keep it deliberately narrow.

Supported production tasks:

```text
energy-like DBOC property
closed-shell CCSDT(Q) energy
```

Optionally expose ordinary energies required for testing or user overrides.

Do not claim general gradient support in version 1 unless it is explicitly implemented and validated.

## 9.2 Fixed filenames

CFOUR commonly uses fixed working filenames such as `ZMAT`.

Therefore every calculator invocation must run in its own isolated ASE calculator directory.

Never execute two CFOUR tasks in the same directory.

## 9.3 DBOC

Treat DBOC as a named result/provenance field, not necessarily ASE's standard `energy`.

Possible calculator results:

```python
results["energy"] = ...
results["dboc"] = ...
```

If ASE disallows arbitrary result keys in a path, keep DBOC parsing in a backend-specific result helper called after calculator execution.

The ANL task node should explicitly request:

```text
property = dboc
```

## 9.4 Closed-shell CCSDT(Q)

Use CFOUR as requested.

Parse a stable final energy marker or CFOUR's documented method-specific line.

Test against real output fixtures.

---

# 10. MRCC ASE calculator

## 10.1 Initial scope

Start as energy-only:

```text
open-shell CCSDT(Q)
closed-shell CCSDTQ(P)
open-shell CCSDTQ(P)
```

This is enough for the ANL higher-order graph.

Do not implement optimization/forces simply for architectural symmetry.

## 10.2 Driver mode

Support a backend parameter:

```text
driver = direct
driver = molpro_interface
```

only if both are genuinely needed.

Recommended first implementation:

```text
direct MRCC
```

unless reproducing a specific literature value requires the Molpro-MRCC interface.

The task provenance must record the driver.

## 10.3 Reference

Recipe can force:

```text
UHF
```

for higher-order open-shell tasks.

This must not inherit the normal ROHF reference by accident.

---

# 11. ASE/Sella job runner

The current template system works, but adding multiple profiled backends to many templates would duplicate routing logic.

The minimal new-path solution is one generic runner.

Add:

```text
kinbot/ase_job.py
```

## 11.1 Job specification

KinBot writes:

```text
<job>.kinbot.json
```

containing:

```json
{
  "schema": 1,
  "job": "...",
  "profile": "l2",
  "calculator": "gaussian",
  "calculator_parameters": {},
  "atoms": {},
  "charge": 0,
  "multiplicity": 1,
  "task": "optimize",
  "stationary_order": 0,
  "constraints": [],
  "frequency_mode": "native_hessian",
  "sella": {},
  "database": "/abs/path/kinbot.db"
}
```

The scheduler Python file can be tiny:

```python
from kinbot.ase_job import run_job
run_job("job.kinbot.json")
```

This avoids one template per calculator per task.

## 11.2 Supported runner modes

Initial modes:

```text
energy
optimize
transition_state
hindered_rotor
frequency
```

Optional later:

```text
irc
```

Do not move current specialized native IRC logic into this runner unless necessary.

## 11.3 Optimization

For `N >= 3`:

```python
Sella(...)
```

For diatomics:

```python
BFGS(...)
```

matching existing behavior.

Atoms:
- one energy evaluation;
- no optimization;
- zero vibrational modes.

## 11.4 Sella is the optimizer, calculator supplies energy/forces

The profiled path must not embed a program-native geometry optimizer.

That produces consistent behavior across:
- Gaussian;
- Molpro;
- future calculator backends.

## 11.5 HIR constraints

Use current Sella `Constraints` semantics.

The generic runner accepts:

```json
{
  "kind": "dihedral",
  "atoms": [0, 1, 2, 3],
  "value": ...
}
```

Do not reproduce program-specific ModRedundant syntax for the profiled path.

---

# 12. Frequency and Hessian architecture

Optimization and frequency calculation should be separate concepts.

Add:

```text
kinbot/hessian.py
```

or a similarly small module.

## 12.1 Frequency modes

Support:

```text
auto
ase_forces
native_hessian
```

### `ase_forces`

Use current:

```python
ase.vibrations.Vibrations
```

which finite-differences calculator forces.

This already integrates correctly with KinBot's frequency projection logic.

### `native_hessian`

Run the program's fixed-geometry Hessian/frequency job, parse a Cartesian Hessian, and then pass that Hessian through the same KinBot mode-analysis code.

### `auto`

Choose the most sensible validated route for the calculator/method.

## 12.2 Default L1 UMA frequency mode

If an L1 Hessian is actually required:

```text
ase_forces
```

using the FairChem force calculator is appropriate.

Once L2 succeeds, L2 properties supersede those L1 properties for state counting.

## 12.3 Default Gaussian L2 frequency mode

Use:

```text
native_hessian
```

by default.

Reason:
- Gaussian can calculate an analytic B2PLYP Hessian;
- finite-differencing Gaussian gradients would cost roughly `6N` gradient evaluations;
- the native Hessian avoids unnecessary work.

Still pass the parsed Cartesian Hessian through KinBot's `get_frequencies(...)` machinery so:
- linear molecules;
- external modes;
- low frequencies;
- HIR projection

use one consistent convention.

## 12.4 Molpro conventional CCSD(T) Hessian

Use a method-specific policy.

Options:
1. Molpro native numerical frequency procedure;
2. ASE finite-difference of Molpro analytic gradients.

Choose the validated cheaper/stabler path in the recipe.

Do not blindly finite-difference already numerical gradients.

## 12.5 F12 Hessian warning

A scaled-F12 energy surface plus numerical gradient plus ASE finite-difference Hessian can become nested numerical differentiation.

The validator should reject or strongly gate:

```text
scaled F12 energy
+ numerical force
+ ASE finite-difference Hessian
```

unless the user explicitly allows it.

---

# 13. Theory profile model

Add:

```text
kinbot/theory.py
```

Example:

```python
@dataclass(frozen=True)
class TheoryProfile:
    name: str
    calculator: str
    method: str = ""
    basis: str = ""
    command: str = ""
    calculator_kwargs: Mapping[str, Any] = field(default_factory=dict)
    optimizer: str = "sella"
    frequency_mode: str = "auto"
    label: str = ""
```

Profiles:

```text
l1
l2
scan
semi_emp
aie
barrierless
barrierless_high
vts
```

Keep specialized contexts because current global-backend semantics would otherwise accidentally route them to FairChem.

---

# 14. Default preset

Proposed:

```python
THEORY_PRESETS["uma-b2plyp-anl"] = {
    "l1": {
        "calculator": "fairchem",
        "task_name": "omol",
    },
    "l2": {
        "calculator": "gaussian",
        "method": "B2PLYP",
        "basis": "cc-pVTZ",
        "calculator_kwargs": {
            "EmpiricalDispersion": "GD3BJ",
            "Symm": "None",
            "scf": "xqc",
        },
        "optimizer": "sella",
        "frequency_mode": "native_hessian",
        "label": "B2PLYP-D3(BJ)/cc-pVTZ",
    },
}
```

Do not include a hard-coded UMA checkpoint path.

That remains user/site input.

---

# 15. Specialized contexts

The router must prevent these from accidentally becoming UMA.

## 15.1 Semiempirical conformer search

If enabled, preserve a dedicated semiempirical profile.

Do not ask FairChem to interpret `AM1`.

## 15.2 AIE

Current AIE logic uses CBS-QB3/Gaussian semantics.

Route it explicitly to its compatible backend.

## 15.3 Scan/MP2 paths

Current `scan_method` / `scan_basis` semantics may imply conventional QC.

Give scan a separate profile.

## 15.4 Barrierless saddle

Preserve dedicated barrierless method/basis profiles.

## 15.5 VRC-TST

Do not convert existing VRC-TST scan or long-range potential calculations to UMA merely because ordinary L1 is UMA.

VRC remains its own theory context.

## 15.6 IRC

Do not rewrite IRC in the first implementation unless necessary.

If the current L1 backend-specific IRC path works and UMA supports the intended IRC workflow, preserve it.

For the new mixed path, TS refinement can be B2PLYP while IRC may remain a dedicated profile.

Make that explicit rather than implicit.

---

# 16. Persistent job metadata

Mixed backend runs require persistent provenance.

## 16.1 Before result exists

Maintain:

```text
kinbot_jobs.json
```

or a similarly small manifest.

Example:

```json
{
  "123_well_high": {
    "profile": "l2",
    "calculator": "gaussian",
    "method": "B2PLYP",
    "basis": "cc-pVTZ",
    "label": "B2PLYP-D3(BJ)/cc-pVTZ"
  }
}
```

Write atomically.

## 16.2 After result exists

Store the same metadata in ASE DB row `data`.

Example:

```json
{
  "kinbot_profile": "l2",
  "kinbot_calculator": "gaussian",
  "kinbot_method": "B2PLYP",
  "kinbot_basis": "cc-pVTZ",
  "kinbot_theory_label": "B2PLYP-D3(BJ)/cc-pVTZ",
  "kinbot_frequency_mode": "native_hessian"
}
```

## 16.3 Resolver order

`backend_for_job(job)`:

```text
1. latest matching DB row metadata
2. job manifest
3. deterministic legacy inference
4. global legacy backend only as final fallback
```

---

# 17. Flag-guarded integration strategy

This is the most important software-maintainability decision.

## 17.1 Legacy path remains untouched

Default:

```python
self.profiled_theory = False
```

Existing functions continue through existing bodies.

## 17.2 New path delegates at function entry

Conceptually:

```python
def qc_opt(...):
    if self.profiled_theory:
        return self.profiled.qc_opt(...)
    # existing function body below, unchanged
```

Do the same only where necessary.

Recommended delegated functions:

```text
get_qc_arguments
qc_conf
qc_opt
qc_hir
qc_freq
check_qc
read_qc_hess
hessian_is_massweighted
```

Do not mass-refactor all 1500+ lines of `qc.py`.

## 17.3 New profiled implementation

Add:

```text
kinbot/profiled_qc.py
```

This module owns:
- profile selection;
- generic ASE job-spec generation;
- job registration;
- result checking for profiled jobs;
- Hessian retrieval for profiled jobs.

This keeps `qc.py` edits tiny and merge-friendly.

---

# 18. ANL composite package

Use a package rather than continuing to grow the old `molpro.py`.

```text
kinbot/anl/
    __init__.py
    model.py
    recipes.py
    workflow.py
```

Avoid separate program parsers here.

Program parsing belongs inside the calculator interfaces.

## 18.1 `model.py`

Define:

```text
TaskSpec
ExpressionTerm
ExpressionSpec
CompositeRecipe
TaskResult
CompositeResult
```

## 18.2 `recipes.py`

Contains only declarative method definitions.

No scheduler logic.

## 18.3 `workflow.py`

Responsible for:
- dependency resolution;
- task de-duplication;
- task signatures;
- restart;
- expression evaluation;
- provenance;
- final result file.

---

# 19. ANL task graph

Example node IDs:

```text
geometry
harmonic

reference_qz
reference_5z
reference_6z

f12_tz
f12_qz
f12_5z

cv_tz_fc
cv_tz_ae
cv_qz_fc
cv_qz_ae

rel_nr
rel_dkh

dboc

ccsdtq_dz
ccsdtq_tz
ccsdtqp_dz

vpt2_harm
vpt2_anh

spin_orbit
```

Not every recipe uses every node.

## 19.1 De-duplicate exact nodes

Task identity should hash:

```text
geometry hash
charge
multiplicity
calculator
program version policy
method
basis
reference
correlation form
frozen-core mode
relativistic mode
F12 ansatz
scale_trip
```

If two expressions require the same task:
- compute it once;
- reuse it.

---

# 20. Cross-code higher-order expressions

Represent literally.

Example:

```python
higher_order = ExpressionSpec(
    terms=(
        (+1.0, "uuccsdtq_dz"),
        (-1.0, "ruccsdt_dz"),
    )
)
```

For an ANL1-like extension:

```text
+ CCSDTQ(P)/DZ
- CCSDT(Q)/DZ
+ CCSDT(Q)/TZ
- CCSD(T)/TZ
```

Each absolute term can come from the program dictated by the recipe.

No artificial same-program correction rule.

---

# 21. Program ownership defaults

## Molpro

Default for:
- ANL geometry where canonical recipe calls for conventional CCSD(T);
- conventional harmonic frequencies;
- conventional CCSD(T) basis calculations;
- all ordinary F12b energy terms;
- core-valence;
- scalar relativistic;
- ordinary RUCCSD(T) denominator/reference nodes.

## CFOUR

Default for:
- DBOC;
- closed-shell/singlet CCSDT(Q).

## MRCC

Default for:
- open-shell CCSDT(Q);
- all CCSDTQ(P).

## Gaussian

Default for:
- VPT2;
- KinBot L2 B2PLYP-D3(BJ)/cc-pVTZ.

---

# 22. ANL geometry versus KinBot L2 geometry

They are separate.

KinBot:

```text
UMA L1
↓
B2PLYP-D3(BJ)/cc-pVTZ L2
```

ANL may then do:

```text
Molpro CCSD(T)/TZ geometry
```

or another recipe-defined high-level geometry.

The L2 geometry is the input/start structure to ANL.

Do not redefine KinBot L2 to mean ANL geometry.

---

# 23. Final energy abstraction

Add:

```text
kinbot/energy.py
```

Define:

```python
@dataclass
class FinalEnergy:
    electronic_hartree: float | None
    zero_k_hartree: float
    source: str
    method_label: str
```

Rules:

## L2 only

```text
E0 = E_L2 + ZPE_L2
```

## Existing legacy L3

```text
E0 = E_L3(electronic) + ZPE_L2
```

Preserve current behavior.

## ANL

```text
E0 = E_ANL,0K
```

Do **not** add L2 ZPE.

Use this resolver in both:
- PES;
- direct MESS generation.

---

# 24. State-counting properties versus final energy

Default MESS/statistical-mechanics behavior:

```text
energy zero:
    ANL E0

geometry:
    B2PLYP-D3(BJ)/cc-pVTZ L2

harmonic frequencies:
    B2PLYP-D3(BJ)/cc-pVTZ L2

hindered rotors:
    B2PLYP-D3(BJ)/cc-pVTZ L2
```

This is deliberate.

Do not mix the ANL thermochemical ZPE bookkeeping with the MESS HIR representation without a separate, explicitly defined thermochemical correction model.

---

# 25. Hindered rotors

All default HIR scans should use L2.

For the default preset:

```text
Gaussian B2PLYP-D3(BJ)/cc-pVTZ
```

The HIR zero-angle energy and constrained scan energies must be on the same level.

If HIR finds a lower minimum:
1. invalidate the accepted L2 result;
2. reoptimize at L2;
3. regenerate L2 Hessian;
4. invalidate dependent HIR points if required;
5. rerun HIR;
6. only then launch/reuse ANL.

Do not run ANL on a superseded L2 conformer.

---

# 26. Atoms, diatomics, and linear species

## Atom

```text
no optimization
no frequency
ZPE = 0
no HIR
```

ANL electronic corrections and SO still may apply.

## Diatomic

Use BFGS rather than Sella if that is the established robust current path.

One vibrational mode.

## Linear polyatomic

Mode count:

```text
3N - 5
```

Do not assume `3N-6`.

Reuse current KinBot linear/near-linear safeguards.

---

# 27. TS validation

At each refinement stage:

## L2 TS
- exactly one chemically meaningful imaginary mode;
- use current small-mode tolerance behavior.

## ANL geometry TS
- preserve intended connectivity;
- verify one intended imaginary mode;
- compare mode displacement to reaction center when possible.

Do not fan out expensive ANL single points until the high-level stationary point is validated.

---

# 28. Open-shell and multireference edge cases

## 28.1 Radical ordinary CC

Use:

```text
ROHF/restricted-spin reference
unrestricted CC
```

where supported.

## 28.2 Higher-order open-shell

Use recipe-pinned UHF + unrestricted CC when required.

## 28.3 Broken-symmetry singlets

Do not silently treat every multiplicity-1 state as safe RHF.

If an open-shell singlet/biradical is detected:
- flag;
- require explicit treatment;
- do not label an ad hoc broken-symmetry calculation exact ANL.

## 28.4 Diagnostics

Record when available:
- T1;
- D1;
- `<S^2>`;
- HOE magnitude;
- CC convergence quality.

Large HOE should generate a multireference warning.

---

# 29. Fluorinated systems and heavier elements

The Ram et al. work demonstrates why basis/component definitions must be recipe data rather than hard-coded “CHNO ANL.”

For F-containing species:
- diffuse augmentation rules matter;
- CV basis choices matter;
- DKH basis choices matter;
- very large F12/conventional basis calculations may be feasible only for smaller fragments.

Do not pretend an original CHNO ANL recipe automatically has the same benchmark uncertainty for F systems.

A fluorinated extension should have its own recipe/provenance label.

For elements lacking required:
- orbital basis;
- F12 basis;
- CABS/OPTRI;
- CV basis;
- DKH basis;

fail preflight unless the user supplies an explicit custom recipe.

---

# 30. Basis handling

Do not rely on ambiguous shorthand such as `aug'` without defining exactly what it means.

Represent basis construction in recipe data.

For example:

```text
H: augmentation rule A
C/N/O: augmentation rule B
F: recipe-specific rule
```

Generated Molpro input should be snapshot-tested.

---

# 31. Spin-orbit

SO is not a generic calculator fallback.

Represent it as a task/provider with status:

```text
known_zero
table_value
calculated
manual_override
missing_required
```

Key it by electronic state, not formula alone.

Never treat “lookup not found” as zero for a state where SO can matter.

---

# 32. Geometry identity and hashes

Every derived ANL node records the geometry hash.

If the geometry changes:
- invalidate all dependent energy tasks;
- invalidate harmonic task;
- invalidate geometry-dependent VPT2 task.

For minima/TS/vdW structures also store:
- connectivity signature;
- fragment count;
- selected conformer identity.

---

# 33. vdW complexes

Check after every high-level optimization:
- fragment identities retained;
- no unintended dissociation;
- no collapse into covalent minimum.

Very low intermolecular modes should be flagged rather than silently deleted.

---

# 34. Barrierless/VRC-TST channels

Do not run the full stationary-point ANL composite at each VRC geometry.

Use ANL for:
- separated reactants;
- separated products;
- actual stationary wells/TSs where appropriate.

Keep VRC interaction potentials on their dedicated VRC theory path.

---

# 35. SCF retry policy

Permitted automatic retries that do not change theory:
- increased iterations;
- damping;
- level shift;
- projected orbitals;
- basis continuation;
- orbital ordering;
- tighter/adjusted convergence.

Not permitted silently:
- ROHF -> UHF;
- multiplicity change;
- different electronic state;
- different method.

---

# 36. Coupled-cluster retry policy

Permitted:
- restart amplitudes;
- iteration increase;
- stable integral/threshold options consistent with method.

Not permitted:
- dropping `(Q)`;
- dropping `(P)`;
- substituting CCSD(T);
- changing frozen-core convention;
- changing reference.

Incomplete is incomplete.

---

# 37. Frequency failures

## Gaussian L2
If native Hessian fails:
- retry same theory with safe numerical/settings changes;
- optionally allow ASE-force finite difference as an explicit fallback.

Do not silently change functional/basis.

## Molpro conventional
Allow validated numerical Hessian/gradient fallback.

## F12
Avoid nested numerical differentiation by default.

## VPT2
Resonance/pathology failure should make exact recipe incomplete unless a documented fallback belongs to that recipe.

---

# 38. Conformer behavior

Recommended sequence:

```text
UMA conformer generation/ranking
↓
selected conformer(s)
↓
B2PLYP L2 reoptimization
↓
HIR can discover lower B2PLYP structure
↓
final accepted L2 conformer
↓
ANL
```

First implementation should not launch an ANL conformer search.

Optional future:
- ANL/F12 re-ranking of a small L2 conformer set.

---

# 39. Multi-conformer TST

Preserve current behavior that can disable 1-D HIR when multi-conformer TST is active.

The theory preset must define methods, not override that workflow decision.

---

# 40. Resource profiles

Keep resources separate from theory.

Example:

```text
uma_l1
gaussian_l2
molpro_geometry
molpro_hessian
molpro_f12_small
molpro_f12_large
cfour_dboc
cfour_ccsdtq
mrcc_ccsdtq
mrcc_ccsdtqp
gaussian_vpt2
```

Resource settings:
- cores;
- memory;
- walltime;
- scratch;
- queue.

Changing them does not alter protocol exactness.

---

# 41. Local mode

For:

```text
queuing = local
```

new paths should:
- parse existing output;
- update DB/result manifests;
- report missing jobs;
- not submit anything.

This is especially useful for manually run CCSDTQ(P) jobs.

---

# 42. File organization

Recommended final production structure:

```text
kinbot/
├── theory.py
├── profiled_qc.py
├── ase_job.py
├── hessian.py
├── energy.py
│
├── anl/
│   ├── __init__.py
│   ├── model.py
│   ├── recipes.py
│   └── workflow.py
│
└── ase_modules/
    └── calculators/
        ├── __init__.py
        ├── factory.py
        ├── gaussian.py       # existing, minimal edit
        ├── qchem.py          # existing
        ├── orca.py           # existing
        ├── molpro.py         # new
        ├── cfour.py          # new
        └── mrcc.py           # new
```

Do not add many new templates.

The profiled path uses one generic JSON-driven ASE runner.

Legacy templates stay exactly where they are.

---

# 43. Estimated production diff

These are engineering estimates, not exact promises.

## New files

| File | Estimated LOC |
|---|---:|
| `theory.py` | 180–230 |
| `profiled_qc.py` | 220–300 |
| `ase_job.py` | 220–300 |
| `hessian.py` | 120–180 |
| `energy.py` | 70–110 |
| `calculators/factory.py` | 80–120 |
| `calculators/molpro.py` | 220–300 |
| `calculators/cfour.py` | 130–190 |
| `calculators/mrcc.py` | 140–200 |
| `anl/model.py` | 120–170 |
| `anl/recipes.py` | 220–320 |
| `anl/workflow.py` | 240–340 |
| **new production total** | **1,960–2,760** |

A disciplined implementation should target the lower-middle of that range.

## Existing hot files

| File | Estimated edit |
|---|---:|
| `parameters.py` | +80–130 |
| `qc.py` | +70–130, mostly guards/delegation |
| `optimize.py` | +15–30 |
| `frequencies.py` | +10–30 |
| `calculation.py` | +0–15 |
| `pes.py` | +20–40 |
| `mess.py` | +15–30 |
| `utils.py` | +10–25 |
| `gaussian.py` | +0–25 |
| `pyproject.toml` | 0–5 |
| **existing-file churn** | **~220–460 LOC** |

The key objective is not merely low total LOC.

It is:

> keep changes to existing files small enough that current upstream branches can merge with minimal conflict.

## Tests

Expect:

```text
1,500–2,500 LOC
```

of tests/fixtures over time.

Tests do not count as undesirable code bloat here; they are the protection against corrupt composite energies.

---

# 44. Exact targeted edits to existing files

## 44.1 `parameters.py`

Only:
1. new profile/composite defaults;
2. preset application before explicit input overrides;
3. validation;
4. automatic `profiled_theory=1` activation.

Do not reorganize the rest of the defaults dictionary.

## 44.2 `qc.py`

At `__init__`:
- create `ProfiledQuantumChemistry` only when flag is active.

At selected methods:
- one early guard/delegate.

Example:

```python
def qc_hir(...):
    if self.profiled_theory:
        return self.profiled.qc_hir(...)
    # current body unchanged
```

Likewise for:
- `qc_conf`;
- `qc_opt`;
- `qc_freq`;
- relevant status/Hessian methods.

Do not rewrite every old backend block.

## 44.3 `optimize.py`

Only:
- pass job identity into mass-weight/Hessian helper if needed;
- dispatch composite workflow at current L3 hook;
- use final-energy/provenance hook where appropriate.

Preserve current result-publication ordering.

## 44.4 `pes.py`

Replace only final energy acquisition with the shared resolver under the new flag/composite path.

Leave legacy path intact.

## 44.5 `mess.py`

Use:
- profile human-readable L2 label;
- shared final energy resolver.

Do not redesign MESS generation.

---

# 45. Why this is less conflict-prone than refactoring `qc.py`

Current master is active and has recently merged:
- FairChem support;
- HIR fixes;
- coherent result publication;
- frequency/linear-molecule fixes.

Those are exactly the regions likely to continue changing.

A large rewrite of `qc.py` would repeatedly conflict.

A guarded delegation path means:
- current legacy implementation stays intact;
- new functionality lives in new modules;
- upstream fixes to the old path can be merged independently.

---

# 46. Composite result files

Per stationary point:

```text
anl/<species>/
    manifest.json
    result.json
    geometry/
    harmonic/
    reference/
    f12/
    cv/
    relativistic/
    dboc/
    higher_order/
    vpt2/
```

`manifest.json` contains task state.

`result.json` contains final composite result.

Write both atomically.

---

# 47. Task states

Use:

```text
pending
prepared
submitted
running
complete
failed
invalidated
skipped
```

A file's existence is not proof of completion.

A completed task must pass:
- signature match;
- parser validation;
- normal termination;
- method/reference validation.

---

# 48. Result provenance

Store:
- KinBot git SHA;
- recipe version;
- calculator backend;
- program version;
- executable command;
- method;
- basis;
- reference;
- CC form;
- frozen-core flag;
- F12 ansatz;
- scale_trip;
- geometry hash;
- cores/memory;
- raw absolute energy;
- parser evidence;
- timestamps.

For UMA:
- FairChem version;
- model checkpoint path or hash/identifier;
- task name.

---

# 49. Validation tests — calculator layer

## Molpro
Fixtures:
- RHF CCSD(T);
- ROHF/RUCCSD(T);
- F12b unscaled;
- F12b `SCALE_TRIP=1`;
- analytic derivative mode;
- numerical FORCE mode;
- CV;
- DKH.

Validate energy and force units.

## CFOUR
Fixtures:
- DBOC;
- closed-shell CCSDT(Q).

## MRCC
Fixtures:
- open-shell CCSDT(Q);
- closed-shell CCSDTQ(P);
- open-shell CCSDTQ(P).

## Gaussian
Reuse existing tests and add:
- B2PLYP-D3(BJ) gradient;
- native Hessian parse.

---

# 50. Validation tests — optimizer layer

Test one simple surface with a mock calculator:
- minimum Sella;
- TS Sella;
- dihedral-constrained Sella;
- diatomic BFGS;
- atom shortcut.

Then real small-program smoke tests where software is available.

---

# 51. Validation tests — mixed profile workflow

Must prove:

```text
L1 conformer -> FairChem
L2 optimization -> Gaussian
L2 Hessian -> Gaussian
HIR -> Gaussian
restart -> all jobs still resolve correct backend
```

A fresh Python process must successfully resolve old jobs from persistent metadata.

---

# 52. Validation tests — ANL arithmetic

Use synthetic values.

Test:
- ANL0 CBS;
- ANL0-F12 substitution;
- ANL1 CBS;
- CV;
- DKH;
- DBOC;
- cross-code HOE;
- CCSDTQ(P) extension;
- harmonic ZPE;
- VPT2 delta;
- final E_electronic;
- final E0.

Include a regression showing cross-program subtraction is allowed.

---

# 53. Validation against literature

Use small reference species where published component tables exist.

Suggested progression:

```text
H2
H2O
CH4
CH3
OH
O2
HO2
H2O2
C2H6
C2H5
```

For every species compare:
- geometry where applicable;
- harmonic ZPE;
- each absolute high-level node if published;
- each correction;
- final composite.

Do not validate only the final sum.

---

# 54. Default configuration example

```json
{
  "theory_preset": "uma-b2plyp-anl",
  "fc_model_path": "/path/to/uma.pt",

  "conformer_search": 1,
  "high_level": 1,
  "rotor_scan": 1,

  "composite_method": "ANL0-F12"
}
```

Equivalent expanded conceptual configuration:

```json
{
  "profiled_theory": 1,

  "l1_profile": {
    "calculator": "fairchem",
    "model_path": "/path/to/uma.pt",
    "task_name": "omol"
  },

  "l2_profile": {
    "calculator": "gaussian",
    "method": "B2PLYP",
    "basis": "cc-pVTZ",
    "calculator_kwargs": {
      "EmpiricalDispersion": "GD3BJ"
    },
    "optimizer": "sella",
    "frequency_mode": "native_hessian"
  },

  "composite_method": "ANL0-F12"
}
```

---

# 55. Custom configuration example

```json
{
  "profiled_theory": 1,

  "l1_profile": {
    "calculator": "gaussian",
    "method": "M06-2X",
    "basis": "6-31+G(d,p)"
  },

  "l2_profile": {
    "calculator": "gaussian",
    "method": "B2PLYP",
    "basis": "cc-pVQZ",
    "calculator_kwargs": {
      "EmpiricalDispersion": "GD3BJ"
    }
  },

  "composite_method": "ANL0-F12",

  "l3_overrides": {
    "f12_qz": {
      "scale_trip": 1
    }
  }
}
```

This must automatically mark the result custom if the L2 or L3 change affects the named canonical recipe.

---

# 56. Preflight validation

Before submitting an expensive composite:

Check:
- all requested executables exist;
- FairChem model exists;
- calculator supports requested property;
- basis exists for all elements;
- F12/CABS/OPTRI availability;
- charge/multiplicity/electron count;
- recipe element domain;
- derivative policy is viable;
- no forbidden nested numerical derivative combination;
- memory/scratch paths;
- job manifest writable.

Fail before cluster time is consumed.

---

# 57. Method-specific edge cases that should be encoded

## 57.1 F12 energy vs gradient mismatch

Never let:

```text
energy = scaled F12b
gradient = unscaled analytic F12b
```

pass as one exact ASE potential surface without documenting the inconsistency.

For an optimization, energy and gradient must correspond to the same declared task policy.

If `scale_trip=1` requires numerical forces:
- compute force from that same scaled energy.

If using analytic unscaled F12 derivative:
- the task's energy should be the corresponding unscaled derivative-method energy during Sella.

A final scaled F12 single point can then be evaluated at the optimized geometry as a separate node.

This separation is theoretically cleaner.

## 57.2 Sella convergence with numerical gradients

Use tighter/appropriate noise-aware thresholds only if validated.

Do not loosen convergence silently just because a backend is expensive.

## 57.3 Hessian at geometry level different from final energy

This is normal in composite methods.

Provenance must say:
- geometry/hessian level;
- single-point final energy level.

---

# 58. Recommended way to handle future ANL1-QZF / extended recipes

Do not add arbitrary branches to `workflow.py`.

Add a new declarative recipe.

Example conceptual:

```python
ANL1_QZF = CompositeRecipe(
    geometry=...,
    harmonic=...,
    tasks=(...),
    expressions=(...),
)
```

A Ram-style feasibility ladder can provide conditional task alternatives:

```text
try preferred highest recipe tier if feasibility predicate passes
else fall back to pinned next tier
```

But the selected tier must be written into result provenance.

Do not call a fallback result exact highest-tier theory.

---

# 59. Feasibility predicates

For extended recipes only, allow rules based on:
- electron count;
- correlated electron count;
- basis-function estimate;
- multiplicity;
- estimated memory;
- configured walltime;
- user hard caps.

Never base a theoretical downgrade on “the job crashed once.”

A crash is not a scientific feasibility criterion.

---

# 60. Branch and implementation instructions

The following is the recommended exact development sequence.

---

## Step 0 — create/update your fork clone

If already cloned:

```bash
cd KinBot
git remote -v
```

Ensure upstream exists:

```bash
git remote add upstream https://github.com/zadorlab/KinBot.git
```

if needed.

Fetch current upstream:

```bash
git fetch upstream --prune
git fetch origin --prune
```

---

## Step 1 — branch from current upstream master

```bash
git switch master
git merge --ff-only upstream/master
```

Record baseline:

```bash
git rev-parse HEAD
git log -1 --oneline
git status
```

As of this review GitHub reported master at:

```text
924faba
```

but use whatever SHA `upstream/master` actually has when you start.

Create branch:

```bash
git switch -c feature/profiled-ase-anl
git push -u origin feature/profiled-ase-anl
```

---

## Step 2 — create an optional worktree

Useful if you want your normal KinBot checkout untouched:

```bash
cd ..
git -C KinBot worktree add KinBot-profiled feature/profiled-ase-anl
cd KinBot-profiled
```

If already on the branch in the main checkout, skip this.

---

## Step 3 — create clean Python environment

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[fc]'
python -m pip install pytest
```

Sanity:

```bash
python -c 'import kinbot; print(kinbot.__file__)'
python -c 'import ase; print(ase.__version__)'
python -c 'import sella; print(sella.__version__)'
python -m compileall kinbot
```

Run existing tests:

```bash
pytest -q
```

Record baseline failures, if any.

---

## Step 4 — verify external programs independently

```bash
which g16
which formchk
which molpro
which xcfour
which dmrcc
```

Verify UMA checkpoint:

```bash
test -r /path/to/uma.pt && echo "UMA checkpoint OK"
```

Do not commit executable paths.

---

## Step 5 — add baseline regression tests before code

Create tests asserting old behavior with:

```text
profiled_theory = 0
```

Specifically:
- Gaussian remains Gaussian;
- old high_level still uses old backend;
- old HIR route unchanged;
- old `check_qc` semantics unchanged;
- old Hessian reader unchanged.

Commit:

```bash
git add tests
git commit -m "test(qc): lock legacy QC behavior before profiled theory"
```

---

## Step 6 — add calculator factory and capability model

Create:

```text
kinbot/ase_modules/calculators/factory.py
```

Initially register existing:
- Gaussian;
- Q-Chem;
- ORCA;
- FairChem adapter.

Do not add new program wrappers yet.

Test factory only.

Commit:

```bash
git add kinbot/ase_modules/calculators tests
git commit -m "refactor(ase): add calculator capability factory"
```

---

## Step 7 — add theory profiles and preset parsing

Create:

```text
kinbot/theory.py
```

Modify only necessary parameter lines.

Implement:
- `TheoryProfile`;
- default preset;
- user override merge;
- activation flag.

Test precedence.

Commit:

```bash
git add kinbot/theory.py kinbot/parameters.py tests
git commit -m "feat(theory): add opt-in L1/L2 theory profiles"
```

---

## Step 8 — add generic ASE/Sella runner

Create:

```text
kinbot/ase_job.py
```

Start with mock-calculator tests.

Implement:
- atom;
- minimum;
- TS;
- constrained HIR;
- DB write;
- metadata write.

Do not wire into `qc.py` yet.

Commit:

```bash
git add kinbot/ase_job.py tests
git commit -m "feat(ase): add generic profiled Sella job runner"
```

---

## Step 9 — add guarded profiled QC adapter

Create:

```text
kinbot/profiled_qc.py
```

Then add only early delegation guards to `qc.py`.

First wire:
- `qc_conf`;
- `qc_opt`;
- `qc_hir`.

Keep legacy bodies unchanged below guards.

Commit:

```bash
git add kinbot/profiled_qc.py kinbot/qc.py tests
git commit -m "feat(qc): route opt-in profile jobs through ASE runner"
```

---

## Step 10 — persist per-job metadata

Add:
- atomic job manifest;
- DB row metadata.

Implement resolver.

Then guard/delegate:
- `check_qc`;
- `read_qc_hess`;
- `hessian_is_massweighted`.

Change massweight signature only where the new path needs it, or support:

```python
hessian_is_massweighted(job=None)
```

to avoid breaking legacy callers.

Commit:

```bash
git commit -am "fix(qc): persist and resolve calculator identity per job"
```

---

## Step 11 — implement L1 UMA + L2 B2PLYP smoke workflow

Use a small molecule.

First:

```text
conformer_search = 1
high_level = 1
rotor_scan = 0
```

Verify:
- L1 FairChem;
- L2 Gaussian via Sella;
- coherent L2 DB row.

Then restart from a new process and verify result status.

Commit fixes before moving on.

---

## Step 12 — add native Hessian abstraction

Create:

```text
kinbot/hessian.py
```

First implement Gaussian native Hessian.

Pass parsed Hessian through KinBot frequency projection.

Wire profiled `qc_freq` / L2 finalization.

Test:
- nonlinear;
- linear;
- TS.

Commit:

```bash
git add kinbot/hessian.py kinbot/qc.py kinbot/frequencies.py tests
git commit -m "feat(freq): add profile-aware native Hessian driver"
```

---

## Step 13 — validate L2 hindered rotors

Turn on:

```text
rotor_scan = 1
```

Verify every `hir/` calculation is:
- Gaussian;
- B2PLYP;
- cc-pVTZ;
- GD3BJ;
- Sella constrained.

Test lower-minimum restart.

Commit:

```bash
git commit -am "fix(hir): keep profiled rotor scans on L2 surface"
```

---

## Step 14 — audit specialized contexts

Before ANL, verify:
- semiempirical conformer;
- scan/MP2;
- AIE;
- barrierless;
- VRC.

Add profile mappings only where needed.

Commit:

```bash
git commit -am "fix(theory): preserve specialized QC contexts under profiles"
```

---

## Step 15 — implement Molpro calculator

Create:

```text
kinbot/ase_modules/calculators/molpro.py
```

Implement in this order:
1. energy-only RHF;
2. ROHF/RUCCSD(T);
3. forces;
4. F12b energy;
5. `SCALE_TRIP`;
6. DF-F12 derivative mode;
7. CV;
8. DKH.

Use real fixture outputs.

Commit:

```bash
git add kinbot/ase_modules/calculators/molpro.py tests
git commit -m "feat(ase): add Molpro calculator"
```

---

## Step 16 — implement CFOUR calculator

Start with:
- DBOC;
- closed-shell CCSDT(Q).

Commit:

```bash
git add kinbot/ase_modules/calculators/cfour.py tests
git commit -m "feat(ase): add narrow CFOUR calculator"
```

---

## Step 17 — implement MRCC calculator

Start with:
- open-shell CCSDT(Q);
- CCSDTQ(P).

Commit:

```bash
git add kinbot/ase_modules/calculators/mrcc.py tests
git commit -m "feat(ase): add narrow MRCC calculator"
```

---

## Step 18 — add ANL data model and recipes

Create:

```text
kinbot/anl/model.py
kinbot/anl/recipes.py
```

No execution yet.

Test all equations with synthetic values.

Include:
- cross-program terms;
- F12b normalization;
- scale_trip policy.

Commit:

```bash
git add kinbot/anl tests
git commit -m "feat(anl): add declarative composite recipes"
```

---

## Step 19 — add ANL workflow/restart engine

Create:

```text
kinbot/anl/workflow.py
```

Implement:
- task DAG;
- signature;
- deduplication;
- restart;
- task directories;
- atomic manifest/result;
- expression evaluation.

Commit:

```bash
git add kinbot/anl tests
git commit -m "feat(anl): execute and resume composite task graph"
```

---

## Step 20 — connect ANL at existing L3 hook

Edit `optimize.py` only around current `L3_calc` block.

Pseudo:

```python
if self.par.get("composite_method"):
    run_or_check_composite(...)
elif self.par["L3_calc"]:
    # current legacy Molpro/ORCA code unchanged
```

Commit:

```bash
git add kinbot/optimize.py tests
git commit -m "feat(anl): connect composite provider at L3 boundary"
```

---

## Step 21 — add final energy resolver

Create:

```text
kinbot/energy.py
```

Patch:
- PES;
- MESS.

Explicit test:

```text
ANL E0 must not receive + L2 ZPE
```

Commit:

```bash
git add kinbot/energy.py kinbot/pes.py kinbot/mess.py tests
git commit -m "feat(energy): unify L2 L3 and ANL zero-K energies"
```

---

## Step 22 — literature regression

Do not call the feature complete until:
- at least one closed-shell original ANL result;
- one open-shell original ANL result;
- one ANL0-F12 result;
- one ANL1 higher-order result

have been reproduced component-by-component.

Then add fluorinated/later-family tests.

---

# 61. Rebase strategy while developing

Because upstream master is currently active:

```bash
git fetch upstream
git rebase upstream/master
```

regularly, especially before modifying:
- `qc.py`;
- `parameters.py`;
- `optimize.py`;
- `pes.py`;
- `mess.py`.

Resolve conflicts while commits are still narrow.

Do not wait until the whole ANL implementation is one giant branch diff.

---

# 62. Recommended PR split

Even if development begins in one branch, the clean upstreamable sequence is:

## PR 1
```text
calculator factory
theory profiles
generic ASE/Sella runner
mixed L1/L2
```

## PR 2
```text
Molpro/CFOUR/MRCC calculators
native Hessian support
```

## PR 3
```text
ANL composite recipes/workflow
final-energy resolver
```

This makes review and conflict resolution much easier.

If you do not intend to upstream immediately, the same commit boundaries still make debugging safer.

---

# 63. Things not to do

Do not:

1. replace global `qc` behavior for all users;
2. rewrite the entire `qc.py`;
3. add one template for every program × task combination;
4. perform geometry optimization inside Molpro/Gaussian while claiming Sella owns optimization;
5. expose `forces` from MRCC/CFOUR before actually parsing them correctly;
6. call all F12 variants simply `-F12`;
7. apply `SCALE_TRIP=1` to analytic F12 gradients silently;
8. subtract L2 ZPE from/add L2 ZPE to make ANL fit the old L3 API;
9. require same-program HOE terms;
10. rerun an existing denominator only to make software names match;
11. silently switch ROHF to UHF;
12. silently downgrade CC rank;
13. silently use a lower basis;
14. label a modified recipe canonical ANL;
15. run ANL on every VRC geometry;
16. launch high-cost ANL nodes before stationary-point validation;
17. infer a completed method solely from a filename;
18. rely on one global backend after mixed profiles are enabled.

---

# 64. Final recommended architecture

The final conceptual layering should be:

```text
KinBot chemistry workflow
        |
        +-- legacy QC path ------------------------------+
        |      active when profiled_theory == 0         |
        |                                                |
        +-- profiled QC adapter                          |
               |                                         |
               +-- TheoryProfile                         |
               +-- calculator factory                    |
               +-- generic ASE/Sella runner              |
               +-- Hessian/frequency driver              |
               |                                         |
               +-- FairChem / UMA                        |
               +-- Gaussian                              |
               +-- Molpro                                |
               +-- CFOUR                                 |
               +-- MRCC                                  |
                         |                                |
                         +-- ANL task graph               |
                         +-- ANL expression graph         |
                         +-- result.json                  |
                                  |                       |
                                  +-- FinalEnergy resolver+
                                             |
                                      PES / MESS
```

The default user experience becomes:

```json
{
  "theory_preset": "uma-b2plyp-anl",
  "fc_model_path": "/path/to/uma.pt",
  "composite_method": "ANL0-F12"
}
```

while the expert user can independently replace:
- L1;
- L2;
- optimizer/frequency mode;
- any named L3 node;
- any resource profile.

That provides a short default interface without sacrificing theoretical transparency or control.

---

# 65. Source references used to define/review this plan

## KinBot

Repository:
https://github.com/zadorlab/KinBot

Current source files reviewed include:
- `kinbot/qc.py`
- `kinbot/parameters.py`
- `kinbot/calculation.py`
- `kinbot/optimize.py`
- `kinbot/frequencies.py`
- `kinbot/pes.py`
- `kinbot/mess.py`
- `kinbot/ase_modules/calculators/gaussian.py`
- `kinbot/ase_modules/calculators/orca.py`
- `kinbot/tpl/ase_sella_opt_well.tpl.py`
- `kinbot/tpl/ase_sella_hir.tpl.py`
- `pyproject.toml`

## Original ANL

S. J. Klippenstein, L. B. Harding, B. Ruscic,
“Ab Initio Computations and Active Thermochemical Tables Hand in Hand: Heats of Formation of Core Combustion Species”
J. Phys. Chem. A 2017, 121, 6580–6602.
DOI: 10.1021/acs.jpca.7b05945

Open manuscript:
https://www.osti.gov/servlets/purl/1389058

## Perfluoroalkyl ANL-style application

H. Ram, Y. Georgievskii, S. N. Elliott, S. J. Klippenstein,
“Association Kinetics for Perfluorinated n-Alkyl Radicals”
J. Phys. Chem. A 2025, 129, 555–569.
DOI: 10.1021/acs.jpca.4c07388

ACS:
https://pubs.acs.org/doi/10.1021/acs.jpca.4c07388

## ANL0-F12 practical implementation example

“Nascent energy distribution of the Criegee intermediate CH2OO from direct dynamics calculations of primary ozonide dissociation”
J. Chem. Phys.
DOI: 10.1063/1.5028117

## Molpro F12 and gradients

Molpro explicit correlation manual:
https://www.molpro.net/manual/doku.php?id=explicitly_correlated_methods

Molpro energy gradients:
https://www.molpro.net/manual/doku.php?id=energy_gradients

Molpro geometry optimization / F12 CABS-gradient discussion:
https://www.molpro.net/manual/doku.php?id=geometry_optimization_optg

Key implementation fact from the current Molpro manual:

```text
DF-CCSD(T)-F12 analytical gradients exist,
but scaled (T) is not currently available for the analytic gradient.
```

Therefore the implementation should default to:

```text
energy-only F12b:
    SCALE_TRIP=1

analytic derivative F12b:
    SCALE_TRIP=0 unless recipe says otherwise

scaled derivative surface:
    explicit numerical derivative policy
```

---

# 66. Final implementation rule set

If only the most important decisions are retained, retain these:

1. **New functionality is opt-in and flag-guarded.**
2. **Legacy KinBot path stays untouched by default.**
3. **Every profiled optimization uses ASE + Sella/BFGS; programs supply energies/forces.**
4. **Gaussian, Molpro, CFOUR, and MRCC live as KinBot ASE calculator interfaces.**
5. **L1 and L2 are independent profiles.**
6. **Default L1 = UMA/FairChem.**
7. **Lower-tier L2 = Gaussian B3LYP/cc-pVTZ; higher-tier L2 = B2PLYP-D3(BJ)/cc-pVTZ.**
8. **HIR and VPT2 use the selected tier's L2 surface and accepted geometry.**
9. **Frequencies/Hessians are capability-aware, not one-size-fits-all finite differences.**
10. **Generic F12 energy means F12b.**
11. **Energy-only F12b defaults to `SCALE_TRIP=1`.**
12. **Analytic F12 gradients cannot silently use scaled (T).**
13. **Normal open shell = restricted-spin HF/ROHF reference + unrestricted CC.**
14. **Higher-order UHF exceptions are recipe-pinned.**
15. **CFOUR = DBOC + closed-shell CCSDT(Q).**
16. **MRCC = open-shell CCSDT(Q) + all CCSDTQ(P).**
17. **Molpro does most remaining ANL work.**
18. **Gaussian does VPT2.**
19. **Cross-program corrections are allowed and should reuse existing nodes.**
20. **ANL is an expression graph, not a giant L3 template.**
21. **ANL returns complete E0; never add KinBot L2 ZPE again.**
22. **Every expensive result carries exact method/reference/program/geometry provenance.**
23. **A modified recipe must not retain a canonical ANL label.**
24. **Keep existing-file churn under roughly 300–450 LOC if possible.**
25. **Build mixed L1/L2 first; ANL comes only after routing/restart/Hessian behavior is proven.**

---

# 67. Local and external QC validation addendum (2026-09-17)

The development machine has Miniforge but does not have Gaussian, Molpro,
CFOUR, or MRCC. This is a project constraint, not a reason to assume that a
rendered input or parsed output is correct. Create a Python 3.11 Miniforge
environment locally and run all tests that do not require those executables.
The branch's `docs/composite_qc_validation.md` is the detailed, versioned
runbook for the two-site test procedure.

Before implementing each new ASE calculator, research the corresponding
vendor manual and example inputs and outputs. Verify the exact invocation,
working-directory and filename rules, charge and multiplicity, Cartesian
coordinate units, basis and reference keywords, requested derivative or
property, normal-termination marker, final energy or correction label, and
force sign and units. Use official examples plus real, versioned fixtures
from a licensed external site. Sella tests with toy ASE calculators and
mocked execution validate KinBot's plumbing; they do not validate a vendor
program's syntax or scientific result.

The external site must run each program's own example, then the KinBot-
generated input directly and through the ASE wrapper, and compare every
component. The required matrix includes Gaussian B2PLYP-D3(BJ) gradients,
native Hessian and VPT2; Molpro conventional and F12b closed/open-shell
energies and permitted gradients; CFOUR DBOC and closed-shell CCSDT(Q);
MRCC open-shell CCSDT(Q) and both closed/open-shell CCSDTQ(P). Finally run
a mixed L1/L2/HIR restart and compare ANL component energies with published
values. Record program versions, inputs, output hashes, units, tolerances,
references, and exact task provenance. Do not claim the feature complete
until the external results pass; unsupported capabilities must fail clearly.

Early documentation findings affect the implementation: Molpro writes its
own `.out` file, so ASE must capture launcher stdout separately. A native
Molpro `OPTG` geometry optimization also writes `.log` geometry information;
request `SAVEXYZ` to preserve the final `.xyz` and numbered intermediate
geometries. In the KinBot ASE/Sella path, Sella controls the optimization and
Molpro supplies per-step energy and forces, so save the accepted geometry as
an ASE `.xyz` as well. CFOUR reads fixed `ZMAT` and `GENBAS` files, so every
run needs an isolated directory; MRCC reads `MINP` and is invoked with
`dmrcc`; and CFOUR's DBOC is limited to HF/CCSD with RHF/UHF references.
The plan's CFOUR closed-shell CCSDT(Q) node must use and verify the installed
`xncc` capability or be explicitly reassigned by a recipe override. The
detailed sources and checks are in the branch runbook.

---

# 68. L3 geometry barrier and bounded parallel dispatch (2026-09-17)

For each ANL stationary point, submit and accept the recipe's L3 geometry
first. Persist its coordinates and geometry hash before releasing any task
that uses that geometry. Once accepted, make all independent single-point,
harmonic-frequency, and VPT2 tasks ready at the same time. VPT2 should perform
its own prerequisite frequency work when needed so that it need not wait for
a separate harmonic task. Add a dependency between post-geometry tasks only
when the pinned recipe or program actually requires it. Evaluate the final
composite energy only when all required nodes have passed validation.

Bound dispatch with the user-provided maximum number of concurrent exclusive
L3 nodes. Discover node CPU and memory for automatic per-task sizing; optional
task-specific caps remain available after scaling tests. Count submitted and
running jobs against the node limit across all species in a KinBot run.
Reserve resources atomically before submitting and release them on completion
or failure. A task larger than the configured pool must fail preflight with
an actionable error. On restart, recover submitted jobs and reservations
from persistent metadata rather than submitting duplicates. Any change in
accepted geometry invalidates dependent task signatures.

Mock-scheduler tests must prove the geometry barrier, actual overlap of ready
post-geometry tasks, maximum job/core/memory occupancy, restart recovery,
and final result publication only after all required components complete.

---

# 69. Exclusive nodes for every L3 task (2026-09-17)

Every L3 ANL task, including geometry optimization, harmonic frequency,
VPT2, DBOC, and every constituent single point, must have its own exclusive
node allocation. For Slurm, every generated L3 batch script must contain
`#SBATCH --exclusive`; no L3 submission may fall back to a shared node.
Allocate one node per task initially, even when the task uses only part of
the node's cores. Automatic resources request all node memory; explicit
memory caps remain possible for controlled tests. Distinct concurrent L3 tasks use distinct
exclusive allocations. If another queue backend is used, preflight must
verify its equivalent exclusive-node request or reject L3 submission.

Interpret the user-provided L3 concurrency limit as the maximum number of
exclusive nodes KinBot may have submitted or running across the run. Validate
each task's core count and per-rank memory against a single node's capacity. Reserve
one node slot for each submitted/running task, including the geometry task;
release the slot only when that allocation ends. The geometry-first barrier
and post-geometry parallel fan-out from section 68 still apply. Local mode
remains a read/ingest mode for externally run calculations.

Scheduler tests must inspect every generated Slurm L3 script for
`#SBATCH --exclusive`, verify one task per node allocation, and prove that
parallel fan-out never exceeds the configured exclusive-node limit.

---

# 70. Existing KinBot interfaces and MRCC driver review (2026-09-17)

Before writing new QC wrappers, inspect and test the existing KinBot
integration points. `kinbot/tpl/ase_sella_opt_well.tpl.py` already attaches
the custom Gaussian ASE calculator to `Atoms`, runs Sella/BFGS, and exports
XYZ trajectory data. `kinbot/tpl/rotdPy_calc.tpl` hands Molpro method, basis,
memory, and queue settings to external rotdPy; compare a real rotdPy-generated
Molpro input/output pair on the external site because rotdPy's renderer is
not present in this repository. `kinbot/molpro.py` and the legacy Molpro
template can guide filename and variable handling, but their generic F12
call stores `energy(1)` (F12a), whereas the new ANL recipe requires explicit
F12b. Keep legacy behavior stable; do not copy that selection into ANL.

Prefer Molpro→MRCC where the installed software supports the exact task.
The Molpro MRCC manual documents `mrcc,method=ccsdt(q)` and closed-shell
example results. It also says the interface is serial and perturbative
methods are restricted to closed-shell, although its open-shell O2 example
shows some iterative MRCC methods. The MRCC manual describes RHF, UHF,
ROHF, and MCSCF orbitals through its Molpro interface. Because these
statements do not prove that open-shell CCSDT(Q) or CCSDTQ(P) work through
Molpro on the target installation, those nodes must pass direct input/output
tests before using that driver. If unavailable, use direct `dmrcc` or a
separately verified CFOUR→MRCC interface, without changing the method or
reference. Record the chosen driver and program versions in provenance.
This refines section 10's initial direct-driver recommendation.

Official references:
- Molpro MRCC: https://www.molpro.net/manual/doku.php?id=the_mrcc_program_of_m._kallay_mrcc
- MRCC manual: https://www.mrcc.hu/MRCC/manual/pdf/manual.pdf
- Molpro F12 variables: https://www.molpro.net/manual/doku.php?id=quickstart

---

# 71. Prioritized first offsite gate: general dispatch with CH4 (2026-09-17)

The initial priority was a molecule-independent dispatch framework followed
by a CH4-only test. `kinbot/anl/dispatch.py` accepts a declarative molecule
and task graph; CH4 existed only in temporary task data, never as a
conditional in the dispatcher. That runnable
smoke-test generator was retired after the successful offsite run; its general
regression coverage now uses `tests/anl_fixture.py`.
This first gate tests file generation, actual execution, resource bounds,
restart, and geometry sequencing before full ANL recipe arithmetic.

The CH4 fixture first runs Gaussian B2PLYP-D3(BJ)/cc-pVTZ optimization
through the existing ASE calculator and Sella. The accepted geometry feeds
a Molpro CCSD(T)/cc-pVTZ geometry optimization also controlled by Sella
through ASE. A KinBot Molpro ASE calculator runs `RHF; CCSD(T); FORCE,NUMERICAL`
for every requested energy/force evaluation, reads Molpro's printed numerical
gradient in `.out`, then uses `PUT,XYZ` to save Molpro's native geometry for
atom-order mapping. Each step retains its `.inp`, `.out`, `.xyz`,
launcher stdout/stderr, and a detailed Molpro `.log` requested with `-g`;
Sella retains its trajectory,
log, and final `.xyz`. Real licensed output remains required to validate
the documented units, sign, and atom mapping against the installed version.

Only after the accepted Molpro/Sella XYZ passes identity and termination checks,
the dispatcher releases independent Molpro harmonic, F12b/TZ, F12b/QZ,
CCSD(T)/DZ, CFOUR DBOC, and Gaussian VPT2
tasks. The higher-tier VPT2 task uses the accepted B2PLYP-D3(BJ)/cc-pVTZ
L2 geometry and performs `Freq=Anharmonic` at that same level, without
Gaussian `Opt`. It retains an L3-completion dependency to join the post-L3
fan-out while reading coordinates from L2. They each get their own isolated
task directory and one-node Slurm
script with `#SBATCH --exclusive`; at most the user's configured number
of exclusive nodes is in flight. MRCC is deferred from the first offsite
CH4 test at the user's request; the reported site module list does not include
it. The general dispatcher retains MRCC support for a later test once a licensed
executable and the method-specific driver are confirmed.

This test is **not** a complete ANL0-F12 result. It does not assemble a
composite energy or write MESS. The external site's QC outputs must be
inspected and captured as versioned fixtures before scientific parsers and
recipes are enabled. `execution.json` currently records process success,
required artifacts, geometry identity, and artifact hashes; for CFOUR it
does not certify the DBOC value.

The offsite operator should run the exact procedure in
`docs/composite_qc_validation.md`, review discovered site setup/resources,
and return the program versions, generated inputs, stdout/stderr, native
outputs, Slurm logs, and execution manifests or permitted redacted fixtures.
Compare the Molpro input with a rotdPy-generated case on that site.

The first code review found and corrected two concrete bugs in the initial
slice: switching an L1/L2 calculator could retain a previous preset's
method/frequency settings, and translating requested Slurm MB directly to
a Molpro `MEMORY` card would multiply that amount by the process count.
The Molpro 2024 audit supersedes that initial memory choice: this one-node
workflow uses Molpro's default disk implementation and per-process `-m`,
budgeted across the requested MPI ranks with program and node headroom.

---

# 72. Subsequent validation and thermochemistry ladder (2026-09-17)

After the CH4 dispatch workflow passes, use local tests for orchestration
and file handoff, then use the external site for each distinct QC interface:

1. Keep one accepted, provenance-checked electronic-plus-ZPE energy per
   stationary point. A named ANL tier is available only when its required
   components and native output checks all pass. A process exit code alone
   does not establish an accepted energy.
2. Generate CBH-0/1/2/3 reactions directly from the species' SMILES graph,
   with connected capped fragments, elemental balance, and explicit
   electronic states. Read gas-phase 0 K ATcT references from a pinned,
   cached official release. Require an explicit state ID when matching is
   ambiguous or a non-singlet reference is used.
3. Solve 0 K heats of formation from the CBH reaction and accepted 0 K
   energies. Try ANL1-F12, ANL1, ANL0-F12, ANL0, a preselected L3 partial
   composite, then bare L2. Use one exact accepted method for every species
   in a CBH reaction; species elsewhere in the kinetic network may settle
   on different tiers. A tier name whose full equation is not yet pinned
   must not be dispatched as though it were a canonical ANL method.
4. Keep L2 geometries, harmonic frequencies, VPT2 corrections where
   available, and L2 hindered-rotor scans for MESS state counting. Use
   accepted ANL E0 values exactly once for MESS relative zero energies.
   Store each CBH Hf(0), ATcT version, reaction rung, and source identities
   alongside the MESS input. Derive Hf(T) from the MESS thermal increment
   only when a validated MESS thermochemistry reader is available.
5. Extend the graph generator to radicals and other validated states,
   then run CH3 and a small reaction on the external site. Verify TS and
   product restart, L2 rotor propagation, PES MESS assembly, and response
   to a failed constituent calculation. These runs test that the workflow
   can execute and resume, not numerical agreement with a reference heat.

The dispatcher is a reusable staging and execution layer. The remaining
calculator wrappers, scientific output parsers, composite recipes, KinBot
result routing, multi-species resource accounting, and MESS integration
remain separate implementation gates. Later feature work must not treat
process exit status as scientific validation.

ATcT source: https://atct.anl.gov/Thermochemical%20Data/

---

# 73. Pre-smoke audit and reproducible HPC handoff (2026-09-17)

The pre-smoke audit checked the generic dispatcher, its CH4 task graph,
legacy KinBot interfaces, and vendor input syntax. It corrected these
execution risks before the first licensed run:

- A completed geometry or output artifact could be changed after acceptance
  without blocking dependent work. The driver now checks staged input,
  Slurm-script, geometry, and completed-artifact hashes on every advance.
- An `execution.json` could be accepted without a matching task identity or
  complete artifact list. Acceptance now checks the task, schema, geometry
  hash, every declared artifact hash, and the final XYZ hash.
- Site modules, the CFOUR basis file, and required executables had no
  dedicated pre-submission check. `preflight` now sources the same setup file
  as the batch scripts and checks program paths, `CFOUR_GENBAS`, Python
  imports, exclusive Slurm directives, and `sbatch --test-only` for staged
  jobs. Compute-node behavior still
  requires the actual licensed run.
- A failed calculation could only be repeated by starting over. `retry`
  archives a failed attempt and restages that task without rerunning its
  successful siblings; changed chemistry or resources require a new run.
- The CFOUR CH4 input now uses the documented `COORD=CARTESIAN` spelling.
  The Molpro ASE calculator now requests `FORCE,NUMERICAL` after CCSD(T),
  reads the undisplaced total energy, and maps the printed numerical-gradient
  rows back to ASE's coordinates using `PUT,XYZ`. Molpro's F12b `ENERGY(2)`, `SCALE_TRIP=1`, and per-rank
  `-m` memory option were checked against the Molpro 2024 manual; direct MRCC
  `MINP`/`dmrcc` and CFOUR `ZMAT`/`GENBAS`/DBOC were checked against their
  manuals. These checks establish input plausibility, not runtime success.
- The direct MRCC XYZ block now includes the blank line required by its
  documented `geom=xyz` format and explicitly records frozen-core and
  spherical-basis choices. CFOUR must print its DBOC label, and MRCC must
  print normal termination without documented fatal/termination-error text;
  these are interface checks until full output parsers are validated.
- The earlier environment command could install conda-forge Sella 2.1.0 even
  though KinBot requires 2.6.0. The completed HPC procedure pinned ASE
  3.29.0 and installed Sella 2.6.0 from PyPI before editable KinBot.

The historical cluster procedure was retired after the completed first run.
Future site checks should use a recipe-generated workflow on a tested commit
of `composite`. First collect successful native output
and method/version banners. Then compare scientific values and build real
parsers before enabling an ANL label or a MESS handoff. The local regression
suite passes without the four licensed codes; the external result remains the
gate for those program interfaces.

The first licensed Molpro force step must be inspected before releasing all
downstream jobs: compare its CCSD(T) energy and parsed numerical-gradient
forces with the native `.out` table, verify the sign and units by a
small independent finite displacement, confirm atom mapping, and archive the
installed Molpro version. The local fake Molpro test verifies Sella's repeated
calculator calls and artifact handling, but its synthetic output does not
replace this scientific check.

---

# 74. Molpro 2024.1 examples inspected before the ASE/Sella smoke test

The available local reference is
`~/Documents/research/zador_group/peroxy/molpro24_autocbh_conformer_f12_workflow_CURRENT_BACKUP/DZ-F12/work/`.
The CH4 `c000/CH4_c000.inp`, `.out`, and `.log` were inspected directly,
together with its numbered and final `.xyz` files, the CH3 input, and the
H2O2 `c000` and `c001` input/output/log/final-XYZ sets. They are
Molpro's own `OPTG` runs at UCCSD(T)-F12b/VDZ-F12; the CH4 ASE/Sella smoke
test uses conventional CCSD(T)/cc-pVTZ, so they verify Molpro file behavior
and syntax without claiming to validate that different energy label.

- The CH4 input uses an XYZ block with atom count, comment, and coordinates,
  explicit charge/spin, RHF, and a correlated method before
  `OPTG,...,NUMERICAL,...,SAVEXYZ`. The ASE calculator now uses the same XYZ
  block and explicit neutral-singlet settings, followed by RHF and CCSD(T).
- The actual `.out` identifies Version 2024.1, prints the correlated energy
  and `Molpro calculation terminated`, and says its long optimization output
  is in `.log`. The `.log` contains `Numerical gradient for UCCSD(T)-F12B`,
  `Atom dE/dx dE/dy dE/dz` rows, displacement counts, and optimization steps.
  The final and numbered `.xyz` files have the documented count/comment and
  coordinate rows. Both H2O2 conformers have distinct final energies and
  their own `.out`, `.log`, and final `.xyz`. Do not assume the numerical
  gradient table is in `.out` or conflate conformer files.
- The inspected tree contains native `OPTG` inputs rather than standalone
  `FORCE,NUMERICAL` inputs. The first live Sella/Molpro attempt showed that
  `PUT,XYZGRAD` fails after that numerical force command on Molpro 2024,
  even though the complete numerical-gradient table was printed in `.out`.
  The calculator now parses that table, changes Hartree/Bohr gradients into
  ASE eV/Å forces, and uses `PUT,XYZ` for atom-order mapping. A successful
  offsite step remains the acceptance gate for this revised path.
- The Molpro 2024 parallel manual says one-node disk mode is the default and
  recommends `-m` without `-M`/`-G` there. Both the Sella calculator and
  the independent Molpro nodes now divide the requested memory into a
  bounded per-rank stack after reserving at least 200 MW per rank and node
  headroom. MPI ranks run with one OpenMP/MKL thread each.

---

# 75. Portable site discovery before the CH4 run (2026-09-17)

`prepare` now generates `site_setup.sh` from the login-node environment. It
records the loaded module for each QC backend, pins executable directories,
finds Gaussian's adjacent `bsd/g16.profile` and `g16root`, and finds CFOUR's
`basis/GENBAS` beside its installation. It checks both a PATH entry and its
symlink target for adjacent support files. A user-provided `CFOUR_GENBAS`
takes priority. This avoids embedding one cluster's paths or versions in the
generic dispatcher; the operator should still review the generated script
and run preflight before submitting. The actual licensed compute-node run
remains the acceptance gate.

Gaussian scratch selection checks writability of `GAUSS_SCRDIR`, then
`SLURM_TMPDIR`, `SCRATCH`, `TMPDIR`, and finally the task directory. This
handles sites where `/scratch` is exposed in the login environment but not
writable there. The task directory under the user's working area (often
`$HOME`) holds persistent inputs and outputs; scratch is temporary.
The first offsite preflight exposed a Gaussian profile that exited silently
under the batch script's `set -e`; the same setup returned zero when checked
without `errexit`, and the Python imports succeeded. The generated script
now sources the profile in a checked conditional, then reports a nonzero
final status explicitly. This still needs confirmation through the next
offsite preflight and licensed Gaussian job.

For tasks without an explicit partition, preparation reads Slurm's available
partition, CPU, node-memory, and time-limit fields and chooses the shortest
partition that fits the task. An explicit task partition wins. The reported
Blodgett `short-cpu` default has a 30-minute limit, so the CH4 example's
multi-hour tasks should select `day-long-cpu` on that site. Preflight uses
`sbatch --test-only` to catch scheduler restrictions that `sinfo` does not
expose, including account or QoS access. The completed HPC procedure also
inspected the generated setup and selected partition.

---

# 76. First licensed CH4 results and Molpro MPI retry (2026-09-17)

The Gaussian L2 ASE/Sella optimization completed on Blodgett at
`-40.487174956111` Hartree. Its saved final geometry has four C-H distances
between 1.0863 and 1.0865 Å, and the maximum stored atomic force is
0.026075 eV/Å, below the requested 0.03 eV/Å. Sella's final log line shows
0.0349 eV/Å because its logger uses the previous convergence result before
ASE checks the new step. The execution record and saved forces verify this
specific L2 convergence.

The first L3 Molpro task received an exclusive node but failed in four
seconds before creating a Molpro `.out` or `.log`. The launcher `.stderr`
reports `PSM3 can't open nic unit`, OFI endpoint failure, and `PMPI_Init`
abort. Thus the generated Molpro chemistry input and force parser were not
tested by that attempt. For the one-node retry, set
`I_MPI_FABRICS=shm` in the Molpro branch of the run's `site_setup.sh`, based
on Intel MPI's documented shared-memory fabric. The second attempt reached
the numerical gradient, confirming that this setting bypassed the MPI
startup failure on node32. Do not release downstream jobs until Molpro
produces a normal first force step and the input/output/gradient are reviewed.

Blodgett then reported `squeue ... Invalid job id specified` for the already
failed job. The dispatcher previously raised instead of reconciling the
existing failed `execution.json`, so `retry` could not proceed. `_job_active`
now treats only that specific Slurm response as inactive; other scheduler
errors still stop the driver. A regression test exercises reconciliation,
archival, and restaging of the failed task.

---

# 77. Molpro numerical-gradient export correction (2026-09-17)

The second CH4 L3 attempt on Blodgett (`52301032`) ran Molpro 2024 far enough
to finish the CCSD(T)/cc-pVTZ energy and `FORCE,NUMERICAL` evaluation. Its
`.out` contains `SETTING KB_GEOM_ENERGY = -40.43808112 AU` and five rows
under `Numerical gradient for KB_GEOM_ENERGY`, with the expected `dE/dx`,
`dE/dy`, and `dE/dz` columns. It then failed at `PUT,XYZGRAD` with Molpro's
message that the gradient was unavailable for saving. The `.xyz`
was empty, so no forces from that failed attempt were accepted.

The live output reported `Switching to mppx mode, nproc=11` although the
calculator requested `-n 8`. That revision set `MPPX=0` to test whether the
gradient parallel mode caused the process mismatch. The following successful
attempt still launched 12 total processes, disproving that explanation; see
section 79 for the scheduler correction.

The ASE calculator now uses `PUT,XYZ` after `FORCE,NUMERICAL`, reads the
printed numerical-gradient table in `.out`, converts Hartree/Bohr to eV/Å,
negates it to obtain forces, and maps each table row to the ASE atom using
the native XYZ symbols and positions. Missing, duplicate, incomplete,
nonfinite, or unmappable gradient data fail the step. A simulated Molpro
run now exercises this format through multiple Sella steps. The next
offsite retry must confirm Molpro writes the plain XYZ after the numerical
gradient, completes normally, and supplies forces consistent with an
independent finite displacement. The downstream L3 tasks remain held.

---

# 78. Successful Molpro/Sella CH4 geometry and process-count audit (2026-09-17)

The third Blodgett L3 geometry attempt completed. Both Molpro 2024.1 force
evaluations used neutral-singlet conventional CCSD(T)/cc-pVTZ, printed the
`KB_GEOM_ENERGY` numerical-gradient table, wrote a native `.xyz` and `.log`,
and terminated normally. Sella accepted the second geometry. The final
energy in `execution.json` is `-40.43809881` Hartree (the eight-decimal
Molpro `SETTING` value); Molpro's full-precision CCSD(T) line is
`-40.438098814910` Hartree. The final C–H distances are 1.08903–1.08920 Å.
Its maximum atomic force norm is 0.01893 eV/Å, below the requested
0.03 eV/Å. The final force rows match the signs and Hartree/Bohr to eV/Å
conversion of Molpro's second printed gradient. As an additional check,
the energy change between the two Sella geometries is
`-1.7698552e-5` Hartree, while the trapezoidal projection of the two
native gradient tables along that step is `-1.7838087e-5` Hartree.
This confirms the direction and approximate magnitude of the force response
for the observed step; it is not an independent displaced-coordinate test.

The output still reports 11 compute processes plus one helper, with 225 MW
per compute process, despite `#SBATCH --cpus-per-task=8` and the ASE
calculator's `molpro -n 8`. `MPPX=0` prevented the automatic mppx switch
but did not resolve this discrepancy. The site launcher or compute-node
environment needs to be checked before submitting the larger Molpro
harmonic/F12 jobs; otherwise the nominal eight-core and 32 GB budgets
may be exceeded. The task's `--exclusive` directive did reserve a full
96-core node, so the CH4 geometry result itself remains usable. The
downstream tasks can be staged and preflighted, but Molpro jobs should
remain unsubmitted until the process count is explained or corrected.

---

# 79. Exclusive-node resource sizing and Molpro rank correction (2026-09-17)

The Blodgett Molpro 2024 launch script uses Intel Hydra. Under `SLURM_JOB_ID`
it removes its own `-np` launcher option, so `molpro -n 8` alone cannot
control the MPI process count. [Intel's Slurm integration guide](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/2021-10/job-schedulers-support.html)
lists `SLURM_NTASKS` and `SLURM_CPUS_PER_TASK` as Hydra inputs. The old
`#SBATCH --cpus-per-task=8` described
one eight-core Slurm task and left Hydra to infer the count; the CH4 output
reported 12 total processes (11 compute plus one helper). New Slurm scripts
request `--ntasks=<Molpro ranks> --cpus-per-task=1`; Gaussian and CFOUR
retain one task with the requested shared-memory CPU count. Every L3 script
still requests `--exclusive`. Molpro's default mppx numerical-gradient mode
is restored, because disabling it did not correct the launcher count. The
runner now rejects Molpro output that reports more total MPI processes than
the declared rank allocation. A small live DZ CCSD(T) job is the next gate;
the larger Molpro jobs remain held until its reported count agrees.

The operator now needs to provide only `limits.max_nodes` for an automatic
Slurm run. `prepare` reads the selected partition's `sinfo` node CPU count
and memory, uses its smallest eligible node as the safe sizing basis, and
requests all node memory with [`--mem=0`](https://slurm.schedmd.com/sbatch.html).
For Molpro, automatic sizing budgets
only 85% of that memory and reserves the [manual's 200 MW per-process program
overhead](https://www.molpro.net/manual/doku.php?id=running_molpro_on_parallel_computers)
plus a configurable minimum stack of 1024 MW per rank. It fails preparation
if even one rank cannot meet that minimum. The memory-affordable rank count
is reduced to an efficient candidate in 1, 2, 4, 8, 12, or 16, with a default
16-rank performance ceiling. `resources.max_cores` can lower or explicitly
raise that ceiling for a particular method after timing measurements; no
single count is assumed optimal for every CCSD(T), F12, or harmonic job.
Gaussian and CFOUR use an analogous configurable memory-per-core floor.
Explicit cores/memory remain available for reproducible smoke fixtures; the
existing `ch4_run_auto3` retains its prepared 4/8-core, 16/32/48-GB values.

The next Blodgett procedure is: pull `composite`, reconcile the completed L3
geometry with `advance(..., submit=False)`, preflight the newly staged jobs,
then run `drive ch4_run_auto3 --once --only ccsdt_dz`. Inspect its Slurm script
for `--ntasks=4 --cpus-per-task=1`, and its native output for at most four
total processes, the intended CCSD(T)/cc-pVDZ method, and normal
termination. Only then release the independent harmonic, F12, CFOUR DBOC,
and Gaussian VPT2 jobs up to the user's node limit. The CH4-only fixture and
its temporary tests may be removed after this acceptance pass; the general
resource resolver, Slurm layout, and process-count check remain.

---

# 80. CH4 Molpro rank-count probe passed (2026-09-17)

The live Blodgett `ch4_run_auto3/tasks/ccsdt_dz` job exited zero and its
`execution.json` passed dispatch artifact checks. Molpro 2024.1 reports four
total MPI processes (three compute plus one helper), matching the requested
`--ntasks=4`, and 225 MW per compute process. The native input uses the
accepted L3 geometry, neutral-singlet RHF, `basis=cc-pVDZ`, and `ccsd(t)`.
Molpro printed `!CCSD(T) total energy = -40.387076267138` Hartree,
`SETTING KB_CCSDT = -40.38707627 AU`, and its normal-termination marker.
Thus the Slurm/Hydra rank correction is verified for this four-rank job.
Eight-rank F12 and harmonic jobs still need their own live output checks.

The historical smoke-fixture name `ccsdt_dz` and its `KB_CCSDT` variable
falsely suggest full CCSDT; this output is CCSD(T) and cannot be used as an
ANL CCSDT correction. New generated CH4 fixtures use `molpro_dz_sp` and
`KB_DZ_ENERGY`, while the already prepared immutable run retains its old
identifiers. The next step is to reconcile the completed probe, preflight,
and release the remaining independent jobs up to the configured exclusive
node count. Their successful process exits are only dispatch gates; each
native output still needs method, energy, frequency, DBOC, and VPT2 review
before ANL assembly or MESS handoff.

---

# 81. Refresh completed jobs when reporting status (2026-09-17)

After the first post-geometry fan-out, Blodgett's queue was empty while the
saved state still showed CFOUR DBOC and Gaussian VPT2 as `submitted`. This is
expected from the original `status` implementation, which only read
`state.json`; `advance(..., submit=False)` reconciles finished jobs and
updates their accepted status without new submissions. The default `status`
command now performs that reconciliation under `drive.lock` when no driver
owns the run. If a polling driver holds the lock, it reports the driver's
last atomic snapshot. `status --cached` explicitly requests the saved
snapshot without a Slurm query. A regression test verifies that refreshing
an exited job marks it complete and stages its child without calling
`sbatch`. The live Gaussian VPT2 process has since passed the dispatch gate;
its scientific output remains unverified.

---

# 82. Resolve CFOUR's missing compiler runtime (2026-09-17)

The CH4 fan-out completed seven of eight dispatch tasks. CFOUR DBOC failed
before reading `ZMAT`: `xcfour` exited 127 because the compute node's dynamic
loader could not find `libgfortran.so.4`. Its `cfour.out` was empty. This is a
runtime dependency of the site's CFOUR executable, not evidence of a bad DBOC
input or a completed DBOC result. Both login and compute-node `ldd` showed the
same missing library. The site's `/opt/anaconda3/lib/libgfortran.so.4` exists
on the node; a short `srun` probe with that file alone in `LD_PRELOAD` resolved
the CFOUR executable and its `libquadmath`/`libgcc_s` dependencies. A native
CFOUR calculation still has to verify that this runtime works through DBOC.

The dispatcher now checks CFOUR ELF dependencies after site setup and again
on the compute node just before launch. If the module already resolves them,
it does nothing. For a missing `libgfortran.so.N`, it looks for an exact
SONAME in the CFOUR/Conda installation and shallow common software roots,
checks the candidate with `ldd`, and preloads only that file into the CFOUR
child process. It does not add a complete Conda library directory to
`LD_LIBRARY_PATH`, which could replace unrelated libraries for Python, Git,
or Slurm. `CFOUR_LIBGFORTRAN` is an explicit fallback for an unusual site.
The chosen path is recorded in `execution.json`; unresolved dependencies fail
early with a specific message. This also works for the already prepared CH4
run, whose task input and site setup need no edits. After updating the branch,
preflight, archive/retry only `cfour_dboc`, and submit that task. Keep the
failed attempt and the eventual native output for provenance. A successful
process/marker gate still does not certify the numerical DBOC; parse and
review its value before any ANL arithmetic.

---

# 83. Respect CFOUR 2.1's fixed-width keyword input (2026-09-17)

The next Blodgett CFOUR attempt passed the loader and launched `xjoda`, then
failed before SCF. Its echoed `ZMAT` stopped at `...,CHARGE=0,MU` and reported
`Must supply value for keyword string MU`. The generated `*CFOUR(...)` line
was 123 characters; CFOUR 2.1's reader used only its first 80 columns. The
input was therefore truncated before `MULTIPLICITY`, independently of the
earlier shared-library failure. No DBOC result has been accepted.

[CFOUR's keyword-section manual](https://cfour.uni-mainz.de/cfour/index.php?n=Main.CfourKeywordSection)
allows continuation across lines, requires no trailing comma at a continued
line end, and closes the parenthesis only on the last line. The CH4 template
now follows this syntax. More generally, every staged CFOUR `*CFOUR(...)`
section is normalized to one keyword per line, preserving nested method
parentheses such as `CCSD(T)` and rejecting any individual keyword longer
than a conservative 72 columns. This stage-time formatting also repairs a
retry of the already prepared immutable workflow, since `retry` regenerates
`ZMAT` from its original template using the updated code. The previous
attempt's input and output remain archived. The next live retry must confirm
the printed keyword table, SCF, DBOC value, and normal CFOUR termination.

---

# 84. First complete CH4 dispatch and CFOUR result gate (2026-09-17)

The next retry completed: all eight CH4 tasks now report `complete`. CFOUR
2.1 ran on a compute node using an automatically selected
`/opt/anaconda2/lib/libgfortran.so.4.0.0`, exited zero, and archived hashes
of `ZMAT`, `GENBAS`, native output, and execution record. Its echoed input
contains `BASIS=cc-pVTZ` and `MULTIPLICITY=1`; the keyword table reports PVTZ
and 86 basis functions. The final output reports the HF DBOC as
`0.0025887093` Hartree, `568.156016 cm-1`, and `6.797 kJ/mol`. These three
values agree within print precision. The same output also reports an MP1 DBOC
of `0.0026718675` Hartree; this is a distinct level and must not be selected
by matching the last occurrence of the repeated DBOC label. CFOUR's final
electronic energy `-40.210746196622267` Hartree includes DBOC and is not the
DBOC correction itself.

The first method-aware native parser now selects the named HF summary,
cross-checks its reported units, retains the MP1 summary as separately
labeled provenance, and verifies CFOUR's completion and final energy.
Future CFOUR DBOC tasks can declare `result_parser` so an invalid or missing
value fails the task; the CH4 fixture declares HF. The already completed
immutable run predates this declaration, so a read-only parser command will
inspect its existing `cfour.out` without rewriting the accepted
`execution.json`.

---

# 85. Full CH4 native-output review and method-aware result gates (2026-09-17)

The completed CH4 text archive was inspected locally. It contains native
inputs, outputs, scheduler records, workflow state, and execution hashes for
all eight tasks. The original full run remains on the HPC; the review archive
omits large CFOUR binary scratch, `GENBAS`, and ASE trajectories. Every
included file named in an execution artifact map matched its SHA-256 hash.
The audit findings and exact component values are recorded in
`docs/composite_qc_validation.md`.

The Molpro F12 outputs printed both F12a and F12b totals, so a generic
`energy` or last-energy parser would be scientifically ambiguous. The native
F12b values are `-40.454906199189` Hartree at cc-pVTZ-F12 and
`-40.456608306474` Hartree at cc-pVQZ-F12. The reported `ENERGY(2)` variable
matches F12b but is rounded to eight decimal places in the assignment line;
the parser now reads the exact named `!CCSD(T)-F12b total energy` and checks
the final method summary. The DZ output is conventional CCSD(T), not CCSDT,
despite the historical `ccsdt_dz` task ID. Molpro's numerical harmonic task
reports nine positive CH4 modes and ZPE `0.04479801` Hartree. The parser
checks the half-sum of vibrational wavenumbers against the reported ZPE and
separates zero rotation/translation modes.

Gaussian's B3LYP/cc-pVTZ Opt/Freq=Anharmonic task completed at a stationary
point on its own B3LYP surface and reports harmonic ZPE `9783.68667 cm-1`,
total anharmonic ZPE `9646.36482 cm-1`, and correction `-137.32185 cm-1`.
The output also contains two warnings about unreliable cubic force constants
and one rotor/framework warning. Near-degenerate Darling-Dennison resonances
were active. The VPT2 result is parsed and preserved with
`review_required=true`, **not accepted for final thermochemistry yet**. A
future calculation at the selected tier must check the ZPE components and
warning behavior before the correction is used in an ANL expression. The
user does not request another CH4 smoke run now. The run's `complete` status
records dispatch and output
completion, not that convergence assessment.

External tasks can now declare method-aware `result_parser` requests. The
dispatcher validates each parser against its backend, expected basis, and
input method before submission; after normal program completion it extracts
the named result and fails malformed or inconsistent output. The CH4 fixture
declares parsers for both F12b points, conventional DZ, harmonic frequencies,
Gaussian VPT2, and CFOUR HF DBOC. These are reusable parser kinds; no
CH4-specific dispatch branch was added. The previously completed prepared
workflow is immutable and lacks the new declarations, so its native outputs
are parsed read-only with the CLI rather than editing accepted records.

Remaining before an ANL electronic energy or MESS handoff: resolve the VPT2
quality flag; establish full ANL0-F12 and later ladder expressions with their
missing higher-order corrections and open-shell/MRCC checks; then validate
the arithmetic and MESS mapping on small molecules and a small reaction.
The CH4 test alone does not authorize an ANL heat of formation.

---

# 86. Tie VPT2 and hindered rotors to the selected L2 tier (2026-09-18)

The user clarified the intended two-tier surface policy. Lower ANL-tier
work uses B3LYP/cc-pVTZ for L2 geometry, VPT2, and hindered-rotor scans.
Higher ANL-tier work uses B2PLYP-D3(BJ)/cc-pVTZ for those three calculations.
The corresponding theory presets are `uma-b3lyp-anl-low` and
`uma-b2plyp-anl`; absent an explicit preset, ANL0/ANL0-F12 select the lower
profile and ANL1/ANL1-QZF select the higher profile. Explicit L2 profile
overrides remain possible and must propagate consistently to VPT2 and HIR.
An ANL1 energy with the higher-tier B2PLYP correction is a profiled variant
of the published B3LYP-corrected equation, so its final provenance must say
so rather than claiming an exact published ANL1 reproduction.

A Gaussian `Freq=Anharmonic` job does **not** optimize unless `Opt` is present.
The completed first CH4 smoke fixture did include `Opt=(Tight,CalcFC)` at
B3LYP/cc-pVTZ, and its log explicitly shows a stationary-point search. That
old result is a historical B3LYP interface test, not the intended higher-tier
B2PLYP anharmonic correction. The revised test-only graph omits Gaussian
`Opt`; it reads the accepted L2 geometry, requests B2PLYP/cc-pVTZ with
`EmpiricalDispersion=GD3BJ`, and keeps a dependency on L3 geometry completion
so all post-L3 tasks can still be dispatched together. The L2 Sella
threshold is tightened to `0.0005 eV/Å` and both L2 and VPT2 request an
UltraFine integration grid. The general spec validator rejects a
frequency-only VPT2 node if its Gaussian geometry-source method, basis, or
dispersion does not match the declared VPT2 level. The parser records whether
the native VPT2 output optimized in that job.

No replacement CH4 smoke calculation is requested now. Keep the completed
run immutable. A future licensed-site validation of the profiled higher-tier
VPT2 result should use a generated recipe workflow and review its warnings.

---

# 87. Post-smoke cleanup and composite assembly foundation (2026-09-18)

The temporary runnable CH4 generator and step-by-step HPC runbook are removed.
The native-output findings above remain as an audit record. A small task graph
now lives only in the tests; it is not a production preset or a submission
script. General dispatcher, scheduler, Molpro/Sella, and parser regression
tests remain.

`kinbot/anl/model.py` evaluates complete electronic and zero-point expressions
from method-labeled components. `kinbot/anl/recipes.py` declares ANL0,
ANL0-F12, and ANL1 term lists; the chosen scaled-triples F12 implementation
and any B2PLYP-D3(BJ) VPT2 substitution receive distinct profiled labels.
Missing components, wrong method/basis/backend,
stale L2/L3 geometry, mismatched electronic state, unreviewed native warnings,
or missing source hashes block the result. Spin-orbit is an explicit term,
including when a state-specific provider concludes it is zero. Synthetic
tests cover F12b selection, cross-program higher-order subtraction, both
energy levels, and these refusal paths.

`kinbot/anl/workflow.py` can convert a *newly prepared* completed task's
method-aware parsed output into a component only after checking staged files,
artifact hashes, and a fresh parse of the native output. It currently covers
individual Molpro energies/harmonic ZPE, Gaussian VPT2 correction, and CFOUR
DBOC. The first completed CH4 run predates some parser declarations and is
kept as historical validation rather than being silently upgraded.

**Next implementation gates:** build core-valence, scalar-relativistic, higher-order
CCSDT(Q)/CCSDTQ(P), and
state-specific spin-orbit providers with exact references and program versions.
Then connect the recipe graph and restart engine to KinBot's L3 boundary,
route the assembled 0 K result through PES and MESS without an extra L2 ZPE,
and validate on small closed/open-shell species and a small reaction. No
complete ANL energy, heat of formation, or MESS handoff is claimed yet.

---

# 88. Original recipe equation and two-point coefficients verified (2026-09-18)

The [original manuscript, section 2.1, equations 1–3](https://www.osti.gov/servlets/purl/1389058)
was checked on its rendered pages 6–8. ANL0-F12 is *exactly* ANL0 with the
conventional `a'QZ/a'5Z` CCSD(T) CBS reference removed and the
`cc-pVTZ-F12/cc-pVQZ-F12` CCSD(T)-F12b CBS reference added. It retains the
ANL0 CCSD(T)/TZ geometry and harmonic-ZPE framework. ANL1 instead uses
CCSD(T)/QZ geometry, a TZ/QZ harmonic-ZPE CBS term, the `a'5Z/a'6Z`
electronic CBS pair, and both TZ CCSDT(Q) and DZ CCSDTQ(P) increments.

The common core-valence correction is the difference between *all-electron*
and *frozen-core* CCSD(T) energies, each extrapolated from `cc-pcVTZ` and
`cc-pcVQZ`. The DKH term is the CCSD(T) energy difference with and without
Douglas–Kroll one-electron integrals using `aug-cc-pcVTZ-DK`; DBOC is
HF/cc-pVTZ. The paper obtains state-specific spin-orbit corrections from
experiment. The published anharmonic correction uses B3LYP/cc-pVTZ. The
user-selected B2PLYP-D3(BJ) higher-tier VPT2 term remains a *profiled*
variant; the original paper reports numerical instability for several DFT
functionals including B2PLYP-D3, so a normal Gaussian exit alone is not
sufficient to accept it.

For adjacent cardinal numbers `n-1,n`, the paper uses the total-energy rule

```text
E_CBS = E_n + alpha_n (E_n - E_(n-1))
alpha_n = (n-1)^3.7 / (n^3.7 - (n-1)^3.7)
```

| Basis pair | Upper `n` | Derived `alpha_n` | Paper rounded value |
| --- | ---: | ---: | ---: |
| TZ/QZ, TF/QF, cTZ/cQZ | 4 | 0.52654647 | 0.53 |
| a'QZ/a'5Z | 5 | 0.77922802 | 0.78 |
| a'5Z/a'6Z | 6 | 1.03817643 | 1.04 |

The paper also notes approximate RMSD-minimizing values near 0.5, 0.75,
and 1.1, respectively. Those are observations, not the adopted coefficients.
The original paper says its `l^-3.7` rule is used throughout; later
ANL-family variants may pin different coefficients and must receive distinct
recipe provenance. `kinbot/anl/extrapolation.py` implements this arithmetic.
The 2017 article identifies the F12b approximation but does not specify
whether perturbative triples were scaled. [Molpro's manual](https://www.molpro.net/manual/doku.php?id=explicitly_correlated_methods)
defines `SCALE_TRIP=1` as an explicit change. The currently generated F12
inputs use it, so their recipe label remains a profiled ANL0-F12 variant
until original computational inputs or component values resolve that setting.
Using the accepted native CH4 F12b/TZ and F12b/QZ totals only as a numeric
regression gives `-40.457504545051066` Hartree for their CBS reference.
That number is one component, not a complete ANL0-F12 energy.

---

# 89. Verified Molpro pair to CBS component (2026-09-18)

`kinbot.anl.workflow.cbs_task_component` accepts two completed electronic
Molpro task IDs and a recipe requirement. It reparses both native outputs after
checking dispatcher stage and execution hashes, then requires the declared
basis order, method, electronic state, geometry hash, quantity, and native
calculation settings to agree. It rejects unfinished tasks, warnings requiring
review, changed output files, and missing extrapolation parameters. The
derived component records the exact exponent, cardinal number, input paths,
and a deterministic digest of both native output hashes. It can provide the
F12 T/Q reference and, once those native pairs are generated, conventional
CCSD(T) electronic CBS references. A synthetic completed
F12 T/Q task pair checks the accepted CH4 numeric result and rejection paths;
this does not represent a new licensed QC run.

**Next:** produce and parse the conventional `a'QZ/a'5Z` and `a'5Z/a'6Z`
Molpro basis pairs, then add four validated all-electron/frozen-core tasks for the
core-valence difference. The scalar-relativistic, higher-order CC, and
spin-orbit providers, full recipe assembly, KinBot PES/MESS handoff, and
small-species/reaction validation remain open.

---

# 90. Geometry-aware harmonic and anharmonic CBS components (2026-09-18)

The CBS task-pair provider now also accepts Molpro harmonic ZPE and Gaussian
VPT2 anharmonic-correction pairs. Each TZ and QZ zero-point task must use a
*separately completed and verified* ASE/Sella geometry optimization at its
own method, basis, and, for Gaussian, dispersion setting. The provider checks
both native output hashes and both geometry-task artifact hashes, rejects
crossed levels or missing quality review, and records each geometry source in
the derived component digest. Its result geometry is the QZ geometry, which
the recipe evaluator must match to the accepted highest-level L2 or L3
geometry. Electronic CBS reference single points must share one completed
Molpro geometry optimization at the recipe's declared method and geometry
basis: CCSD(T)/TZ for ANL0 and ANL0-F12 or CCSD(T)/QZ for ANL1. Other
electronic single point providers must use the highest available geometry
when implemented. This keeps the TZ and QZ ZPE calculations
on their own level-matched structures.

ANL1's published TZ/QZ harmonic-ZPE CBS requirement now declares that
geometry mode. `recipe(..., vpt2_cbs=True)` opts into a separately labeled
TZ/QZ extrapolation of the B3LYP or B2PLYP-D3(BJ) anharmonic correction;
the original 2017 recipes use a TZ-only anharmonic correction. The two-point
power and basis pair remain explicit in the profile so alternate extrapolation
choices can be identified. Synthetic completed task graphs cover both ZPE
kinds and source-geometry tampering. They are parser and provenance tests,
not a claim that the new TZ/QZ frequency calculations have run on Blodgett.

For either type of ZPE term, let `X_TZ` and `X_QZ` be values parsed from
separately optimized TZ and QZ calculations. The extrapolated value is
`X_QZ + alpha_4 (X_QZ - X_TZ)`, where `alpha_4 = 0.5265464668` for the
configured `l^-3.7` profile. For harmonic ZPE, `X` is the native harmonic
ZPE. For VPT2, `X` is the within-run anharmonic minus harmonic ZPE, so the
correction is extrapolated once and then added once to the independent
CCSD(T) harmonic term.

---

# 91. Local CBH, ATcT, ladder, and MESS handoff fixtures (2026-09-18)

The CBH implementation builds capped graph neighborhoods directly from
KinBot's SMILES connectivity, using the atom, bond, and larger-neighborhood
construction of the [original CBH paper](https://pubs.acs.org/doi/10.1021/ct200279q).
It does not import reaction-generation code from a CBH library. It currently
generates neutral, closed-shell, nonaromatic C/N/O/F/Cl CBH-0 through CBH-3
reactions, deduplicates fragments, and checks element balance. Unsupported
radicals, ions, aromatic states, and unusual valences fail explicitly.
The 0 K solver can accept an externally specified radical reaction with
explicit spin states and an explicit ATcT ID for a non-singlet reference;
systematic radical fragment generation is a separate remaining step.

The ATcT reader fetches and caches a pinned official release, verifies its
release header, retains the source SHA-256, and requires a unique gas-phase
0 K state match or an explicit reference ID. It uses the versioned
[ATcT tables](https://atct.anl.gov/Thermochemical%20Data/), not a moving
value embedded in the recipe. A completed calculation keeps its original
reference snapshot even if a newer ATcT release appears.

Local table fixtures check four of the supplied schemes: ethane CBH-0,
CH3CF3 CBH-1, C2F6 CBH-2, and C3F8 CBH-3. The CH3 radical CBH-0 scheme
checks explicit doublet state handling in the 0 K solver. The image table
prints component totals with the opposite sign to the displayed forward
reaction arrows. Internally, reactants are negative and products positive:

```text
C2H6 + H2 -> 2 CH4
table component total: +15.345 kcal/mol
forward reaction delta H(0): -15.345 kcal/mol
Hf,0(C2H6) = 2 Hf,0(CH4) - delta H(0) = -16.461 kcal/mol

CH3CF3 + 3 CH4 -> C2H6 + 3 CH3F
table reverse-signed reaction total: -44.938 kcal/mol
forward reaction delta H(0): +44.938 kcal/mol
Hf,0(CH3CF3) = Hf,0(C2H6) + 3 Hf,0(CH3F)
                 - 3 Hf,0(CH4) - delta H(0) = -176.728 kcal/mol
```

The tests use the table's rounded ATcT 1.202 kcal/mol numbers as synthetic
fixtures, not a downloaded ATcT 1.202 release or native electronic
calculations. They successively remove one reference species' accepted
energy and verify fallback through ANL1-F12, ANL1, ANL0-F12, ANL0,
a preselected L3 profile, and L2. The selected method is common to all
species of one CBH reaction. A separate local MESS network test confirms that
different stable species may use different accepted tiers after each species
has been converted to the common physical formation-enthalpy reference. A
stationary TS barrier still requires one consistent accepted method for the
reactant and TS.

KinBot's direct MESS writer writes a `me/formation_0k.json` sidecar when a CBH
Hf(0) is attached. The sidecar records the rung, 0 K value, method, ATcT
release/hash, reference IDs, energy sources, network reference, and stationary
barrier provenance. The kinetic energy axis uses differences of Hf(0), rather
than differences of unrelated absolute electronic energies. A fake local MESS
handoff checks the sidecar, energy axis, and L2 hindered-rotor writer path.

**Remaining implementation gates:** pin and implement every advertised ANL
recipe, extend radical/charged CBH graph rules, propagate uncertainty, compare
the local NASA7 fitter with PAC99 on real MESSPF output, and run a small
TS/reaction restart test on the external site. The local fixtures establish
reaction algebra and file handoff only.

---

# 92. Ram MESS kinetics and partition-function procedure (2026-09-18)

The supplied `jp4c07388_si_001` files from the
[Ram et al. association study](https://pubs.acs.org/doi/10.1021/acs.jpca.4c07388)
were compared with KinBot's MESS writer.
They separate the thermochemical roles that the earlier handoff had conflated:

1. A standalone MESSPF input gives each species `ZeroEnergy = 0` and contains
   its geometry, corrected frequencies, electronic degeneracy, and explicit
   hindered rotors. It produces partition functions, entropy, and heat
   capacity data.
2. A kinetic MESS input repeats the same state-counting model, but its well
   `ZeroEnergy` and fragment-channel `GroundEnergy` values are differences of
   CBH/ANL Hf(0) values on one arbitrary reference axis.
3. A stationary TS is placed at `Hf,0(reactant) + barrier_0K`, where the
   barrier is an accepted same-method ANL E0 difference. A barrierless VRC
   transition state is placed at its fragment-channel threshold and does not
   acquire a fictitious stationary-saddle energy.

For example, the SI C2F6 input uses the two-CF3 channel as zero:

```text
Hf,0(C2F6) - 2 Hf,0(CF3)
= -318.333 - 2(-111.164)
= -96.005 kcal/mol
```

KinBot now activates this route whenever a composite method is requested or a
CBH formation result is attached. It fails if a stable network species lacks
a validated Hf(0), or if a stationary barrier cannot be derived from a
consistent accepted reactant/TS pair. Per-species ladder fallback is valid
after Hf(0) conversion because every stable species is then on the same
formation scale. The old direct difference of absolute electronic energies is
bypassed in this mode. The legacy MESS route is unchanged for non-ANL runs.

For every stable network species KinBot also writes
`me/partition_functions/<chemid>.inp`. These inputs use zero energy, include an
explicit 298.15 K point, and call the same geometry, frequency, electronic
degeneracy, and L2 rotor serializers as the kinetic input. If an accepted
harmonic-plus-VPT2 frequency array is attached it is used in both places;
otherwise both use the selected L2 frequency array.
`run_messpf.sh` runs the configurable `messpf_command` (or the
`KINBOT_MESSPF_COMMAND` environment override) on every input. The reader uses
the official executable's `<chemid>.dat` output name.

The SI frequency comments such as `VTZ-F + B2D3/VQZ ANH` mean a high-level
harmonic frequency plus a mode-specific VPT2 shift, not the unmodified DFT
fundamental. KinBot therefore parses every Gaussian `Fundamental Bands` row
and constructs

```text
nu_i(corrected) = nu_i(target harmonic)
                + nu_j(VPT2 anharmonic) - nu_j(VPT2 harmonic).
```

The mode assignment minimizes the total harmonic-frequency mismatch and has a
hard mismatch tolerance; a reviewed explicit one-based map can replace the
automatic assignment. The target list is `reduced_freqs`, after KinBot has
projected the internal rotations represented as MESS hindered rotors. Thus the
Gaussian torsional oscillator is not counted again. Gaussian quality warnings,
missing modes, nonpositive fundamentals, bad mode assignments, and nonpositive
corrected frequencies stop this handoff. The verified native output hash and
mode map are recorded in the partition-function manifest. When the target
harmonics and VPT2 force field are the same L2 method and geometry, this formula
reduces to the Gaussian VPT2 fundamentals requested for the profiled L2 model.

The [official MESSPF source](https://github.com/Auto-Mech/MESS/blob/master/src/partition_function.cc)
reports `Z_1 = d ln(Q)/dT`, entropy, and constant-pressure heat capacity. It
also appends 298.2 K automatically, so KinBot inserts and selects a separate
298.15 K point. With the MESSPF species
ground set to zero, KinBot evaluates

```text
H_species(T) - H_species(0) = R T^2 Z_1(T) + R T
```

where the final `RT` is the ideal-gas `pV` term. Formation enthalpy at 298.15 K
then follows from

```text
Hf,298(species) = Hf,0(species)
                 + [H298-H0](species)
                 - sum_e n_e [H298-H0](element e, standard state per atom).
```

The parser requires a real 298.15 K MESSPF row. The elemental increments are
pinned to NIST-JANAF reference-state tables: H2 8.467, graphite 1.051, N2
8.670, O2 8.683, F2 8.825, and Cl2 9.181 kJ/mol from 0 to 298.15 K; molecular
values are divided by two per atom.

KinBot now fits two NASA7 heat-capacity polynomials to the MESSPF grid. The
low-range integration constants reproduce Hf(298.15) and the MESSPF entropy at
298.15 K. The high-range integration constants reproduce the low-range
enthalpy and entropy at 1000 K. The JSON record contains both seven-coefficient
sets, temperature limits, anchors, and Cp root-mean-square errors. This is an
internal PAC99-compatible fit that still needs comparison against PAC99 before
it is treated as a replacement for that program. Hf(0) remains the kinetic
MESS energy anchor; Hf(298) anchors the thermochemical fit.

---

# 93. Full workflow audit and first profiled external-site test (2026-09-20)

A full pass found that theory profiles could be parsed but could not execute:
`QuantumChemistry` deliberately raised `NotImplementedError`. The opt-in
`ProfiledQuantumChemistry` router now creates independent legacy backends,
uses FairChem UMA for L1 reaction/conformer work, sends accepted stationary
points and hindered rotors to the selected Gaussian L2 surface, and persists
the backend plus scheduler ID for every job. Legacy inputs still construct the
original `QuantumChemistry` class. Scheduler scripts now use Slurm
`--partition` and the exact Python interpreter that submitted the job, which
is required for FairChem installed in a repository-local environment.

The first external-site example is
`examples/anl/ethane_profiled_hpc`. It performs a restricted ethane C-C homolysis
search with conformer processing at UMA L1, B2PLYP-D3(BJ)/cc-pVTZ L2
refinement and hindered rotors, then exports the accepted parent L2 record to
a general non-MRCC interface graph. That graph runs the Molpro
CCSD(T)/cc-pVTZ ASE/Sella geometry first and releases Molpro harmonic/F12/DZ,
CFOUR DBOC, and frequency-only Gaussian VPT2 nodes concurrently under the
user's node cap. All dispatcher nodes are exclusive and Molpro ranks are
bounded by both per-rank memory and method scaling caps.

This is not yet the requested production ANL1-F12/CBH/MESS test. Its required
terminal label is `interface_complete_recipe_incomplete`. The branch has no
citable canonical `ANL1-F12` equation; MRCC and its CCSDTQ(P)/DZ increment are
disabled; core-valence, scalar-relativistic, state-specific spin-orbit, and
closed-shell higher-order providers are not all connected; and L3 graph
creation is not yet automatic for every accepted well/product/TS. The next
implementation step is to complete those providers and attach the dispatcher
at the stationary-point acceptance boundary, then run CBH/ATcT assembly and
the existing ANL MESS handoff for every network state. A stationary-TS/IRC
restart test is required after the ethane interface run.

The exact install, pull, run, restart, and acceptance procedure is in
`docs/ANL_full_HPC_validation.md`.
