# KinBot Profiled ASE/Sella + ANL Composite Methods
## Final revised implementation plan

**Repository reviewed:** `zadorlab/KinBot`, current `master` as of 2026-09-16
**Observed upstream state during review:** KinBot 2.4.1; current master head reported by GitHub activity as `924faba`; no open pull requests at the time of review.
**Primary objective:** implement a compact, opt-in, backward-compatible hierarchy

```text
UMA / FairChem L1
        ↓
B2PLYP-D3(BJ)/cc-pVTZ L2
        ↓
ANL0 / ANL0-F12 / ANL1 / future pinned ANL-family recipe
```

while retaining full user control over L1, L2, and every constituent L3/composite calculation.

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
+ CV_CBS(cTZ,cQZ)
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
+ CV_CBS(cTZ,cQZ)
+ REL_DKH
+ DBOC_HF/TZ
+ SO
```

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
7. **Default L2 = Gaussian B2PLYP-D3(BJ)/cc-pVTZ.**
8. **Default HIR = L2.**
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

Bound dispatch with user-provided global L3 limits for simultaneous jobs,
cores, and memory, in addition to per-task resource profiles. Count submitted
and running jobs against those limits across all species in a KinBot run.
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
the node's cores or memory. Distinct concurrent L3 tasks must use distinct
exclusive allocations. If another queue backend is used, preflight must
verify its equivalent exclusive-node request or reject L3 submission.

Interpret the user-provided L3 concurrency limit as the maximum number of
exclusive nodes KinBot may have submitted or running across the run. Validate
each task's core and memory request against a single node's capacity. Reserve
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

The user's current priority is a molecule-independent dispatch framework,
followed by a CH4-only test of it. `kinbot/anl/dispatch.py` accepts a
declarative molecule and task graph; CH4 is data in
`examples/anl/ch4_dispatch.py`, never a conditional in the dispatcher.
This first gate tests file generation, actual execution, resource bounds,
restart, and geometry sequencing before full ANL recipe arithmetic.

The CH4 fixture first runs Gaussian B2PLYP-D3(BJ)/cc-pVTZ optimization
through the existing ASE calculator and Sella. The accepted geometry feeds
a Molpro CCSD(T)/cc-pVTZ geometry optimization also controlled by Sella
through ASE. A KinBot Molpro ASE calculator runs `RHF; CCSD(T); FORCE,NUMERICAL`
for every requested energy/force evaluation, then `PUT,XYZGRAD` to save
Molpro's native gradient XYZ. Each step retains its `.inp`, `.out`, `.xyz`,
launcher stdout/stderr, and a detailed Molpro `.log` requested with `-g`;
Sella retains its trajectory,
log, and final `.xyz`. Real licensed output remains required to validate
the documented units, sign, and atom mapping against the installed version.

Only after the accepted Molpro/Sella XYZ passes identity and termination checks,
the dispatcher releases independent Molpro harmonic, F12b/TZ, F12b/QZ,
CCSD(T)/DZ, CFOUR DBOC, and Gaussian VPT2
tasks. The VPT2 task performs its own B3LYP optimization before frequency
analysis, so its anharmonic correction is defined at a stationary point on
that surface. They each get their own isolated task directory and one-node Slurm
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
`docs/composite_qc_validation.md`, edit only site module setup/resources,
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

After CH4 dispatch passes, use the same graph and result schema in this
order. Each stage must fit practical cluster time and memory limits:

1. Parse and validate CH4 native output for every task, including termination,
   method/reference identity, F12b selection, frequencies, VPT2 zero-point
   correction, DBOC, and higher-order energy. Implement complete ANL0-F12
   expression and compare its components with a pinned reference. Connect
   the accepted molecular result to KinBot's existing energy and MESS data
   path, preserving a single consistent geometry and ZPE convention.
2. Run CH3 as the first open-shell molecule. It must test charge/multiplicity,
   ROHF/UHF decisions, Molpro open-shell coupled cluster, MRCC compatibility
   or direct fallback, and any spin-orbit term. Reuse the same offsite
   dispatch procedure and collect versioned fixtures.
3. Expand to a small set of stable molecules and radicals, then a small
   reaction. Prefer systems whose complete selected ANL ladder is reasonably
   quick on the target site. Add KinBot reaction routing, per-species
   provenance, restart across species, and final MESS handoff only after the
   single-species values are correct.
4. Add CBH-0/1/2/3 and higher connectivity-based reaction generation,
   chemically valid reference selection, atom and bond balance checks,
   deduplication, and uncertainty propagation. Let the user choose the ANL
   and lower-level methods for each rung. Compute heats of formation from
   the validated reaction ladder, then pass them to MESS with their source
   and uncertainty.
5. Integrate ATcT reference species for CBH-0/1 and other supported rungs
   through an officially supported machine-readable source if one is
   available. Store the ATcT network/version, species identity, phase,
   temperature, value, uncertainty, citation, and retrieval date. Cache a
   pinned snapshot for reproducibility; do not silently change a completed
   reaction result when the live database updates. The public ATcT site
   exposes versioned thermochemical tables, but API availability and usage
   terms need verification before implementing a live connector.

The dispatcher is a reusable staging and execution layer. The remaining
calculator wrappers, scientific output parsers, composite recipes, KinBot
result routing, multi-species resource accounting, and MESS integration
remain separate implementation gates. Later feature work must not treat
process exit status as scientific validation.

ATcT source: https://atct.anl.gov/

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
  reads the undisplaced total energy, and maps `PUT,XYZGRAD` rows back to
  ASE's coordinates. Molpro's F12b `ENERGY(2)`, `SCALE_TRIP=1`, and per-rank
  `-m` memory option were checked against the Molpro 2024 manual; direct MRCC
  `MINP`/`dmrcc` and CFOUR `ZMAT`/`GENBAS`/DBOC were checked against their
  manuals. These checks establish input plausibility, not runtime success.
- The direct MRCC XYZ block now includes the blank line required by its
  documented `geom=xyz` format and explicitly records frozen-core and
  spherical-basis choices. CFOUR must print its DBOC label, and MRCC must
  print normal termination without documented fatal/termination-error text;
  these are interface checks until full output parsers are validated.
- The earlier environment command could install conda-forge Sella 2.1.0 even
  though KinBot requires 2.6.0. The HPC runbook now pins ASE 3.29.0 and
  installs Sella 2.6.0 from PyPI before the editable KinBot install.

Run the exact cluster procedure in `docs/CH4_HPC_smoke_test.md`, using a
tested commit on the `composite` branch. First collect successful CH4 output
and method/version banners. Then compare scientific values and build real
parsers before enabling an ANL label or a MESS handoff. The local regression
suite passes without the four licensed codes; the external result remains the
gate for those program interfaces.

The first licensed Molpro force step must be inspected before releasing all
downstream jobs: compare its CCSD(T) energy and XYZGRAD forces with any
native `.log` or `.out` gradient table, verify the sign and units by a
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
  `FORCE,NUMERICAL` plus `PUT,XYZGRAD` inputs. The latter commands and the
  `XYZGRAD` force convention come from the Molpro manual and are exercised
  locally against a simulated executable. The first offsite Sella/Molpro
  step remains the acceptance gate for the actual Molpro 2024 `XYZGRAD`
  layout, force sign, units, atom order, and `-g` `.log` content.
- The Molpro 2024 parallel manual says one-node disk mode is the default and
  recommends `-m` without `-M`/`-G` there. Both the Sella calculator and
  the independent Molpro nodes now divide the requested memory into a
  bounded per-rank stack after reserving at least 200 MW per rank and node
  headroom. MPI ranks run with one OpenMP/MKL thread each.
