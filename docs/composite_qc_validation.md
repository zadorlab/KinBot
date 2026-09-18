# Composite QC interface validation

The `composite` branch can be developed and tested on a machine without
Gaussian, Molpro, CFOUR, or MRCC. A licensed external site must run the final
program checks before any ANL method is described as validated. The local
checks and the external checks should use the same commit and saved task specs.

## Local environment

Use a Miniforge Python 3.11 environment. The FairChem extra is needed only for
UMA integration tests; the core interface tests do not need it.

```bash
mamba create -y -c conda-forge -p .venv python=3.11 pip numpy scipy ase=3.29.0 pytest networkx rmsd openbabel
.venv/bin/python -m pip install 'sella==2.6.0'
.venv/bin/python -m pip install -e . --no-deps
MPLCONFIGDIR=/tmp/kinbot-mpl .venv/bin/python -m pytest -q --ignore=tests/test_kinbot.py
```

Keep `.venv` out of Git. Record `python`, ASE, Sella, and KinBot versions in
external validation reports. No executable or model checkpoint path belongs in
the repository. The excluded legacy `test_kinbot.py` imports a hard-coded
development path and a top-level `thread_kinbot` module during collection;
it cannot be collected in a normal checkout. The new dispatch tests are run
by this command.

## First offsite test: general dispatcher with CH4

The first test runs one molecule through the molecule-independent
`kinbot.anl.dispatch` graph. The CH4 fixture is in
`examples/anl/ch4_dispatch.py`; no CH4-specific branch exists in the
dispatcher. It exercises these gates:

1. Gaussian B2PLYP-D3(BJ)/cc-pVTZ optimization through KinBot's ASE
   calculator and Sella.
2. Molpro CCSD(T)/cc-pVTZ geometry optimization through the new ASE calculator
   and Sella, starting from the accepted Gaussian XYZ. Molpro computes
   `FORCE,NUMERICAL` at each Sella step. The calculator reads the printed
   numerical gradient from `.out`; `PUT,XYZ` writes native coordinates used
   to map the gradient rows to ASE atom order. Retain each step's `.inp`, `.out`, `.xyz`, any
   `.log`, and Sella's trajectory and final `.xyz`. The calculator invokes
   Molpro with `-g` so each force evaluation has a detailed `.log`.
3. After that XYZ passes atom-order, finite-coordinate, and Molpro termination
   checks, stage the independent harmonic, F12b/TZ, F12b/QZ, conventional
   CCSD(T)/DZ, CFOUR HF/cc-pVTZ DBOC, and Gaussian VPT2 tasks. The revised
   higher-tier fixture runs VPT2 at B2PLYP-D3(BJ)/cc-pVTZ on the accepted
   L2 geometry, with no Gaussian `Opt`; it still joins the post-L3 fan-out.
   Up to `max_nodes` run concurrently.

The CH4 graph is an **interface smoke test**, not a complete ANL0-F12
electronic-energy expression or a MESS input. It samples Gaussian, Molpro,
and CFOUR before full recipe arithmetic, MRCC, open-shell CH3, and reactions.
For newly prepared runs, `result_parser` also checks the named Molpro,
Gaussian, and CFOUR output values. The completed first run predates those
declarations, so its native outputs were parsed read-only. A task marked
`complete` is still a dispatch result; any `review_required` parser flag
must be resolved before thermochemical use.

For the complete clone, setup, preflight, submission, monitoring, retry, and
result-collection commands, use [the CH4 HPC runbook](CH4_HPC_smoke_test.md).
Prepare on the cluster: generated scripts pin the Python executable used for
preparation and absolute run-directory paths. The dispatcher checks that
prepared inputs and accepted outputs have not changed before releasing
dependent jobs. A failed attempt can be archived with `retry` after fixing
site setup; changing chemistry inputs requires preparing a new run.

Inspect each task directory's exact input, `execution.json`, vendor output,
stdout/stderr, and Slurm logs. Confirm the printed method, basis, reference,
charge, multiplicity, F12b variable, normal termination, atom order,
frequencies, and DBOC value against the vendor output. Preserve program
versions and output hashes. Each Molpro ASE force step captures launcher
stdout/stderr separately from its native `.out`; external single points use
`launcher.stdout`. The Molpro command uses per-process `-m` after reserving
the manual's 200 MW per MPI process plus node headroom. The calculation is
single-node disk mode, which Molpro 2024 uses by default; its manual says to
avoid `-M`/`-G` in this mode because they can preallocate unused GA memory.
`-m` multiplies across `-n` processes. See its
[memory allocation rules](https://www.molpro.net/manual/doku.php?id=general_program_structure)
and [parallel memory guidance](https://www.molpro.net/manual/doku.php?id=running_molpro_on_parallel_computers).

Inputs are grounded in the [Molpro XYZ/NOORIENT and PUT,XYZ](https://www.molpro.net/manual/doku.php?id=molecular_geometry),
[FORCE,NUMERICAL](https://www.molpro.net/manual/doku.php?id=energy_gradients),
[FREQUENCIES](https://www.molpro.net/manual/doku.php?id=harmonic_vibrational_frequencies_frequencies),
and [F12 variable](https://www.molpro.net/manual/doku.php?id=quickstart)
examples, plus the [CFOUR Cartesian ZMAT](https://cfour.uni-mainz.de/cfour/index.php?n=Main.MolecularGeometryInput)
and [DBOC input and output example](https://cfour.uni-mainz.de/cfour/index.php?n=Main.CalculationOfDBOC)
pages. The completed CFOUR run now confirms the DBOC label and its HF value.
MRCC is deferred until its executable is available on the target site.
Licensed Gaussian confirmed the CH4 Sella force route and VPT2 execution;
its anharmonic correction still requires the convergence review below.

The local Molpro 2024.1 CH4 `DZ-F12/work/CH4/c000` reference supplied with
the peroxy workflow was inspected. Its native `OPTG` input confirms the
count/comment/coordinate XYZ block, charge and spin settings, RHF, and
correlated-method syntax. Its `.out` confirms the version, energy label,
normal-termination marker, and the location of long output in `.log`; the
`.log` contains numerical `dE/dx` rows, and `SAVEXYZ` plus `PUT,XYZ` produce
the expected numbered/final XYZ files. That calculation is
UCCSD(T)-F12b/VDZ-F12, whereas this first Sella calculator requests
conventional CCSD(T)/cc-pVTZ. The supplied tree has no standalone
`FORCE,NUMERICAL` example. The second Molpro 2024 CH4 attempt on Blodgett
computed CCSD(T) energy and the full `Numerical gradient for KB_GEOM_ENERGY`
table in `.out`, then failed at `PUT,XYZGRAD` because Molpro reported the
gradient unavailable for saving. The calculator now reads that table in
Hartree/Bohr and uses `PUT,XYZ` for atom-order mapping. The live output also
entered mppx mode with 11 processes despite the requested `-n 8`. The third
attempt set `MPPX=0` and completed the Sella geometry in two force evaluations.
It produced native `.out`, `.log`, and `.xyz` files with normal termination
and a final maximum atomic force norm of 0.01893 eV/Å. `MPPX=0` stopped the
mppx switch, but the third attempt still reported 11 compute processes plus
one helper. The Blodgett Molpro launcher drops its `-np` argument in Slurm,
so Intel Hydra used the Slurm task layout, which had only
`--cpus-per-task=8`. Prepared Molpro scripts now use `--ntasks=<ranks>` and
`--cpus-per-task=1`, and Molpro's default mppx numerical-gradient path has
been restored. A small live DZ job must verify the new count before the larger
jobs are released. The runtime rejects an output whose reported total MPI
process count exceeds the declared allocation.

## Local test matrix

| Layer | Test without QC executables |
|---|---|
| Profile and factory | Preset and override precedence, lazy optional imports, capability checks, per-job directory and command selection. |
| Input renderers | Compare exact generated input for closed-shell and radical examples; verify charge, multiplicity, atom order, units, basis, method, reference, frozen-core setting, F12 policy, and requested property. |
| Program runners | Mock process execution; assert the actual argument vector, working directory, fixed filenames, stdout/stderr destinations, exit-code handling, and isolation between concurrent tasks. |
| Output parsers | Parse versioned, redacted real output fixtures; test the intended energy or property line, normal termination, failure and truncation, duplicate markers, Fortran exponent notation, units, force sign, and atom order. |
| ASE and Sella | Use analytic toy ASE calculators for minimum, TS, constrained dihedral, atom, and diatomic cases. Check Sella gets energy and forces from the selected calculator on every step. |
| KinBot workflow | Mock the calculator backend and scheduler. Test L1/L2/HIR routing, coherent geometry/energy/frequency/Hessian selection, manifests, restart in a new process, and no legacy behavior change when the flag is off. |
| Composite | Use synthetic component energies to test every recipe equation, duplicate-node reuse, exact/custom labels, complete E0, and no second ZPE addition. |
| L3 scheduling | A mock scheduler proves no child starts before accepted geometry, all ready single-point/harmonic/VPT2 nodes can overlap on distinct exclusive nodes, global node/core/memory limits are never exceeded, and restart does not resubmit completed nodes. |

Fixtures should be collected from successful licensed runs and kept only if the
site permits redistribution. Otherwise, keep short redacted excerpts with
recorded program version and a SHA-256 of the complete private output. A
synthetic fixture does not establish that a program output parser is correct.

## Findings that constrain wrapper design

* ASE `GenericFileIOCalculator` calls `write_input`, `execute`, and
  `read_results` in that order; `BaseProfile.run` launches an argument vector
  in the calculator directory and captures stdout. This is verified against
  the installed ASE 3.29 source and its
  [calculator interface](https://ase-lib.org/ase/calculators/calculators.html).
  The wrappers for Molpro, CFOUR, and MRCC must supply their own templates.
* [Molpro's runner](https://www.molpro.net/manual/doku.php?id=running_molpro)
  reads an input file and normally writes a separate `.out` file. Its ASE
  template must capture launcher stdout separately and read Molpro's actual
  output file; pointing `BaseProfile.run` stdout at the same `.out` path risks
  overwriting the program output. Use a fixed safe input basename inside each
  task directory. In KinBot's ASE/Sella path, Sella owns the optimization and
  Molpro supplies per-step energies/forces. The input saves `ENERGY` into
  `KB_GEOM_ENERGY` immediately after `CCSD(T)`, then uses
  `FORCE,NUMERICAL,VARIABLE=KB_GEOM_ENERGY,STARTCMD=RHF` so the force and
  parsed energy refer to the same method. It reads Molpro's printed numerical
  `dE/dx`, `dE/dy`, and `dE/dz` table in Hartree/Bohr, then returns its
  negative in eV/Å as ASE forces. `PUT,XYZ` saves coordinates according to
  the [Molpro geometry manual](https://www.molpro.net/manual/doku.php?id=molecular_geometry).
  The calculator maps gradient rows back to ASE atom order by element and
  coordinate in that XYZ, reads the first `SETTING KB_GEOM_ENERGY` emitted after the
  undisplaced CCSD(T) calculation and before the numerical displacements,
  and requires normal termination. Keep Molpro's `-g` `.log` and Sella's
  trajectory and final XYZ. A separate native
  `OPTG` reference may be compared later, using Molpro's documented `.log`
  and `SAVEXYZ` artifacts; it is not the geometry driver here.
  [Molpro's F12 documentation](https://www.molpro.net/manual/doku.php?id=explicitly_correlated_methods)
  distinguishes F12a/F12b and limits analytic gradients to supported DF
  approximations. [The energy-gradient documentation](https://www.molpro.net/manual/doku.php?id=energy_gradients)
  requires `FORCE` after the corresponding energy calculation. A scaled F12b
  energy must never be silently paired with an unscaled analytic gradient.
* [CFOUR runs `xcfour` in the directory containing `ZMAT` and `GENBAS`](https://cfour.uni-mainz.de/cfour/index.php?n=Main.RunningCfour).
  Each ASE call therefore needs its own directory. Its
  [keyword documentation](https://cfour.uni-mainz.de/cfour/index.php?n=Main.ListOfKeywordsInAlphabeticalOrder)
  limits `DBOC=ON` to HF and CCSD with RHF/UHF references. Closed-shell
  CCSDT(Q) can use CFOUR's `xncc` module according to its
  [CC modules documentation](https://cfour.uni-mainz.de/cfour/index.php?n=Main.CCModules);
  the external test must identify the installed module/version and distinguish
  this from `CC_PROG=MRCC`. DBOC is a separate property node, not the
  CCSDT(Q) energy. The first Blodgett CH4 attempt found the site's `xcfour`
  needed `libgfortran.so.4` before it could read `ZMAT`. KinBot now checks
  ELF dependencies and supplies a validated matching Fortran runtime only to
  the CFOUR child when the site module does not already provide one. The
  retry then reached `xjoda` and revealed that CFOUR 2.1 truncated a long
  `*CFOUR(...)` line at 80 columns. KinBot now stages one keyword per line,
  following the [official continuation syntax](https://cfour.uni-mainz.de/cfour/index.php?n=Main.CfourKeywordSection),
  and rejects a keyword that cannot fit. The next retry must still verify
  actual native DBOC output and its numerical value. The subsequent CH4 run
  completed with HF DBOC `0.0025887093` Hartree and a separately printed MP1
  result. The result parser selects the named HF summary, cross-checks the
  printed units, and never substitutes the final energy including DBOC for
  the correction.
* [MRCC's manual](https://www.mrcc.hu/MRCC/manual/pdf/manual.pdf) specifies a
  `MINP` file in the run directory, invoked with `dmrcc`; it documents
  `CCSDT(Q)` and `CCSDTQ(P)` as distinct `calc` options and `geom=xyz` with an
  explicit coordinate unit. Use method-specific output parsing and pin the
  reference (`scftype`) in the task spec. Verify output termination and the
  exact final-energy label on a real fixture before enabling a recipe node.
* KinBot's Gaussian wrapper already writes `.com` and parses `.log` through
  ASE's Gaussian IO. The local renderer verifies the B2PLYP route,
  `EmpiricalDispersion=GD3BJ` representation, and charge/multiplicity. The
  licensed CH4 run confirmed the L2 force route and VPT2 execution. Native
  Hessian integration and VPT2 convergence remain separate validation gates.

## Existing KinBot interfaces to compare

* `kinbot/tpl/ase_sella_opt_well.tpl.py` already attaches KinBot's Gaussian
  calculator to ASE `Atoms`, runs Sella (or BFGS for a diatomic), and exports
  an XYZ trajectory. Preserve this behavior as the L2 reference when adding
  profile routing. Its broad exception catches are not a success signal;
  the new path should require validated result data and provenance.
* `kinbot/tpl/rotdPy_calc.tpl` configures rotdPy with `code='molpro'`, method,
  basis, memory, and queue. This proves KinBot passes a Molpro request into
  rotdPy, but the rotdPy package is external to this repository. Obtain one
  rotdPy-generated Molpro input and output on the external site and compare
  syntax and energies with KinBot's new renderer before using it as evidence.
* `kinbot/molpro.py` and `kinbot/tpl/molpro.tpl` provide a working legacy
  single-point pattern and parse `SETTING <key>` from Molpro `.out` files.
  The legacy template's generic `CCSD(T)-F12` call stores `energy(1)` as
  `mytza`; Molpro's [documented F12 variables](https://www.molpro.net/manual/doku.php?id=quickstart)
  make this the F12a result, while `energy(2)` is F12b. Keep the legacy path
  stable, but render and parse the requested F12b value explicitly for ANL.
  The legacy parser does not check normal termination, so it must not be
  reused as the ANL completion check.
* The current `slurm_molpro.tpl` has no `--exclusive` directive, and the
  generic `slurm.tpl` executes a shell command before the Molpro body. A new
  L3 script generator must put `#SBATCH --exclusive` among the header
  directives, before any executable shell line; appending it to the body
  would not create an exclusive allocation.

## MRCC driver policy

Prefer the Molpro→MRCC interface for a task only when the installed versions
support its exact reference and correction method. The
[Molpro MRCC chapter](https://www.molpro.net/manual/doku.php?id=the_mrcc_program_of_m._kallay_mrcc)
shows a closed-shell `mrcc,method=ccsdt(q)` input and example energies. It
also states that its interface is serial and restricts perturbative methods
to closed-shell. The same chapter has an O2 open-shell CCSDT/CCSDTQ example,
while the [MRCC manual](https://www.mrcc.hu/MRCC/manual/pdf/manual.pdf)
describes Molpro interface support for RHF, UHF, ROHF, and MCSCF orbitals.
These statements do not establish open-shell CCSDT(Q) or CCSDTQ(P) support
through Molpro on a particular installation. Test those exact methods and
their output labels on the external site before selecting that driver.

If a Molpro→MRCC task is unsupported, use direct `dmrcc` with explicit
`MINP`, or a separately validated CFOUR→MRCC driver, recording the driver in
provenance. Never silently switch the reference, method, or correction variant
when changing driver. A serial Molpro→MRCC job still receives its own
exclusive node because all L3 tasks do.

These are interface facts and design constraints. They do not substitute for
real execution or establish published ANL component agreement.

## External-site procedure and gates

1. Check out the tested commit and create the same environment. Record the
   program versions, module names, `PATH`, model checksum, and executable
   availability. Run each vendor's own small example before KinBot. Never
   commit licenses, executable paths, or private output.
2. For each backend, generate one input with KinBot and run it in a fresh
   directory. Save the exact input, stdout/stderr, program-created output, and
   exit status. Compare the parsed energy in Hartree against the program's
   printed value. For gradients, compare every Cartesian component and sign
   after conversion to ASE eV/Å.
3. Molpro: test RHF CCSD(T), ROHF/UCCSD(T), F12b with and without
   `SCALE_TRIP=1`, the permitted F12 derivative path, core-valence, and DKH.
   Confirm the selected F12b energy variable and normal termination. For a
   direct native `OPTG` reference, retain `.out`, `.log`, and `SAVEXYZ` output;
   compare its final geometry with the ASE/Sella result at identical theory.
4. CFOUR: test a separate DBOC case at a documented supported level and a
   closed-shell CCSDT(Q) case. Confirm `ZMAT`/`GENBAS`, `xcfour` invocation,
   module selection, and the reported correction or method-specific energy.
5. MRCC: test open-shell CCSDT(Q), closed-shell CCSDTQ(P), and open-shell
   CCSDTQ(P) with explicit references. Check exact output labels and
   convergence before using any value in a difference expression.
6. Gaussian: test B2PLYP-D3(BJ)/cc-pVTZ energy and forces, a native Hessian,
   constrained HIR through Sella, and a separate VPT2 job.
7. Re-run the same inputs through the KinBot ASE wrapper and compare every
   component with the direct program run. Then run L1→L2→HIR and a restarted
   mixed-backend workflow. Run each ANL variant against literature component
   values, including one closed-shell and one open-shell species. Record the
   geometry, basis, reference, correlation options, units, tolerances, and
   provenance for every comparison.

The feature remains experimental until all applicable external-site gates
pass. Missing licensed output fixtures or unavailable program capabilities
must cause a clear unsupported-feature error rather than an assumed energy.

## L3 dispatch and resource limits

After an ANL geometry result is accepted and its geometry hash is persisted,
release every independent single-point, harmonic-frequency, and VPT2 node.
The VPT2 task must contain its own prerequisite frequency work so it can run
alongside the harmonic node. A recipe may add a true dependency only when its
method requires one. Evaluate the composite expression only after all required
nodes complete. A changed geometry invalidates every dependent node.

Every L3 task, including geometry, harmonic frequency, VPT2, and each single
point, must run on its own exclusive node allocation. For Slurm, emit
`#SBATCH --exclusive` for every L3 job; do not batch independent nodes into
one shared-node job or omit exclusivity when a task requests fewer cores than
the node has. Use one node per task initially. A different queue system must
have a verified equivalent exclusive-allocation directive, or L3 submission
fails preflight. KinBot's `queuing=local` mode only reads existing results.

The user supplies the global maximum number of concurrent exclusive L3 nodes.
For automatic resources, preparation discovers the selected partition's node
CPU and memory through `sinfo`, sizes from the smallest eligible node, and
requests all node memory with `--mem=0`. Molpro's default rank selection
reserves 200 MW program overhead plus at least 1024 MW stack per rank and 15%
node headroom, then caps rank count at 16 to avoid blindly using all cores on
a memory-bound coupled-cluster calculation. The count is selected from
1/2/4/8/12/16. Task-specific `max_cores` and `min_stack_mw` permit a
benchmark-based policy; the default is a conservative starting point, not a
claim about optimal scaling for every system. Explicit per-node limits remain
optional for controlled fixtures. Reserve a node slot before submission and
release it on completion or failure. Count submitted and running jobs against
the node pool across all species in the KinBot run.
Persist reservations and job IDs so a restart recovers in-flight work without
duplicate submissions. Test this with a fake scheduler before using a cluster
queue, and assert the `--exclusive` directive in generated Slurm scripts.

## L2, hindered rotor, and VPT2 surface policy

For the lower ANL tier, use B3LYP/cc-pVTZ for all three. For the higher tier,
use B2PLYP-D3(BJ)/cc-pVTZ for all three. Hindered-rotor optimizations and
scans use the resolved L2 profile. A frequency-only VPT2 task consumes an
accepted L2 geometry at the identical method, basis, and dispersion level;
the dispatcher checks this relationship before staging. Gaussian's `Freq`
keyword alone does not optimize. The first completed CH4 VPT2 job *did*
optimize because its input explicitly contained `Opt=(Tight,CalcFC)`; the
revised fixture removes `Opt` and records that distinction in parser output.
The revised higher-tier L2 optimization uses Sella with a tighter force
threshold and an UltraFine grid before VPT2 uses its coordinates. The
published ANL1 equation's B3LYP anharmonic term and the user-selected
B2PLYP higher-tier variant must carry distinct provenance.

[B2PLYP analytic second derivatives and VPT2](https://pubs.acs.org/doi/10.1021/ct100212p)
have been reported in the primary literature; a later
[Gaussian-based VPT2 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8280743/)
also treats B2PLYP with D3(BJ). Exact Gaussian 16 input/output for the
revised higher-tier route still needs licensed-site validation when that
method is used. No new CH4 run is requested now.

## Completed CH4 native-output audit

The first Blodgett run finished all eight dispatch tasks. The text-only
archive preserved the native inputs and outputs, `execution.json`, and Slurm
records. Every archived file that appears in an execution artifact map matches
its SHA-256 value. The archive intentionally omits CFOUR's binary scratch,
`GENBAS`, and ASE trajectory files; the full run remains on the cluster.
These results validate the **CH4 interface smoke test**, not a complete ANL
energy or a MESS file.

| Task | Native-output finding |
| --- | --- |
| Gaussian L2 / Sella | B2PLYP-D3(BJ)/cc-pVTZ geometry converged; final B2PLYP energy `-40.487174956111` Hartree. |
| Molpro L3 / Sella | CCSD(T)/cc-pVTZ geometry converged in two Molpro force evaluations; final energy about `-40.43809881` Hartree. Both `.out` and `.log` are retained for each step. |
| Molpro harmonic | CCSD(T)/cc-pVTZ numerical frequencies have nine positive vibrational modes, six zero rotation/translation modes, and ZPE `0.04479801` Hartree (`9832.03 cm-1`). The reference-geometry gradient norm was `0.47427D-03` Hartree/bohr. |
| Molpro F12b/TZ | `SCALE_TRIP=1`, `ENERGY(2)`, `cc-pVTZ-F12`; exact F12b total energy `-40.454906199189` Hartree. The separately printed F12a energy must not be selected. |
| Molpro F12b/QZ | `SCALE_TRIP=1`, `ENERGY(2)`, `cc-pVQZ-F12`; exact F12b total energy `-40.456608306474` Hartree. |
| Molpro DZ | Conventional **CCSD(T)**/cc-pVDZ energy `-40.387076267138` Hartree. The old task ID `ccsdt_dz` did not perform CCSDT. |
| CFOUR DBOC | HF/cc-pVTZ correction `0.0025887093` Hartree (`568.156016 cm-1`), with a separate MP1 value `0.0026718675` Hartree. The final electronic energy includes DBOC and is not the correction. |
| Gaussian VPT2 in the completed first run | B3LYP/cc-pVTZ **with explicit `Opt`** optimized on its own surface; harmonic ZPE `9783.68667 cm-1`, total anharmonic ZPE `9646.36482 cm-1`, correction `-137.32185 cm-1`. Three native `WARNING:` lines make this value **provisional** and it is not the revised higher-tier B2PLYP component. |

The Gaussian warnings concern the rotor/framework classification and two
unreliable cubic force constants. Its resonance analysis also found six
active 2-2 and seven active 1-1 Darling-Dennison resonances. The job
terminated normally, but this VPT2 correction needs a numerical convergence
review before it enters ANL thermochemistry. This is a quality flag, not
evidence that Gaussian failed to run. The user has not requested another CH4
smoke run. A future calculation must use the tier-matched L2 geometry and
compare the ZPE components and warnings before thermochemical acceptance.

Newly prepared external tasks can declare parsers for Molpro's exact
CCSD(T)-F12b and conventional CCSD(T) energies, Molpro's harmonic modes and
ZPE, Gaussian's named anharmonic ZPE components, and CFOUR's named DBOC
level. A missing label, inconsistent units/summary, wrong echoed method or
basis, nonpositive vibrational mode, or absent final termination fails the
task. Gaussian warnings are preserved as `review_required` in the parsed
result. Existing completed task records are immutable; inspect their native
outputs with `python -m kinbot.anl.results` after updating the branch.

The [Molpro F12 manual](https://www.molpro.net/manual/doku.php?id=explicitly_correlated_methods)
defines `ENERGY(2)` as F12b, and the
[Molpro frequency manual](https://www.molpro.net/manual/doku.php?id=harmonic_vibrational_frequencies_frequencies)
defines the numerical frequency/ZPE output. Gaussian's
[vibrational analysis guide](https://gaussian.com/wp-content/uploads/dl/vib.pdf)
describes harmonic and anharmonic frequency analysis; the warning assessment
above is based on the native CH4 output and remains an inference pending the
repeat calculation.
