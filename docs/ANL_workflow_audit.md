# CBH/ANL, VRC-TST, ROTD_py, MESS audit — 14 October 2026

**Work resumed.** See [the current implementation and acceptance procedure](WORKFLOW_READINESS_2026-10-14.md).
The earlier checkpoint is historical. Preserve successful native pilots and
the existing eight-basis campaign. The latest 5Z failure occurred at Molden
export after RHF completed; native XML/WFU recovery is published in `c14913f`.

The current software is **not yet certified for the full CH3S/CH3SS production
campaign**. The automated interfaces and diagnostic reconnaissance have been
extended and tested locally. Licensed native calculations and the chemistry
acceptance checks below have not been performed in this environment. There
is no local Molpro, Slurm, or MESS executable. SSH is handled by the user.

The latest instruction is authoritative: use Molpro's default F12 basis sets.
KinBot supplies the requested native basis name, without a forced revision,
vendored definitions, generated diffuse exponents, or auxiliary substitutions.
No large reconnaissance or production jobs have been submitted by this audit.

## Scope and compatibility

The review covers the current `composite` tree, its changes relative to local
`master` and `origin/master`, and the CBH/ANL → VRC/ROTD_py → MESS contracts.
The starting tree already contained substantial uncommitted implementation;
it has been preserved. A passing regression suite is evidence for the tested
contracts, not a proof that every possible difference from master is correct.

Legacy single-basis ROTD inputs, multiplicity-only MESS blocks, archived
singlet/triplet correction field names, and formation records lacking new
uncertainty evidence remain supported. Changed reviewed WFU records or changed
sampling inputs are deliberately rejected for cache reuse. No previously
generated stock-F12 calculation is rejected for lacking a revision tag.

## Findings addressed in this audit

| Area | Defect or missing behavior | Result |
| --- | --- | --- |
| VRC levels | Expensive CBS sampling and inconsistent orbital records | Preset samples at augmented TZ-F12; TZ/QZ-F12 CBS is retained for the correction/reference energies; expert overrides remain available |
| State routing | Connectivity-only S2 could be treated as a singlet and exported with an incompatible Lewis representation | User-declared electronic models apply before QC naming; validated same-connectivity aliases preserve the intended open-shell identity |
| General spins | Fixed singlet/triplet assumptions did not cover arbitrary fragment coupling | Spins follow accepted parent/products; automatic reference uses maximum spin; an already-high-spin target requires direct correction or an explicit expert reference; overrides cannot relabel the accepted channel |
| Spatial roots | Multiplicity alone omitted fragment orbital degeneracy | Root counts follow declared spin-free terms; CH3SS S–S cleavage retains six spatial product components when CH3S X(2E) and S(3P) are included |
| MESS states | Fixed one-level templates and spin-times-orbital degeneracy could not represent resolved J levels | Stable species support explicit level lists and independent 2J+1 degeneracies; reactive surfaces count their selected multiplicity once |
| SO convention | A zero correction had been described as a universally cancelling CBH gauge | Missing corrections are explicitly spin-free approximations; sourced corrections and complete-term centroid corrections propagate into ANL and fallback records |
| Active spaces | Counts alone did not identify the orbital subspace | Reviewed basis-specific WFU records, hashes, START/DIAB propagation, electron parity, available orbital holes, averaged spins, and selected roots are checked |
| ROTD sampling | Upstream worker could not evaluate a basis pair or the reviewed state model | KinBot adapter evaluates both native inputs, restores parent atom order, selects the requested native MR energy, and extrapolates in hartree |
| Provenance | Correction comparison omitted actual channel context and JSON changed root-key types | Parent atoms/charge/spin, fragment spins, roots, basis models and mappings are carried through the correction and ROTD boundary |
| Provenance | Basis-list mutation and changed WFU/input reuse could contaminate caches | Lists are copied; WFU hashes enter generated inputs; incompatible execution and sampling caches fail clearly |
| Output routing | ROTD generation ignored the configured output directory | The selected directory is honored |
| CBH | A strict-subset check could miss missing energies if extra entries were present | Every stoichiometric participant must have an energy |
| Uncertainty | New uncertainty requirements rejected otherwise valid older recursive records; CBS formation sensitivity had the wrong sign | Unknown uncertainty stays unknown; signed shared-source contributions retain cancellation without blocking the enthalpy solve |
| Exact states | A zero-electron atomic special case depended on a particular ChemID | Zero-electron atomic energy follows atomic number and charge |
| Slurm | Reconnaissance directives followed an executable shell line and could be ignored | All SBATCH directives precede shell execution; ranks, threads, memory and walltime are explicit |
| Runtime | New native probe/sample paths bypassed isolated Molpro scratch and MPI handling | They use the existing production runtime resolver and cleanup |
| Recon restart | Duplicate controllers, overwritten input fingerprints, and reassigned predecessors could invalidate evidence | A manifest lock, immutable input checks, submission journal and preserved completed checkpoint chains prevent these cases |
| Recon evidence | Normal termination alone did not prove completion of the requested diagnostic block | A final computed sentinel, all requested root energies, natural populations, nonempty artifacts and their hashes are required |
| Recon coverage | A few fixed distances, or a flat parent root alone, could hide an unsettled asymptote | Dense near-bond probes, plateaus for every requested root, and adaptive energy/root/occupation refinement replace the fixed three-anchor example |
| General input | Light-atom cleavage could omit its independent projection; ionic splitting could silently assume neutral fragments | Light-atom coordinates include their valence projections; explicit charged-fragment states are validated against connectivity, charge conservation and spin |
| Operations | User-wide kill patterns or stale PID files could stop unrelated work | A run-scoped cleanup command checks live Slurm WorkDir and Linux processes, protects the invoking shell, checks PID identity, and can archive a stopped tree |

The new `anl1-f12` VRC preset uses arithmetic averaging of the two gap models.
The legacy distance-switch policy remains readable/configurable; its arbitrary
weights are not claimed to be a universally validated chemistry model.

## What the reconnaissance now supplies

The examples request eight native orbital bases: augmented conventional D/T/Q/5Z
and augmented F12 D/T/Q/5Z. Arbitrary safe Molpro basis names may be supplied;
installed availability is decided by Molpro, not a KinBot list of elements.
The initial full grids contain 248 CH3S and 600 CH3SS jobs, each on one node.
There are another 24 native basis/channel pilots before the full scans.
The grids include independent rigid-fragment orientation probes at several
separations and at the initial asymptote.
The example submission limit is ten. Refinement can add points and stops at a
configured pass limit with an explicit unresolved-coverage status.

Coordinates are graph-defined separable bonds, or `coordinates: "auto"` for
all bridge bonds. An ionic channel can declare two `fragment_states` entries
with `atom` (a distinguishing parent atom), `charge`, `multiplicity`, and
`smiles`. Ring opening requires an explicit geometry/multibond definition;
the driver does not pretend that removing one ring edge separates fragments.

With `probe_active_space: true`, each geometry produces a full-valence
heavy-atom AVAS/CASCI diagnostic space, natural populations, all requested
root energies, orbital listings, native XML (or legacy Molden) and WFU records. A light-atom
coordinate also includes independent light-atom projections. `probe_space`
can request all light atoms or extra atom-labelled valence/polarization
centres. The report lists occupation candidates and electron populations per
geometry; sorted occupations do not establish orbital identity.

The parent spin is first in the probe's energy array. Each array entry is
labelled by multiplicity and root in the report. These are energy-ordered
labels, not proof of diabatic state continuity. Checkpoint propagation
is independent for each basis and channel. The controller submits a successor
only after its predecessor is complete and stages that exact WFU. Refinement
does not rewrite the history of completed probes.

## Chemistry acceptance still required

1. **Native pilot evidence.** Validate the installed Molpro version, all
   requested basis families, AVAS/CASCI input, multi-spin root arrays, WFU
   restart/DIAB behavior, MR-F12 energy selection, and numerical gradients.
   Confirm scheduler resources before a large batch. HF reconnaissance does
   not itself validate F12 correlation/CABS or MRCI/CASPT2 execution.
2. **Active-space convergence.** CAS(13,8) and CAS(9,7) in the production
   examples are candidate requests. Review occupied and weakly occupied
   valence orbitals, bonding/antibonding partners, diffuse/3d/Rydberg candidates,
   alternative subspaces and larger-space results. A natural occupation cutoff
   alone cannot certify the required active space. Keep consistent orbital
   character and selected roots across bases, separations and both CH3SS
   channels. WFU hashes establish provenance, not chemical correctness.
3. **Geometry and state coverage.** Confirm optimized parent and fragment
   structures, relevant conformers/torsions and fragment orientations, and
   energy/root/occupation convergence throughout the reaction region and at
   separation. The automatic grid covers radial separation at the supplied
   geometry; it does not exhaust all conformations or automatically prove
   wavefunction overlaps through crossings. Follow apparent state changes with
   denser probes and explicit overlap/CI-character analysis.
4. **SO and vibronic model.** Atomic S J-level populations and a sourced
   centroid correction are supplied. CH3S has Jahn–Teller and SO structure;
   its spatial-term declaration does not supply the resolved molecular
   partition function or a compatible molecular SO correction. S2 fine
   structure and molecular CH3SS corrections also need evidence. Do not treat
   a missing correction as an exact physical zero or double-count the same
   levels in electronic degeneracy and separate surfaces.
5. **State-resolved statistical kinetics.** General orchestration now prepares
   independent spin/root potentials and correction models using existing
   scalar ROTD_py integration, sums physical fluxes on aligned energy grids,
   and renders single or statistical/minimum inner–outer MESS bottlenecks.
   Stage-specific active spaces/averaged roots and explicit statistical levels
   are supported. State-average roots are not assumed reactive contributions.
   NH2O SI-derived contract checks and regression tests pass; further MESS
   integration auditing continues. The sulfur reactive subsets, energy offsets,
   and any equilibration assumptions require the reconnaissance evidence.
   Homolytic vdW-complex discovery is deferred under the latest user scope;
   multistate accounting does not force a complex into the network.
6. **Electronic correlation.** Check single-reference diagnostics and
   high-spin reference quality, CASPT2 intruders/level-shift sensitivity,
   MRCI+Q size consistency and separated-fragment gaps. Converge basis/CBS,
   fitting/CABS and correlation choices rather than interpreting a normal
   exit as an accuracy certificate. Every reduced ANL recipe remains labeled.
7. **Kinetic join.** Validate at least one native correction → ROTD sample →
   number-of-states → MESSPF/MESS handoff per channel. Confirm Monte Carlo,
   pivot, dividing-surface and energy/angular-grid convergence; electronic and
   optical statistical counting; zero-energy alignment; rotor/harmonic
   fallbacks; and pressure/temperature/collision-energy-transfer sensitivity.
   The example Lennard-Jones and collider settings are inputs needing physical
   validation, not fitted sulfur-system constants.

These checks remain before a chemistry production claim. No fully automatic
minimal-CAS certification is claimed. The reconnaissance can now run while
MESS work continues; the full KinBot production campaign follows review.

## NH2O literature cross-check

[Klippenstein, Mulvihill and Glarborg (2023), section 2.2](https://doi.org/10.1021/acs.jpca.3c05181)
([author manuscript](https://backend.orbit.dtu.dk/ws/portalfiles/portal/337661919/nh2_o_revised_final.pdf),
pp. 5–9) separates inner and outer bottlenecks. Inner reactive contributions
are two doublet states for NH2 + O and one for NH + OH; outer capture includes
all electronic states, under rapid long-range intersystem crossing. SO
splittings are retained at long range and neglected at short range. The
state-average and active-space choices also differ between inner and outer
calculations. Their comparison found spin splitting unnecessary for these
channels; it is an energy correction strategy, not a state-transition model.

This supports using existing ROTD_py integration with separate state-specific
runs and a statistical MESS network. Its vector energy arrays hold separate
flux estimates; KinBot's scalar extra indices represent correction variants.
Those variants must not be summed as if they were electronic states. Neither
the NH2O state counts nor its equilibration assumption should be copied into
the sulfur inputs without evidence. No changes to the pinned ROTD_py package
are required to evaluate the independent scalar fluxes.

The actual publisher SI MESS network and nine Ne tables were obtained and
checksum verified. [The derived report](NH2O_SI_verification.json) records
barrier/state/level counts and an isolated electronic-population integral
check (relative errors below 0.05%). No full native-rate reproduction is
claimed. The SI's negative-energy Ne thresholds are preserved and their
MESS ground shifts are included when explicit ZeroEnergy values are written.
Legacy origin insertion remains supported; a sampled origin is not mandatory.

## Acetaldehyde interruption

The user-provided monitor reports no live run jobs after a cluster restart,
but saved branches still read `running`. Those are historical state labels,
not live liveness evidence. The restart killed the driver/jobs before a final
status could be published. The HIR messages explicitly select harmonic
fallbacks; the vdW optimization is optional. Neither is evidence of a fatal
driver exception.

The last logged VRC action retried a failed first `relaxed_scan` attempt.
Its native cause cannot be established from the monitor alone. The archived
`execution.json`, launcher diagnostics, native outputs and Slurm accounting
are needed to distinguish node loss, MPI, scratch, convergence or parser
failure. Do not attribute that earlier attempt to the later restart without
timestamps/native evidence. See [Blodgett notes](Blodgett_workflow_runbook.md)
for the exact recorded run path and scoped stop/archive commands.

## Evidence and references

The final stopping-checkpoint suite passes **1231 tests and 391 subtests**, with two warnings
(an existing deprecated `imp` import and denied access to Sella's external
compilation-cache directory). `git diff --check` passes. The module-selection fix additionally has an
executed batch-shell regression test and the campaign has a no-submission
module/scheduler preflight. Tests use native-format fixtures/mock schedulers;
they do not substitute for licensed chemistry calculations.

An additional run of the **unmodified `origin/master` tests** against the
current source gives 602 passes, 41 failures, and 369 passing subtests. The
failures are fully accounted for by these older fixture assumptions:

| Older fixture | Count | Current requirement |
| --- | --- | --- |
| MESS execution fixtures create no input files | 27 | Actual input files are hashed before execution/reuse |
| HIR replacement workers supply no scheduler-generation completion marker | 12 subcases | Replacement completion must be bound to the submitted generation; the old Popen mock also cannot model subprocess.run |
| VRC selector-only fixture omits executable/pivot configuration | 1 | HOMO pivots preflight formchk/cubegen; geometric fixtures declare their intended mode |
| Species-routing fixture patches removed private which/Popen symbols | 1 | Native conversion uses resolved executable paths and subprocess.run, and requires resulting artifacts |

The current suite retains the corresponding behavior tests with complete
artifacts and updated native/scheduler mocks. The unmodified older suite is
therefore **not reported as passing**, and this audit does not claim universal
compatibility beyond the tested contracts. Rejecting unbound replacement
results and input-less MESS fixtures is intentional provenance enforcement.

- [Molpro MULTI, state averaging, natural and diabatic orbitals](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi)
- [Molpro native F12 basis and auxiliary selection](https://www.molpro.net/manual/doku.php?id=explicitly_correlated_methods)
- [CH3S multistate/SO study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5995305/)
- [CH3S spin–vibronic spectroscopy](https://www.sciencedirect.com/science/article/abs/pii/S0009261411008153)
- [NIST atomic sulfur levels](https://physics.nist.gov/PhysRefData/Handbook/Tables/sulfurtable5_a.htm)
- [Slurm SBATCH directive processing](https://slurm.schedmd.com/sbatch.html)

## Native Blodgett pilot feedback

The first CH3S pilot produced a Molpro `UNMATCHED }` error from an extra
ROHF closing brace. Local fixture assertions had checked only substrings
and missed it. The generator is corrected and every generated deck now has
a structural delimiter check before preparation/execution. The job command
also explicitly sets its WFU repository (`-W`) rather than inheriting
Molpro's observed `/home/hram/wfu` default. An executed batch-shell test,
all-example deck sweep, malformed-deck regression, and WFU-command test cover
these contracts. A new campaign first validates one native probe before
widening its queue. A native execution record reports Molpro errors directly.
The failed campaign must be archived and regenerated; it is not relabeled as
a successful native test. Native success remains to be demonstrated.

The replacement Blodgett augmented-DZ probe has now completed natively with
all six requested states and all three orbital artifacts. Its first audit
was a false negative because Molpro prints `HARTREE` instead of `AU` for
the saved variables. Readers now accept both atomic-energy labels, with a
regression based on the supplied native root energies and occupations.
The corresponding ASE/VRC reader and ANL F12 validation were corrected too;
consistent assignment/SHOW duplicates are allowed while mismatches still
fail. This does not validate other bases, geometries or production correlation.

The final stopping checkpoint adds a scheduler-checked `retry-failed` command
that archives known failed native attempts while retaining successful probes.
The user supplied a compute-node `/scratch` mount with about 5.4 TB free;
`KINBOT_MOLPRO_SCRATCH` selects it without embedding a site path in code.
The former required/multistate full-PES guard has now been replaced by shared
direct/PES staging: hash-bound worker/UQ convolution contexts are resolved
after global ROTD execution and final network energies. Both scalar and
multistate actual-assembly tests cover this path. Unselected optional exits
retain PST, and old required worker artifacts lacking context fail with a
regeneration instruction. Workers do not execute a duplicate ROTD campaign.

The resumed audit also preserves reviewed reference-root selections, rejects
an excited CAS root for a ground-state-only SCF/CC anchor, sets charge/spin
before AVAS, replaces production VRC Molden exports with XML after CASSCF,
makes the production asymptote configurable, reports MR method disagreement,
and replaces stale/per-UQ ROTD provenance on model rewrites. True fragment
conformer ensembles require separate sampled fluxes and remain future scope;
explicit ROTD requests reject that combination instead of miscounting it.
