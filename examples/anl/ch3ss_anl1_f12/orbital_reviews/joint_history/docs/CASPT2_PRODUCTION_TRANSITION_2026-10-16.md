# CASPT2 production transition, 16 October 2026

The user selected CASPT2-F12 for both sampling and the multireference spin
gaps. The final MRCI build trial also failed, and the user closed further
MRCI work. Its historical evidence is preserved; it is not a production gate.
Minimum spaces remain CH3S C-S (5,4), CH3SS C-S (7,5), CH3SS S-S (7,5).
Keep IPEA=0.25 and set the real CASPT2 shift to 0.1 Eh. Preserve all native
evidence; do not relabel old shift-0.3 results as shift-0.1 calculations.

## Scope correction: no new two-bottleneck or VRC spin-orbit model

The proposed inner/outer example inputs and new automatic atomic spin-orbit
redistribution/offset code were withdrawn before submission. They were not
part of the agreed sulfur-example model. No VRC-TST flux calculation has been
performed by the fixed-geometry Molpro confirmation jobs. The confirmation
point label `inner` means the 2.1 angstrom geometry, not an established inner
kinetic bottleneck.

The NH2O paper explicitly introduces two VRC bottlenecks, rapid intersystem
crossing/equilibration, and region-dependent spin-orbit approximations for its
own reactions ([Methods 2.2](https://doi.org/10.1021/acs.jpca.3c05181)). That
is evidence for that published model, not sufficient justification to assign
the same model to CH3S/CH3SS automatically. The user's subsequent clarification
allows the existing RotdPy/MESS inner/outer capabilities and a general KinBot
hook to be used. Retain that optional route in the general interface. The latest instruction
selects regular single-bottleneck VRC for these examples, without inner/outer
stages or a rapid-equilibration assumption. Do not make rapid intersystem crossing
or region-dependent SOC general orbital-reconnaissance defaults. State averaging
for orbital optimization does not establish a two-bottleneck kinetic model.

Keep the existing, separately declared sulfur fine-structure data for
stationary-species thermochemistry and electronic partition functions. The
existing `term_centroid` convention gives the ground-J energy correction from
the measured levels; it does not specify a distance-dependent VRC potential,
state-specific capture weights, or an intersystem-crossing rate. No new
quenched/retained SOC treatment is installed in VRC. Production preparation
continues with CASPT2 sampling and the requested spin-gap/relaxation correction;
the removed draft surface files must not be described as production-ready.

### What was inherited and what was added

Comparison baseline: locally recorded `origin/master` at `5c6525a7`
(5 October 2026). The general electronic-state and multistate orchestration
modules were introduced on `composite` in `4793a0f6` (14 October 2026).

| Capability | Existing backend/master behavior | Work added on this branch |
| --- | --- | --- |
| Atomic electronic partition function | MESS already accepts `ElectronicLevels`; KinBot master normally wrote one level at zero with degeneracy equal to spin multiplicity. Master had a separate hardcoded OH fine-structure template. | General `electronic_state_models` JSON, species matching and forwarding of all declared levels to MESS/MESSPF. Sulfur values remain example data. |
| Spin-free energy to atomic ground-J energy | An electronic-level list specifies excitation energies, not an automatic shift of a supplied spin-free ground energy. | Explicit `spin_free_reference="term_centroid"` computes the degeneracy-weighted centroid correction and carries it into ANL species energies, with provenance and an existing-component check. |
| Two VRC flux calculations | RotdPy already integrates specified pivot surfaces, including COM pivots. MESS already supports `Variational`, `OuterRRHO` and `2TSMethod`. KinBot master already used `2TSMethod statistical` for a PST outer model plus an inner RRHO model. | General per-bond/per-stage/per-root orchestration, separate RotdPy tables and their handoff to MESS. This is not a change to RotdPy's integrator or MESS's two-bottleneck formula. |
| Automatic VRC SOC quenching/population redistribution | Neither the sulfur level list nor state averaging specifies this kinetic approximation. | The newly drafted helper and example settings were removed; this is not installed. |

The example sulfur levels are `(energy in cm^-1, degeneracy)` =
`(0,5), (396.055,3), (573.640,1)`, from
[NIST S I](https://physics.nist.gov/PhysRefData/Handbook/Tables/sulfurtable5_a.htm).
They give the atomic electronic partition sum directly. For the explicitly
declared spin-free term-centroid convention, the ground-J correction is
`-(3*396.055 + 573.640)/9 = -195.756111 cm^-1`, about `-0.560 kcal/mol`.
These are separate uses of the same data: the energy correction sets the
ground-state zero; the level list sets thermal excitation populations.
The separate `term_orbital_degeneracy=3` declares the three spin-free spatial
components for orbital/state planning. It is not a second factor multiplying
the resolved J-level degeneracies.

References for native capabilities: [MESS manual, sections 6.2 and 6.5](https://tcg.cse.anl.gov/papr/codes/mess/messmanual_3_23_16.pdf).
Code: `kinbot/electronic_states.py`, `kinbot/anl/automatic.py`,
`kinbot/mess.py`, and `kinbot/vrc_multistate.py`. Successful local rendering
tests are distinct from completed native VRC/MESS kinetics for these examples.

## Matched shift check: completed and independently reviewed

Both native jobs completed with all 24 correlated energies. The deposited
export passes inventory and native checks. The lowest D–Q gap changes by
0.255453 kcal/mol (aDZ), 0.264439 (aTZ), and 0.267019 in the requested
power-3.7 estimate. CASSCF/reference energies and active subspaces match.
See the [full review and reproducible comparison](CASPT2_SHIFT01_NATIVE_REVIEW_2026-10-16.md).
The historical submission commands below need not be run again.

The single geometry is the reviewed CH3SS S-S inner point (2.1 angstrom).
It has the largest shifted first-order total norm in the completed shared-
reference checks: 1.16660671 at aTZ and 1.15169295 at aDZ (SHIFT=0.3).
This is an evidence-based sensitivity proxy, not a proof of which uncomputed
geometry is globally most sensitive. Its fractional active occupations also
make it a useful interacting-state check. Run both bases at this one geometry
because the requested production correction extrapolates their spin gaps.

The plan keeps all six doublet and six quartet roots, the equal same-spin
reference density, IPEA, geometry, basis, core treatment and solver thresholds.
Each job starts from the same completed center checkpoint that fed its old
inner calculation. It makes no additional orbital permutation. The source
is the stopped recovery campaign, independently bound to the deposited audit.
No other campaign is stopped, edited or inspected for cancellation.

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/molpro_confirmation/check_caspt2_shift_on_blodgett.sh
```

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/molpro_confirmation/run_confirmation.py monitor /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_shift01 --all --roots'
```

After completion, export and transfer it to `examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history`:

```bash
/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/ch3ss_anl1_f12/orbital_reviews/joint_history/molpro_confirmation/run_confirmation.py collect /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_shift01 --output /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_shift01.evidence.tar.gz
```

Acceptance examines native completion, root energies/reference character,
first-order norms, selected orbital subspaces and changes in the spin gaps
at each basis and in their CBS estimate. A one-geometry check does not validate
shift independence over a whole PES. Production builds its own shift-0.1
asymptotes and correction scans consistently; it must not mix them with old
shift-0.3 energies. Resources: two jobs maximum, eight ranks each, 1000 MW/rank,
104000 MB/node; existing eleven exclusions inherited.

## Production handoff prepared

The [production package](../../../../ch3s_anl1_f12/orbital_reviews/RUNNING.md) supplies
six basis-specific checkpoints, separate per-bond orbital/state models,
CASPT2-only power-3.7 spin gaps, ordinary single-bottleneck VRC inputs and
the two ten-allocation launch helper. The general lowest-quartet selector
evaluates the entire declared manifold with common reference weights; target
geometry scans retain explicit doublet roots. This handles the observed
correlated-root ordering assumption without claiming universal state tracking.

The final local pass has 571 passing tests and 31 subtests, plus five separate
native RotdPy/MESS/MESSPF checks. Documentation includes every explicit input
parameter and its master/branch provenance. Scoped old-run retirement occurs
only after preparation and on-cluster preflight succeed. Full sulfur kinetics,
their scientific acceptance and deferred repository cleanup remain pending;
see the [maintained checklist](../../../../../../docs/ORBITAL_RECON_DELIVERABLES.md).
