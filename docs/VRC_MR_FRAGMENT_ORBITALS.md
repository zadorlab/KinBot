# MR fragment geometries and fixed-geometry pivot orbitals

The CH3S/CH3SS examples now use the reviewed dimer active space to optimize a
separated supermolecule at the configured VRC asymptote. Fragment coordinates
are extracted using the accepted reaction atom maps. There are no separate
fragment active spaces or isolated-fragment reconnaissance files. Native
orbital/state diagnostics still need review at long range; reusing a reviewed
space does not guarantee that every numerical solution has the intended character.

The earlier Molpro pivot implementation optimized isolated fragments with
UKS/B3LYP3. That implementation removed the Gaussian checkpoint dependency,
but did not provide MR fragment geometries. It remains available as the legacy
mode; it is not used by these examples after the targeted migration.

## Single-input controls

| Input | Default | Meaning |
| --- | --- | --- |
| `vrc_tst_fragment_geometry` | `legacy` | Set `asymptotic_scan` to extract fragments from the MR separated-supermolecule optimization. |
| `vrc_tst_scan_basis` | Existing configured basis | Basis for both the relaxed scan and the MR asymptotic optimization; examples use `avdz-f12`. |
| `vrc_tst_pivot_level` | `L2` | In asymptotic mode, inherit the L2 method, basis, backend, SCF and grid for an orbital-only single point. |
| `vrc_tst_pivot_level: custom` | Opt-in | Use `vrc_tst_pivot_backend`, `vrc_tst_pivot_method`, `vrc_tst_pivot_basis` explicitly. Supported cube backends: Gaussian and Molpro UKS/RKS. |
| `vrc_tst_asymptote_distance` | Existing configured value | Separation for the supermolecule limit; these runs retain 30 angstrom. |

All existing distance constraints remain. No C3v or Cs geometry constraint is
added. The asymptotic calculation uses the same basis-specific reviewed CAS
selection, state average and selected root as that scan leaf. Its saved orbital
record seeds the finite scan. Subsequent points continue the orbital record.
The asymptotic fragment geometries supply the frozen correction geometries and
the rigid fragments passed to RotdPy.

The examples retain CAS(5,4) for CH3S C-S and CAS(7,5) for both CH3SS channels.
There are three chemical channels and four scalar paths: CH3S doublet roots 1
and 2, CH3SS C-S doublet root 1, and CH3SS S-S doublet root 1. The reviewed
basis-specific indices, state weights and spin manifolds remain in the input.

Sampling and correction controls are independent of the geometry basis.
The examples retain MR sampling and aDZ/aTZ-F12 spin-splitting extrapolation
with exponent 3.7, IPEA 0.25 and shift 0.1. This change does not alter the
88-point scan grid, RotdPy flux target, CBH-ANL recipe, or electronic levels.

## Auxiliary orbitals

For these examples L2 is Gaussian B2PLYP-D3BJ/cc-pVTZ. Each polyatomic fragment
gets one fixed-geometry single point at its extracted MR coordinates. Only the
reference HOMO (alpha HOMO for an open-shell fragment) supplies pivot directions.
Its energy is not substituted into the MR correction or sampling potential.
No geometry optimization, Hessian, frequency job or separate fragment CAS is
requested. Monatomic fragments need no orbital job.

Gaussian `formchk` and `cubegen` run within the compute-node local scratch
workspace. Explicit outputs are copied back with provenance. Custom Molpro
UKS/RKS orbital-only jobs export native cubes at the same fixed coordinates.
The Gaussian L2 path rejects composite methods that implicitly optimize and
extra route instructions that would change geometry or request derivatives.

A cube and `.orbitals.json` belong to a channel/root/fragment, not just a
species name. The handoff checks the cube hash, atom order and geometry. A
rigid coordinate-frame rotation is allowed; a stale cube cannot replace the
MR geometry. Pivot directions are rotated into the accepted fragment frame.

The existing pivot ranges overlap: orbital-oriented pivots are used at short
range (default upper limit 6 angstrom), atom-centered pivots at 5-12 angstrom,
and center-of-mass pivots from 10 angstrom. Orbitals are calculated once for
each rigid fragment, not at every separation. This preserves RotdPy's existing
pivot construction; it does not introduce a new MR definition of a unique HOMO.

## Native cube evidence

The user executed the small Molpro H2/CAS(2,2)/cc-pVDZ signed-orbital probe on
Blodgett and reported:

```
PASS: native cube read, atom coordinates, signed field, phase invariance
Grid: (31, 33, 35)
Directions: [0,0,-1], [0,0,1]
/tmp/kinbot-molpro-eifnscg9/probe_orbital_2.1.cube
```

This validates the tested native Molpro cube layout against the installed
reader. It does not certify an entire production VRC run or a chemical orbital
assignment. Local tests additionally cover Gaussian orbital labels, bounded
EOF handling, exact voxel counts, Fortran numeric formats, full voxel vectors,
loop ordering, phase invariance, geometry/hash rejection and frame rotation.
Unsupported multiple-field or angstrom-grid variants are rejected explicitly.
See the [Molpro CUBE specification](https://www.molpro.net/manual/doku.php?id=dump_density_or_orbital_values_cube).

## Targeted controller handoff

Review without changing controllers or scheduler allocations:

```bash
python -u -m kinbot.run_refresh RUN \
  --vrc-scan-basis avdz-f12 \
  --vrc-fragment-geometry asymptotic_scan \
  --vrc-pivot-level L2 --cancel-superseded-vrc-scans
```

Add `--apply` to perform the handoff. It validates the reviewed new-basis
models and utilities first, quiesces only PES/KB controllers, inventories exact
dispatcher ownership again, and cancels only superseded `relaxed_scan` jobs.
Cancellation must be confirmed by terminal scheduler records before launch.
Other running and pending jobs, their task ledgers, completed ANL results and
the live shared node policy remain unchanged. This option cannot be combined
with `--replace-pending`.

Input changes and original bytes are journaled/archived. The manifest is
updated after the input edits. The same command can finish an interrupted
transaction without duplicate drivers or allocations. If Slurm has not yet
confirmed a cancellation, repeat the same command after that job terminates;
an empty queue alone is insufficient evidence.

Old scan directories stay in place, including trajectories, checkpoints and
native diagnostics. `vrc_restart.json` lists the retired graphs and scheduler
outcomes. They are comparison evidence, not production aDZ results. The new
basis/geometry inputs produce new dispatcher fingerprints. Cancellation frees
the nodes after Slurm cleanup and loses in-memory optimizer state. A future aTZ
retry must initialize from saved geometry/checkpoint data, not claim to resume
an interrupted Sella instruction. This handoff does not submit such a retry.

Only completed, validated points with matching channel, root, separation and
constraints should be compared. Partial optimizer steps are useful diagnostic
starting geometries, not converged reference geometries. Changing the fragment
limit as well as the scan basis means an old/new difference cannot automatically
be attributed to basis size alone.
