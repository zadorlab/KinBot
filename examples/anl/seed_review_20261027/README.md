# Retired r7 seed evidence, 27 October 2026

The uploaded `seed-evidence.tar.gz` has SHA-256
`e1d601d46211fc12106df565bc6fb6676dba1ec30c6ec30bed08f3665ad74732`.
All 4,628 captured files match the export manifest. The archive is an audit packet,
not a full calculation backup; binary checkpoint bytes remain in the verified
preservation directory on Blodgett. No source run was restarted or modified.

`report.json` records exact relative paths and hashes. Reproduce it from a safely
extracted packet with `python -B audit.py /path/to/extracted-packet` (KinBot on
PYTHONPATH). The script verifies every captured file before inspecting it.

## PES search history that must survive import

| Case | Accepted searches (-1) | Terminal failures/rejections (-999) | Unfinished constrained searches (0) | Product-stage searches (4) | Database rows |
| --- | ---: | ---: | ---: | ---: | ---: |
| CH3S | 29 | 35 | 0 | 0 | 8,089 |
| CH3SS | 131 | 234 | 101 | 6 | 11,658 |

All 18 worker databases pass SQLite `quick_check`. All 18 search ledgers match
their saved input and theory contracts. These are search-controller decisions,
not counts of unique accepted network reactions or independent calculations.
The 269 terminal failures/rejections must remain closed; the 160 accepted searches
must be reconstructed from retained results and reassessed under the new L1
policy without repeating TS searches. There are three pending-submission records
(in workers ending `s1a4d7914034ddafc`, `sc9062a2d39c6d51a`, and
`s1b0c8a00a2bddd95`) that require reconciliation with retired allocation evidence.
The importer must not clear these records indiscriminately.

Current graph records mark 1,121 CH3S and 205 CH3SS components complete. These
1,326 entries are reuse candidates, not yet validated imports: most complete
native artifacts are deliberately outside this packet. Existing counts include
VRC scan tasks as well as composite components. The older L1 report classifies
71 complete entries under the excluded CH3S worker and 1,255 as retained, shared,
or unknown; this classification is not a proof that all 1,255 are reusable.

## Long geometry calculations

Both use CCSD(T)-F12b/aug-cc-pVQZ-F12. Starting atom order and coordinates match
their staged geometry hashes. The selected frames agree with Sella's explicit
`trjid` and energy, not merely the last trajectory frame.

| Slurm job | Task | Saved evaluations | Last logged step/frame | Force maximum (eV/angstrom) | Target |
| --- | --- | ---: | --- | ---: | ---: |
| 52918345 | `ts_72ecf8afb235bd9e/l3_geometry` | 13 | step 2, frame 12 | 0.00847483 | 0.0005 |
| 52928790 | `ts_0450141d9a891e9c/l3_geometry` | 1 | step 0, frame 0 | 0.02305410 | 0.0001 |

The first contains two post-initial optimizer steps; ten other evaluations are
nearby Hessian probes. Frame 12 is a useful geometry-continuation candidate.
The second has an evaluated starting point only; the packet does not prove any
completed subsequent optimizer step. Both later requeue trajectories are empty.
Neither optimization is converged. Reusing these coordinates must preserve the
original scientific task, re-evaluate forces and rebuild the optimizer Hessian;
this evidence does not prove recovery of a partially evaluated Hessian.

The complete saved coordinates and source hashes are in `geometry_candidates`
in `report.json`. Each trajectory is under the corresponding task's
`.kinbot_slurm_io_retries/.../restart-1/interrupted/` directory.

## CCSDT(Q): preserve the older node3 checkpoint

Job **52929237** belongs to the CH4 fragment task
`ch3s/470831050000000000002/anl_automatic/states/state_684eb53c4637f96a/composite/tasks/ccsdtq_tz`.
The original interrupted attempt is attributable to
`/scratch/hram/kinbot-mrcc-btscbyz1/work` on **node3**, through its archived runtime
record. Its exact `MINP` is identical to the two requeued attempts:

- Input SHA-256: `62534e4e8d6b1ca1a26c7f97ad7abeecdea2f7c3992aefc8d355b7475f4f0484`.
- Saved `fort.16`: 475,979,776 bytes, SHA-256
  `665cf47f3fb8f3a187e6a08913c0398e062be16c8e4869c254d132c3c9b7e8ff`.
- Preserved path:
  `~/kinbot-evidence/r7-outage.OMot4WR0/checkpoints/nodes/node3-b0c30c75/files/df989e252cb534f3/work/fort.16`.
- Native output: CCSDT converged in 14 iterations at **-40.441308731182 hartree**;
  perturbative corrections reached spin case 3 before SIGTERM.
- The node4 and node35 requeues reached CC iterations 6 and 9 respectively.

The original node3 amplitude file is the preferred **candidate**, not a verified
native restart. There is no final CCSDT(Q) energy. We cannot infer preservation
of partial perturbative corrections from `fort.16` or the printed spin-case
progress; a backend restart may have to repeat that correction. The larger
808,812,544-byte node35 file belongs to a different C/H/S calculation and must
not replace this checkpoint.

## Migration readiness

The archive establishes that useful progress and failure history survived. It
does not itself implement an independent seed import or prove native restarts.
The importer still needs source immutability, compatibility/path rebinding,
pending-submission reconciliation, retained L1-policy reassessment, sampling
adoption, and interruption/repeated-invocation tests. The closure orbital panel
also remains pending; no earlier CAS selection is changed by this audit.
