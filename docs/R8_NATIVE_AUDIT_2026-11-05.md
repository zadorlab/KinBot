# R8 native evidence, November 5

This is an audit of the newly supplied export, not a declaration that both PESs
or their MESS models are final. No QC calculations were launched or reset.

## Capture and reproducibility

- Archive: `examples/anl/tmp/seed-evidence.tar.gz`.
- SHA256: `df4c2035dcbe174f81abdceb4ff6f27ad3f1a929fe03e10beb9050c873b50038`.
- All 216,711 captured files passed checksum verification; the manifest reports
  zero collection issues and 18,579,805,441 bytes of source files.
- Preserved extraction: `examples/anl/tmp/r8-native-20261105/`.
- Local read-only analysis: `local_runs/r8_native_audit_20261105/native_audit.py`
  and `check_models.py`; outputs `native_audit.json` and `model_checks.json`.
  SQLite was opened with `mode=ro`. Database energies were converted from eV
  to hartree; saved ZPEs and exported reverse-result energies are already in
  hartree. Geometry comparisons use element-preserving permutations and proper
  rotations, not reflections.

The capture is scientific evidence, not a complete restorable execution tree:
the pool-root shared-QC registry is not included. Producer/consumer receipts,
native files, and pending worker handoffs inside the ordinary workers are present.

## Completion and homolysis

CH3S records successful PES completion at 2026-11-05T16:34:29Z, including its
combined L2 MESS input. CH3SS records a running controller with no recorded
termination or traceback; its worker logs stop around 16:34:43–48Z. The earlier
remote monitor reported no live controllers. The export does not establish the
cause of their disappearance. Three CH3SS optimization-monitor JSON files are
empty; the monitor writer has already been changed to durable atomic writes.
That observation is not proof that those files caused the process exits.

CH3S's ten homolysis checks record six terminal nonfindings, three superseded
parents, and one unretained parent. These are not ten failed QC calculations.
CH3SS's saved check ledger contains nineteen terminal nonfindings, six
superseded parents, seventeen older `unresolved` entries, nine pending entries,
and five unretained parents. Current search state, rather than the stale ledger
alone, has six active companion searches:

| Worker suffix | Companion bond | Saved stage | Native evidence |
|---|---|---|---|
| `s893815ff03c7331d` | 1–2 | 1, step 41 | Normal L1 and L2 records; reconciliation remains |
| `sedd8f7729b69b19b` | 1–2 | 1, step 41 | Normal L1 and L2 records; reconciliation remains |
| `sedd8f7729b69b19b` | 1–6 | 1, step 41 | Normal L1, unsuccessful L2 record; retain its outcome |
| `s1a4d7914034ddafc` | 1–3 | 0, step 19 | Unfinished discovery |
| `s1a4d7914034ddafc` | 2–3 | 0, step 17 | Unfinished discovery |
| `s1a4d7914034ddafc` | 2–6 | 0, step 17 | Unfinished discovery |

A companion accepted with different products does not supersede the original
homolysis channel. A finite monotonic scan, an exhausted optimizer, or SCF
nonconvergence does not prove barrierlessness or that another attempt would help.
Keep terminal nonfindings; resume the owning drivers to reconcile unfinished
work. Do not repeat the all-homolysis or thermochemistry-reset requests.

## Reverse abstraction: the missing tight route was found

CH3S's seven reverse searches are terminal nonfindings. CH3SS has 21
nonfindings and five accepted reverse searches. Two of the five are opposite
approaches to CH2S + SH <-> HCS + H2S:

- `49d0ff5991ad16cb180c/340660020000000000001_450570330000000000002_abstraction_4_2_1`
- `93359f497a7313176574/330330000000000000002_460700680000000000001_abstraction_1_5_3`

Both native L2 frequency logs match their exported geometries and have one
significant imaginary mode: -1157.1665 and -1157.3337 cm^-1. Across the two
geometries, proper atom-mapped RMSD is 0.0005237883 A; normalized absolute
mass-weighted imaginary-displacement overlap is 0.99999065. Electronic energies
differ by 8.90e-9 hartree. Together with the reversed IRC assignments, these
support counting this as one saddle, not two independent pathways. The logs'
SHA256 values and exact mapping are retained in `model_checks.json`.

The scout's E+ZPE saddle is 4.25176 kcal/mol above CH2S + SH and 2.87976 kcal/mol
above HCS + H2S. These are L2 scouting results, not replacements for a higher-level
paper PES. The CH2S + SH entrance complex is 1.30441 kcal/mol below its fragments
at L2. This supplies native evidence for the route; it does not by itself identify
every complex with the paper's geometry or authorize duplicate outer capture edges.

The other three accepted searches are:

| Reverse worker | Reaction suffix | Entrance characterization |
|---|---|---|
| `0995466154c2e06f9a55` | `abstraction_1_6_2` | L1 depth 0.47831 kcal/mol; below the 0.5 retention threshold |
| `26b8db7c87ded6112c4a` | `abstraction_3_1_2` | L1 depth 0.22149 kcal/mol; below threshold; single-well product |
| `b18dcbc8bc8c59de364b` | `abstraction_1_6_3` | L1-bound entrance, but L2 optimization exhausted 750 steps |

The last entrance's late iterations have approximately 0.0212 force maximum
and 0.0001 trust radius with repeated oscillation. Another identical step-budget
retry is not justified by this record. It is not an accepted L2 minimum.

**Remaining integration:** the reverse-search module currently exports accepted
results and blocks replacement of the combined model until integration. It does
not yet merge these entrance complexes, saddles, and product endpoints into the
ordinary PES. This is a real implementation gap. Do not set its manifest to
`integrated` by hand, invent shallow wells, connect two product asymptotes as if
they were a bound well, or present the old CH3SS MESS file as the updated network.

## HIR and the w7/w14/w15 question

The latest base L2 minima retain the geometries and native spectra previously
reviewed. Their proper permutation-aligned RMSDs are 0.7631 A (w7/w14), 0.7902 A
(w7/w15), and 1.0406 A (w14/w15). They are not identical optimized minima that
can simply be deduplicated by the cross-worker QC cache.

Each of w3, w7, w14, w15 now has two 12-point scans with normal native QC rows.
The raw electronic energy ranges below are relative to that worker's L2 minimum;
they are sampled values, not continuous saddle barriers or ZPE-corrected barriers.

| Well | Rotor 0 range, kcal/mol | Rotor 1 range, kcal/mol | Replayed current population check |
|---|---:|---:|---|
| w3 | 0.00016 to 4.47384 | 0.00012 to 5.61375 | No stereochemical exclusion in either scan |
| w7 | 0.00000 to 0.86920 | 0.00000 to 1.63794 | Both cross excluded stereoisomers |
| w14 | -1.21992 to 0.19244 | -1.18291 to 0.62812 | Both cross excluded stereoisomers |
| w15 | -0.00010 to 1.16743 | -0.00006 to 0.82642 | Both cross excluded stereoisomers |

For example, w14 rotor 1 point 08 is 0.0460 A from w7 after proper atom mapping;
w15 rotor 1 point 06 is 0.0527 A from w7. These constrained scan geometries are
nearby configurations, not newly verified unconstrained minima. The data supports
low-energy interconversion; it does not support an assumption that every
torsional route is prohibitively high. It also does not establish the global
lowest path or an independent-rotor representation of coupled motion.

`HIR.test_hir` was replayed against read-only native rows using the actual worker
input references. It excludes both full-domain rotors for w7/w14/w15 because
they cross separately represented configurations; w3 has no such exclusion.
Hessian projection and final model writing are separate checks. This avoids
counting those full rotor domains in addition to all three separate wells.
It leaves harmonic treatment for the excluded rotors. A merged multiwell
anharmonic state-density model is **not implemented or validated by this audit**;
the low-frequency partition-function question remains distinct from QC reuse.

## Reuse and a corrected diagnostic

Current ordinary/reverse worker receipts show 371 adopted shared calculations and
one waiting consumer. That consumer is `781540900660020000001_hir_1_11` in worker
`791671560700020000002-s893815ff03c7331d`; its producer is rotor 0 point 11 in the
same worker, token `134e0e4ae1d34904aa5d59374905abf5`. Producer job 52956302 has
exit code 0, a completed Sella log and an unconsumed `.pkl` handoff. Neither name
has a database row yet. Resume/adopt this evidence rather than rerun the point.
The uncaptured pool registry prevents an exact offline replay of that live claim.

Startup cleanup was also parsing `kinbot.log` itself as Gaussian because that
log quoted an earlier native SCF error, producing recursive misleading warnings.
The general cleanup now excludes controller/optimizer logs from QC parsing and
removes unconditional advice to resubmit. It does not change QC acceptance,
search outcomes, or retry settings. Focused validation: 218 tests passed, one
platform-specific skip, including controller-refresh races/idempotence, live
work/sampling preservation, shared-QC adoption and HIR population checks.
These are local simulated-site tests, not a verified new Blodgett refresh.
