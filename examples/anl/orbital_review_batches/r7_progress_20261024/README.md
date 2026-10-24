# Retained production handoff, 24 October 2026

Run `bash examples/anl/tools/refresh_progress_20261024.sh PRODUCTION_ROOT`
from the updated composite checkout. The script validates and stages the
diagnostic inputs before replacing controllers. It then runs the existing
retained-run handoff with a cap of 22, retaining the three review slots already
configured in this run. Current bad-node states are added to existing exclusions.

Both PES/KB driver trees load an immutable source snapshot. Retained RotdPy
owners keep their sample databases and live jobs; their shared setup file routes
subsequent worker starts through that snapshot, including templates cached by old
owners. A running Python worker cannot load a replacement implementation midway
through a native calculation and is left to finish. The script retains the prior
C + SH3 rejection through the existing live reaction-decision interface.

The diagnostic campaigns live outside production in the sibling directory
`r7_progress_review_20261024`. One campaign owner runs them sequentially so their
combined admission uses at most the three reserved review slots:

1. Three new independent CASCI screens, each at 2, 3.5, 6 and 30 angstrom in
   aDZ-F12 and aTZ-F12: CH2(SH)S C–S cleavage, HCS2 + H2 exit, and CH2S + SH exit.
   The two exit screens use explicitly mapped, saved L2 product-side geometries
   from the accepted reaction records. Their H–H separation coordinates do not
   mean the parent well's H–H distance is stretched. These are broad orbital
   diagnostics; no compact active space or production model is inferred from
   them automatically.
2. Eight independent CASPT2-F12/aDZ-F12 relaxation probes: the CH2SHS C–S and
   S–S channels, each at the same four distances. The earlier verified native
   screens and explicit compact orbital selections are reused. These probes
   start from the actual requests, not the earlier provisional coordinates.
   Only the separation is constrained; internal fragment relaxation remains
   free. The report lists all intrafragment pair distances and flags stretched
   reference bonds. Final structures, native orbital character and state
   continuity still need review before either model is published.

The new three screens need compact-space selection and subsequent relaxation
checks after their native evidence is available. The eight existing-selection
probes do not substitute for those later checks. Completed provisional screens
do not certify relaxed production geometries.

Inspect `orbital_evidence.txt` and `relaxations/report.json`. No diagnostic result
is automatically promoted into the production inputs. The fresh CC iteration
and sampling solver retries are implemented and regression-tested; their success
on the failed production calculations must be established by the remote native
results. See `docs/LIVE_DRIVER_CONTROLS.md` for retry bounds and live reaction
decisions. In particular, the MRCC higher-order engine's compiled iteration
ceiling is not overridden by an ineffective `ccmaxit` edit.

Repeating this exact handoff resumes its saved diagnostic campaigns and avoids
another controller replacement once `refresh_complete.json` exists. A different
checkout revision, changed scientific evidence, or changed plan requires
inspection rather than silently starting a duplicate campaign. The diagnostic
owner holds an exclusive lock; do not start a separate owner for the same
review slots while this script is running.

## CH2(SH)S C-S follow-up

The 24 October native character extraction is retained in
`ch2_sh_s_cs_character.txt`. Atom 2 is the departing SH sulfur; atom 3 remains
with carbon in CH2S. At 3.5-30 angstrom, orbitals 20/21 lie primarily on atom 2
and 19/22 form the CH2S pi pair. At 2 angstrom the 19/22 pair instead contains
substantial breaking C(1)-S(2) bond character and 20/21 substantial S(3)
character. Both bases show this reorganization. This supports a comparison,
not a claim of demonstrated subspace/root continuity. The extra 18/23 pair is
the retained C(1)-S(3) sigma pair, not the S-H pair from the earlier CH2SH case.
Printed AO coefficients are not atomic populations, and their signs can change
with orbital phase; compare chemical subspaces, not signs or sorted occupations
alone.

`ch2_sh_s_cs_candidates.json` defines 24 independent fixed-geometry jobs:
CAS(5,4) from source 19-22 and CAS(7,6) from 18-23, in both bases at 2, 2.75,
3.5, 6, 30 and 42 angstrom. Each optimizes the same two equally weighted
doublet states and computes both CASPT2-F12 roots, retaining shift 0.1 and IPEA
0.25. Basis-matched, geometry-specific CASCI seeds are reused; the new 2.75/42
points use the 3.5/30 seeds. Compare energies relative to each model's own
asymptote and inspect the two roots, rather than choosing by absolute energy.
The 42 point is a finite long-range check, not automatic acceptance of the
previous coverage flag or a changed production asymptote.

`ch2_sh_s_cs_relaxations.json` stages ten independent aDZ-F12 CASPT2-F12
optimizations: both spaces at 2, 2.75, 3.5, 6 and 30 angstrom. Only C(1)-S(2)
separation is fixed. The actual request geometry and verified 2 angstrom seeds
are used; internal fragment changes remain observable. Each has a four-hour
limit and 50 optimizer steps. A failed calculation remains review evidence;
it is not accepted as a relaxed geometry. The `probe5_4` and `probe7_6` files
approve only diagnostic seeds (`diagnostic_only: true`); the production
publisher explicitly refuses these files.

After pulling the follow-up, run the following with the existing environment:

```
.venv/bin/python -B -u examples/anl/tools/continue_ch2_sh_s_review.py PRODUCTION_ROOT
```

Use `nohup` for a background owner. It waits on the original review handoff
lock, adopts any surviving allocations through their existing owners, and
finishes the original screens/eight probes before the new comparison. The
same three reserved slots are used sequentially. No production controller
replacement, sample cancellation, new production model or additional review
capacity is requested. It reuses the verified core-code snapshot from
`refresh_complete.json`; the new code here is confined to example tools.
Repeated invocation retains the staged graphs, completed jobs and live job
identities; a second simultaneous owner is refused. Recipe changes during
the wait are refused before new submissions.

If the fixed-geometry comparisons fail, all their existing allocations are
settled and the new relaxed probes remain staged for review. Otherwise the
ten relaxed probes run automatically. Inspect
`ch2_sh_s_followup/orbital_evidence.txt`, `relaxations/report.json` and
`followup_status.json` under the review directory. Native success, orbital
continuity, relaxed connectivity and the eventual production model still need
review; no diagnostic completion is automatically promoted.
