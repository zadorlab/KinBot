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
