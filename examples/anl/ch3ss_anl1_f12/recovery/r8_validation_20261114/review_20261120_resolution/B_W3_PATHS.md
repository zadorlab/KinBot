# Reuse the completed B–W3 saddle for its two paths

Job **53245807** converged at B2PLYP-D3BJ/cc-pVTZ, charge 0,
multiplicity 2. Its spectrum has one imaginary mode, **154.989i cm⁻¹**,
and eleven positive modes. The C–S approach distance is 2.7173 Å.
The stationary point is established; its two endpoints remain to be assigned.

This handoff imports that completed geometry, energy and native modes, then
runs the existing Sella IRC owner in both directions at the same electronic
level. It stops after the paths and endpoint optimizations. Existing endpoint
spectra can then be matched and reused. It requests two allocations, each with
eight cores, and inherits the shared R8 queue cap and excluded nodes.

```bash
cd "$HOME/KinBot" || exit
test "$(git branch --show-current)" = composite || exit
git pull --ff-only origin composite || exit

KB_VALIDATION="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114"
KB_SOURCE="$KB_VALIDATION/w3_complex_b_saddle_20261119"
KB_PATHS="$KB_VALIDATION/w3_complex_b_paths_20261120"
KB_STAGE="examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/review_20261120_resolution/stage_b_w3_paths.py"

.venv/bin/python -B "$KB_STAGE" "$KB_SOURCE" "$KB_PATHS" || exit
.venv/bin/python -B -m kinbot.saddle_recheck run "$KB_PATHS" --paths-only
```

Repeated invocation adopts completed/running work under the same owner lock.
The source saddle directory remains unchanged. Both directions are bounded at
300 IRC steps with dx 0.05, using the successful earlier direct-L2 strategy.
An unconverged trajectory remains explicitly unconverged even if its later
endpoint optimization succeeds.

Read-only progress from another terminal:

```bash
KB_PATHS="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114/w3_complex_b_paths_20261120"
cat "$KB_PATHS/results.json"
squeue --me -h -o '%i|%T|%R|%j|%Z'
```

Local tests exercise real native Gaussian parsing, imported-saddle reuse,
both path submissions, live-work adoption, completed-work reuse and source
preservation with simulated scheduling. The native path results will come from
Blodgett; the local tests do not establish their endpoints.
