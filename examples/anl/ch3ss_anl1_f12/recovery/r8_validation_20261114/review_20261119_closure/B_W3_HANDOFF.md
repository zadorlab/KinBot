# Complex B–W3 saddle check

One selected calculation tests the shallow maximum in the continuous complex-B
inward scan. Its seed is converged scan point 1, C–S = 2.712320 Å; the source
output SHA-256 and line are recorded in `w3_complex_b_saddle.json`.

The task runs an unconstrained TS optimization and numerical frequencies at
B2PLYP-D3(BJ)/cc-pVTZ, doublet, UltraFine, SCF(XQC,Conver=8). The native
Hessian-sign check and the result parser's single-imaginary-mode check are
retained. Geometry and reactive-mode inspection determine the subsequent
connectivity check.

Run on Blodgett in the existing persistent terminal session:

```bash
cd "$HOME/KinBot" || exit
test "$(git branch --show-current)" = composite || exit
git pull --ff-only origin composite || exit
KB_SPEC="examples/anl/ch3ss_anl1_f12/recovery/r8_validation_20261114/review_20261119_closure/w3_complex_b_saddle.json"
KB_CLOSURE="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114/w3_complex_b_saddle_20261119"
test -f "$KB_SPEC" || exit
if [ ! -e "$KB_CLOSURE" ]; then
  .venv/bin/python -B -m kinbot.anl.dispatch prepare "$KB_SPEC" "$KB_CLOSURE" || exit
fi
.venv/bin/python -B -m kinbot.anl.dispatch preflight "$KB_CLOSURE" || exit
.venv/bin/python -B -m kinbot.anl.dispatch drive "$KB_CLOSURE"
```

The owning dispatcher inherits the run's resource policy and excluded nodes,
retains attempt history, and adopts a saved/live task on repeated invocation.
The drive stays attached until the selected task reaches a terminal state.
An incomplete pre-existing stage is preserved for inspection.
The updated driver prints attachment, task/job changes, and a 60-second
heartbeat. Older running drivers remain silent until completion. In a second
terminal, inspect the cached status and the queue without restarting the owner:

```bash
cd "$HOME/KinBot" || exit
KB_CLOSURE="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114/w3_complex_b_saddle_20261119"
.venv/bin/python -B -m kinbot.anl.dispatch status "$KB_CLOSURE" --cached
cat "$KB_CLOSURE/state.json"
squeue --me -h -o '%i|%T|%R|%j|%Z'
```

Local simulated-site validation covers staging and theory/geometry preservation,
deferred admission, live-result adoption, repeated invocation, and rejection of
minimum or higher-order spectra. Native execution is performed by this handoff.
The task definition contains only `saddle`; further calculations depend on its
reviewed result.
