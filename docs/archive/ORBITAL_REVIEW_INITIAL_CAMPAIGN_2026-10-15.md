# Historical initial-campaign instructions

These commands describe the initial October 15 campaign, not the current
confirmation procedure. Preserve them as provenance; do not restart this old
panel to repeat completed work.

## Blodgett: stop, update, reuse and restart

Run this block in the existing Blodgett SSH session. `set -e` is confined to a
child shell, so a failed command does not exit the SSH login shell. It stops
only the specified campaign trees. It does not remove old results or WFU files.
Pull first to install the corrected cleanup command. Cancellation waits up to
five minutes by default: Slurm can retain a job in `COMPLETING` while its
processes and node epilog finish. The command prints progress and saves
`cleanup_status.json`; a timeout reports the remaining jobs/processes, returns
nonzero and leaves the directory in place. A repeated `--stop` is safe once
the queue has drained. Use a longer `--timeout SECONDS` if needed; do not
restart while scoped jobs remain.
The fixed new campaign path makes a repeated launch resume the same work;
use a new path if you edit the comparison plan.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
OLD=/home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525
NEW=/home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_20261015
PLANS=/home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_20261015_inputs

env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
"$KB_PY" -m kinbot.run_cleanup "$OLD" --stop --timeout 300
if [ -f "$NEW/campaign.json" ]; then
    "$KB_PY" -m kinbot.run_cleanup "$NEW" --stop --timeout 300
fi
"$KB_PY" -m pip install -e . --no-build-isolation --no-deps
"$KB_PY" -m kinbot.orbital_review --help
"$KB_PY" -m kinbot.orbital_campaign archive-excluded "$OLD"
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
mkdir -p /home/hram/.cache/kinbot/molpro

mkdir -p "$PLANS"
if [ ! -f "$NEW/campaign.json" ]; then
    STAMP=$(date +%Y%m%d_%H%M%S)
    "$KB_PY" -m kinbot.orbital_review plan \
        /home/hram/KinBot/examples/anl/ch3s_anl1_f12/orbital_reconnaissance.json \
        "$PLANS/ch3s_$STAMP.json"
    "$KB_PY" -m kinbot.orbital_review plan \
        /home/hram/KinBot/examples/anl/ch3ss_anl1_f12/orbital_reconnaissance.json \
        "$PLANS/ch3ss_$STAMP.json"
    "$KB_PY" -m kinbot.orbital_review prepare "$PLANS/ch3s_$STAMP.json" "$NEW" \
        --additional-plan "$PLANS/ch3ss_$STAMP.json" --reuse-from "$OLD"
else
    "$KB_PY" -m kinbot.orbital_campaign retry-failed "$NEW" --include-interrupted
fi
"$KB_PY" -m kinbot.orbital_campaign preflight "$NEW"
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$NEW" \
    --max-jobs 20 --max-refinements 0 \
    >> "$NEW/controller.log" 2>&1 < /dev/null &
printf '%s\n' "$!" > "$NEW/controller.pid"
"$KB_PY" -m kinbot.orbital_campaign progress "$NEW"
tail -n 30 "$NEW/controller.log"
BASH
```

This first run resolves initialization and optimized reference-orbital
behavior. It does not assume CAS(5,4), CAS(9,6), a particular inner/outer
boundary or a production average. We can use its export with the existing
CASCI archive to choose the smaller candidate comparisons hands-on. Completed
new jobs are reused on resume; an interrupted allocation is explicitly reset.
A native failure still needs its reported diagnostic addressed before retry.
With the current examples, preparation produces 72 finite-panel jobs and
nine native pilots, under one twenty-job cap. Every requested basis is retained;
the examples currently request only TZ-F12. No correlated checks are enabled
in this initial reference-orbital campaign.

Read status at any time:

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
/home/hram/KinBot/.venv/bin/python -m kinbot.orbital_campaign progress \
    /home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_20261015
tail -n 50 /home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_20261015/controller.log
squeue -u hram -o '%.18i %.16P %.48j %.10T %.10M %.6D %R'
BASH
```

After the controller finishes, export:

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
NEW=/home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_20261015
STAMP=$(date +%Y%m%d_%H%M%S)
"$KB_PY" -m kinbot.orbital_review report "$NEW" "$NEW/review_report_$STAMP.json"
"$KB_PY" -m kinbot.orbital_campaign collect "$NEW" \
    "/home/hram/kinbot-runs/ch3s_ch3ss_orbital_review_$STAMP.tar.gz"
BASH
```

Copy the printed archive to the local checkout for review. Keep both campaign
directories on Blodgett because the selection and restart provenance refer
to their WFU paths. Production launch commands follow the user's actual
space/surface selections; the current example inputs deliberately reference
state/surface model files that have not yet been scientifically selected.
