# R8 L1/L2 scouting: CH3S and CH3SS

These inputs use the ordinary `python -m kinbot.pes` driver: full allowed
reaction search, L2 refinement, and automatic PES-wide MESS input assembly.
The paired CH3SS input is at `../../ch3ss_anl1_f12/scouting_r8/input.json`.

| Setting | R8 |
| --- | --- |
| L1 | B3LYP-D3(BJ)/cc-pVDZ, Gaussian/Sella and native Hessian |
| L2 | B2PLYP-D3(BJ)/cc-pVTZ, Gaussian/Sella and native Hessian |
| Conformers | Existing search, one selected conformer for kinetics |
| Rotors | Standard 12-point hindered-rotor scans; no multi-conformer TST |
| Search ceiling | L1 150 kcal/mol relative to initial well; L2 ceiling 150 |
| vdW selection | Matched L1 E+ZPE check, depth 0.5 kcal/mol |
| Barrierless kinetics | Existing PST, coefficient 0.6 au, exponent 6 |
| MESS | `me=1`, `run_me=0`: assemble automatically, analyze afterward |
| Scheduling | Eight processors, exclusive nodes, shared cap of 20 QC jobs |
| Native checkpoint copies | Every 600 s; two generations, 16 GiB ceiling |

The existing MESS workflow requires `rotor_scan=1`; the standard 12-point
scans are retained. This avoids adding a new RRHO-only execution mode.

The historical preset name `b3lyp-b2plyp-anl` supplies only the two Gaussian
profiles. `composite_method=""` and `L3_calc=0` disable composite/single-point
work. `rotdpy_run=0`, empty VRC selectors, and `rotdpy_mess_mode="off"` select
PST without a ROTD_py result or CAS-model gate. The older driver may create
unused `molpro/` and `rotdPy/` bookkeeping directories; it submits no such jobs.

PST is a scouting assumption, not a fitted capture rate. After MESS runs,
inspect its native high-pressure association coefficients (cm^3/s), including
temperature dependence, before pruning by branching. The collision parameters
(He, reactant epsilon=300 K and sigma=4.5 angstrom) retain the earlier example
assumptions; they are not capture-rate parameters. The default PST HPL has not
been quantitatively validated for these channels. MESS's microcanonical output
also supports comparison at photon energy without assuming thermal photolysis.

## Temperature grid

The inputs specify 300, 400, 500, 600, 800, 1000, 1200, 1500, 1800 and 2000 K.
Pressures span 0.01–760 torr. The standard model energy limit is 400 kcal/mol.
The 193 nm photon energy (148.142152 kcal/mol) is a separate microcanonical
comparison, not a temperature. MESS's existing 0–200 kcal/mol microcanonical
output grid brackets it at 148 and 149 kcal/mol. Thermal rates and photon-energy branching answer
different scouting questions and should be labeled separately.

## Stage on Blodgett: no submissions

Use the updated `composite` checkout. This creates a fresh directory and fails
if that name exists. R7 remains untouched. These inputs do not import its QC
results or terminal search history; general cross-run import remains a separate
R9 readiness item. The old high-level inputs and review evidence are preserved
for R9, after manual channel selection from this scouting pass.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8"
mkdir "$KB_RUN"
mkdir "$KB_RUN/ch3s" "$KB_RUN/ch3ss"
cp examples/anl/ch3s_anl1_f12/scouting_r8/input.json "$KB_RUN/ch3s/input.json"
cp examples/anl/ch3ss_anl1_f12/scouting_r8/input.json "$KB_RUN/ch3ss/input.json"
"$PWD/.venv/bin/python" -B -m kinbot.run_policy set "$KB_RUN" \
  --max-jobs 20 --timeout-retries 1 --optimization-step-retries 2 \
  --exclude node1,node9,node10,node13,node14,node16,node20,node21,node22,node23,node24,node29,node30,node33,node34
"$PWD/.venv/bin/python" -B - "$KB_RUN" <<'PY'
import sys
from pathlib import Path
from kinbot.parameters import Parameters
from kinbot.theory import automatic_cbh_anl_requested
for case in ('ch3s', 'ch3ss'):
    path = Path(sys.argv[1])/case/'input.json'
    p = Parameters(str(path)).par
    assert not automatic_cbh_anl_requested(p) and not p['L3_calc']
    assert not p['rotdpy_run'] and not p['vrc_tst_scan'] and not p['vrc_tst_noscan']
    assert p['rotdpy_mess_mode'] == 'off' and p['me'] == 1 and p['run_me'] == 0
    assert p['high_level'] == 1 and p['rotor_scan'] == 1
    assert p['queue_job_limit'] == 20 and p['queue_exclusive'] == 1 and p['ppn'] == 8
    print('STAGED', case, path, 'L1/L2 -> MESS/PST; no chemistry submitted')
PY
BASH
```

Exclusions retain the known site exclusions from the preceding run. Change
them only with site evidence. The cap is not a reservation: still-running
orbital diagnostics may occupy nodes independently; Slurm controls availability.

## Launch when ready

Use the existing Gaussian-enabled Blodgett environment. MESS is needed only
when the generated inputs are executed later. No ANL or rotdPy preflight runs.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
KB_PY="$PWD/.venv/bin/python"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8"
for KB_CASE in ch3s ch3ss; do
  test -f "$KB_RUN/$KB_CASE/input.json"
  test ! -e "$KB_RUN/$KB_CASE/driver.log"
done
for KB_CASE in ch3s ch3ss; do
  (
    cd "$KB_RUN/$KB_CASE"
    nohup "$KB_PY" -B -u -m kinbot.pes input.json > driver.log 2>&1 </dev/null &
    echo "$!" > launch.pid
    printf '%s controller %s; log %s/driver.log\n' "$KB_CASE" "$!" "$PWD"
  )
done
BASH
```

Run the launch block once. For interruption, inspect the controller/job
records and resume through the normal owning driver; do not delete logs to
bypass the fresh-launch check.

Deliverables are each case's combined `pesviewer.inp` and `me/*.inp`, worker
databases, accepted networks and search histories. Energetic/depth selection
runs during exploration; 1%/5% kinetic trimming happens afterward. This cannot
guarantee a smaller discovery tree before rate analysis. Retain all results for
choosing R9's explicit production channels. Standalone orbital diagnostics can
finish independently and be reviewed only for the channels selected for R9.
