# Network revisit and live review handoff

For the subsequent numerical-failure repair and cross-well stable/TS reuse, use
[the October 21 recovery handoff](RECOVERY_AND_PES_REUSE_2026-10-21.md). It preserves
ongoing exploration and does not repeat the blanket rejected-channel revisit.

## What this recovery changes

The recovery is implemented; the command below applies it on the cluster.
It preserves the run directory, native
artifacts, databases and scheduler IDs. It changes `vdW_detection` to zero,
permits fresh L1 search attempts, rechecks accepted paths and gives each recorded
rejected reaction one additional attempt. Valid completed work is parsed/reused;
unchanged running and pending allocations retain their IDs and original runtime
settings. Exploration completes across each PES before new downstream
submissions resume. Existing allocations continue during exploration.

Verified duplicate TS graphs are the only selected cancellations. Reuse needs
matching endpoints, ordered atoms/state, L2 method, proper-rotation geometry,
energy/ZPE and frequencies with one imaginary mode. Same products alone never
establish a duplicate. The canonical graph favors completed/active work and is
then frozen. Duplicate graphs and native results are archived, not deleted.
Unsupported atom permutations/reflections are not guessed. Sharing is confined
to one network, never between the two independent runs.

The archived earlier inventory contained repeated CH3S migration/addition/shift
paths reaching essentially the same saddle, alongside another distinct saddle
to the same product. CH3SS also had two paths converging to one saddle. These are
evidence for selective reuse, not authorization to collapse every reaction with
the same products. The refresh recomputes its selection from current manifests
after quiescing controllers; it does not hardcode those names or old Slurm IDs.

The paper's vdW complex is reported about 4.77 kcal/mol below its separated
products, so the old 0.5 cutoff alone does not establish why it was absent.
Exploration ordering, failed guesses, geometry/method differences and complex
handling all matter. Revisit does not guarantee reproducing every published
stationary point. See Veillet et al., DOI
[10.1039/D4CP04873C](https://pubs.rsc.org/en/content/articlehtml/2026/cp/d4cp04873c).

## Cluster command for the existing two-PES run

These paths identify the user's existing run; the implementation contains no
cluster/account-specific policy. Run on the controller host after updating:

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
git pull --ff-only origin composite
KB_RUN=/home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage
KB_PY=/home/hram/KinBot/.venv/bin/python
"$KB_PY" -m kinbot.run_control enable "$KB_RUN"
"$KB_PY" -u -m kinbot.run_refresh "$KB_RUN" \
  --vdw-detection 0 --revisit-rejected --deduplicate-ts --apply
BASH
```

Do not add pending replacement, VRC basis migration or blanket cancellation.
The existing exclusive-node policy, shared node pool/reserve/exclusions,
node-local scratch selection, CPU/memory settings, walltime/retry controls,
MR native output/audits, aDZ-F12 relaxed scans, reviewed active spaces,
extrapolated correction settings, electronic-state corrections and 1% ROTD
flux target are preserved. Existing scientific inputs are compared, not rebuilt
from the repository's example defaults.

The command reports its archive, rejected-channel count, duplicate allocation
selection and retained Slurm IDs. Unknown scheduler ownership or an in-flight
submission is a stop condition, not permission to duplicate it. If interrupted,
repeat the same command with the same code/options to reconcile the journal.
An already launched identical handoff is a no-op.

```bash
watch -n 10 /home/hram/KinBot/.venv/bin/python -m kinbot.run_monitor \
  /home/hram/kinbot-runs/ch3s_ch3ss_caspt2_production_20261018_r6_salvage
```

Expected stages: exploration/revisits -> all eligible wells ready -> concurrent
ANL and reviewed VRC -> correction assembly -> PES ROTD sampling -> MESS. New
unreviewed channels show `REVIEW HOLD` and retain a live owner. Provide their
reviewed models in `ch3s/input.json` or `ch3ss/input.json`, then use
`kinbot.run_control resubmit INPUT --apply`; inspect `--status` for application.
No further controller restart is needed for supported live changes.

## Verification scope and remaining limits

Regression coverage includes real subprocess exploration ordering, worker
quiescence, late-well discovery, retained running/pending jobs, idempotent
recovery, TS identity rules, append-only CBH plans, complex formation-energy
accounting, input publication/partial-publication recovery, stale-owner
rejection, and reviewed-channel progress during another channel's review hold.
Native cluster calculations and scheduler actions cannot be executed from the
development checkout; inspect the launch report and next monitor on Blodgett.

A successful handoff restores the intended lifecycle and reuse rules. It cannot
guarantee native convergence, discover every physical reaction, approve an
active space, recover lost scratch amplitudes, or prevent a future filesystem
outage. A change to already accepted chemistry/geometry is deliberately refused
unless its dependent work is explicitly migrated. Closed conformer searches
and a PES already in final sampling/solver aggregation do not have a live
network-edit window. Existing retry policies apply; review holds require input.

Final local verification: the focused regression selection passed **121 tests**
with one skipped. The full suite passed **2,426 tests and 429 subtests**, with
seven skipped and the same **11 baseline failures** reproduced before this
change. No new full-suite failures were observed.

The full suite has known baseline failures in older CFOUR, sulfur-example,
CASPT2 handoff, MRCI probe, orbital-campaign and HOMO fixture expectations.
They were reproduced before this change; they must not be described as a clean
full-suite pass. Native Gaussian/Molpro/MRCC/CFOUR/Slurm execution was not part of
this local verification.
