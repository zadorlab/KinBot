# Sulfur TZ-F12 restart and review — 14 October 2026

This is the current procedure, superseding the multi-basis restart blocks in
`Blodgett_workflow_runbook.md`. The user requested stopping the broad recon,
using augmented TZ-F12 for VRC production, and removing routine VRC CBS work.
No cluster action was performed from the development machine.

## What changes

- New sulfur recon configurations use only Molpro's native `avtz-f12`.
  The full grids retain all distances, orientations, spins and roots:
  **31 CH3S + 75 CH3SS = 106 full probes**, plus **1 + 2 = 3 pilots**.
  Adaptive refinement can add probes. This is one sixth of the former
  six-basis initial job count; wall-time reduction is not exactly proportional.
- The `anl1-f12` VRC preset uses `avtz-f12` for relaxed CASPT2-F12b paths,
  sampling, the ROHF-UCCSD(T)-F12b high-spin anchor, both CASPT2-F12b and
  MRCI-F12b+Q gaps, and direct MRCI-F12b+Q high-spin leaves. The two gap
  estimators are still averaged. Routine VRC calculations do not extrapolate.
- Explicit user CBS overrides and historical defaults remain supported.
  The **CBH/ANL1-F12 thermochemistry recipe is unchanged**, including its
  reference CBS and additive corrections. Changing a VRC basis is not a
  change to a formation-enthalpy recipe.
- Existing campaign snapshots remain immutable. `exclude-bases` changes
  selection in their manifests; retained completed TZ-F12 jobs are reused.
  `archive-excluded` moves unused basis directories under `archived_excluded`
  in that same campaign, records their location, and keeps all evidence.
  It does not reclaim disk space or delete shared scratch.

## Resume after the legacy CH3S pilot rejection

The first TZ-only restart exposed a compatibility regression in `ba4beef`.
The completed CH3S pilot `avt_f12__cs_r_1p599154942` used the original Molden
worker: its command ran the prepared `.inp`, and its execution record never
contained `native_input` or `native_input_sha256`. The user verified that the
manifest, execution record and actual prepared input all have SHA-256
`8d2af27a5d011a5f26da685592324f256a4a34b3b0fe25b419dbe2ab67f97035`.
Requiring a later-format `.native.inp` incorrectly rejected this valid result.

The compatibility fix validates that older command/input pair, then applies
the existing termination, Molden, occupations, roots and WFU checks. It never
treats a partial modern record as legacy or replaces a mismatched hash. Valid
completion clears the old publication timeout. The regression test reproduces
the legacy record and verifies that restarting neither submits nor rewrites
the completed pilot. The affected recon/campaign suite passes 60 tests.

The campaign is already TZ-only and its excluded directories are archived.
Use this shorter block after pulling the compatibility fix. It verifies all
three retained pilots before launching; it does not repeat cleanup of artifacts
or rerun successful calculations. Failures stay inside the child Bash process.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
CAMPAIGN=/home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525
test -f "$CAMPAIGN/campaign.json"
"$KB_PY" -m kinbot.run_cleanup "$CAMPAIGN" --stop
test "$(git branch --show-current)" = composite
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
"$KB_PY" -m pip install -e . --no-build-isolation --no-deps

"$KB_PY" - <<'PY'
import json
from pathlib import Path
from kinbot.molpro_reconnaissance import audit
root = Path('/home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525')
for system in json.loads((root / 'campaign.json').read_text())['systems']:
    for key in ('pilot_manifest', 'manifest'):
        report = audit(system[key])
        assert {j['basis'] for j in report['jobs']} == {'avtz-f12'}
        completed = sum(j['status'] == 'complete' for j in report['jobs'])
        print(system['title'], key, completed, '/', len(report['jobs']), 'complete')
        if key == 'pilot_manifest' and completed != len(report['jobs']):
            failures = {j['name']: j.get('validation_errors', [])
                        for j in report['jobs'] if j['status'] != 'complete'}
            raise SystemExit('Pilot validation failed; controller remains stopped: ' + json.dumps(failures))
PY

export KINBOT_SLURM_IO_RETRIES=2
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
export KINBOT_MOLPRO_SCRATCH_FALLBACKS=/scratch2/hram/kinbot-molpro
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
"$KB_PY" -m kinbot.orbital_campaign preflight "$CAMPAIGN"
cd /home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525
if [ -f controller.log ]; then
    mv controller.log "controller.previous.$(date +%Y%m%d_%H%M%S).log"
fi
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$CAMPAIGN" \
    --max-jobs 10 --max-refinements 3 \
    > controller.log 2>&1 < /dev/null &
printf '%s\n' "$!" > controller.pid
printf 'Campaign: %s\nController PID: %s\n' "$CAMPAIGN" "$!"
BASH
```

## Initial stop, update, archive unused bases, and restart on Blodgett

Paste this **whole block on Blodgett**, after the changes have been published
to `origin/composite`. Every path is explicit. `set -e` and `exit` run inside
a child Bash process, so failure returns to the SSH login shell. A failure
after cleanup leaves the campaign stopped, rather than launching old code.
The cleanup checks live WorkDir/PID ownership and protects its invoking shell
and ancestors. It uses no user-wide `pkill`, `killall`, or unscoped `scancel`.

The campaign path below is the recorded current campaign. The 22:46 snapshot
had CH3S pilots 6/6, CH3SS pilots 12/12, six completed CH3S full probes and an
empty queue. Job 52907323 completed natively but its WFU name was silently
truncated to 32 characters. Commands discover current ownership; they do not
assume old job IDs are still live.

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
KB_REPO=/home/hram/KinBot
KB_PY=/home/hram/KinBot/.venv/bin/python
CAMPAIGN=/home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525
test -f "$CAMPAIGN/campaign.json"
test -x "$KB_PY"

"$KB_PY" -m kinbot.run_cleanup "$CAMPAIGN" --stop

test "$(git branch --show-current)" = composite
env -u LD_LIBRARY_PATH -u LD_PRELOAD \
    git -c submodule.recurse=false pull --ff-only origin composite
"$KB_PY" -m pip install -e . --no-build-isolation --no-deps

"$KB_PY" - <<'PY'
import json
from pathlib import Path
from kinbot.orbital_campaign import archive_excluded
from kinbot.molpro import molpro_wfu_name
from kinbot.molpro_reconnaissance import audit
from kinbot.vrc_correction import _VRC_ACCURACY_PRESETS
preset = _VRC_ACCURACY_PRESETS['anl1-f12']
assert preset['vrc_tst_high_basis'] == 'avtz-f12'
spin = preset['vrc_tst_spin_splitting']
assert spin['reference']['basis'] == 'avtz-f12'
assert all(level['basis'] == 'avtz-f12' for level in spin['gap_levels'])
assert len(molpro_wfu_name('avd__cs_orient_f1_bond_a60_r1p776838825')) <= 31
for name in ('ch3s', 'ch3ss'):
    path = Path('/home/hram/KinBot/examples/anl') / f'{name}_anl1_f12/orbital_reconnaissance.json'
    assert [b['basis'] for b in json.loads(path.read_text())['bases']] == ['avtz-f12']
print('TZ-F12 source, checkpoint fix, and archive support verified')
campaign = Path('/home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525/campaign.json')
for system in json.loads(campaign.read_text())['systems']:
    for key in ('pilot_manifest', 'manifest'):
        report = audit(system[key])
        print(system['title'], key, sum(j['status'] == 'complete' for j in report['jobs']),
              'complete of', len(report['jobs']))
PY

"$KB_PY" -m kinbot.orbital_campaign exclude-bases "$CAMPAIGN" \
    --basis aug-cc-pVDZ --basis aug-cc-pVTZ --basis aug-cc-pVQZ \
    --basis aug-cc-pV5Z --basis avdz-f12 --basis avqz-f12 --basis av5z-f12
"$KB_PY" -m kinbot.orbital_campaign archive-excluded "$CAMPAIGN"
"$KB_PY" -m kinbot.orbital_campaign retry-failed "$CAMPAIGN" --include-interrupted

export KINBOT_SLURM_IO_RETRIES=2
export KINBOT_MOLPRO_SCRATCH=/home/hram/.cache/kinbot/molpro
export KINBOT_MOLPRO_SCRATCH_FALLBACKS=/scratch2/hram/kinbot-molpro
export KINBOT_MOLPRO_SCRATCH_MIN_MB=65536
"$KB_PY" -m kinbot.orbital_campaign preflight "$CAMPAIGN"
"$KB_PY" -m kinbot.orbital_campaign progress "$CAMPAIGN"
printf '%s\n' "$CAMPAIGN" > /home/hram/kinbot-runs/last_sulfur_orbital_campaign.txt

cd /home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525
if [ -f controller.log ]; then
    mv controller.log "controller.previous.$(date +%Y%m%d_%H%M%S).log"
fi
nohup "$KB_PY" -u -m kinbot.orbital_campaign run "$CAMPAIGN" \
    --max-jobs 10 --max-refinements 3 \
    > controller.log 2>&1 < /dev/null &
printf '%s\n' "$!" > controller.pid
printf 'Campaign: %s\nController PID: %s\n' "$CAMPAIGN" "$!"
BASH
```

If the source assertion or pull fails, do not bypass it. Local cluster source
edits are preserved; this block does not reset or stash them. If a submission
is recorded as `submitting`, reconcile it with Slurm before retrying: it may
have been accepted even though the controller never received its job ID.
`--include-interrupted` only resets known, no-longer-live allocations with no
terminal record or a `running` record. Successful outputs are preserved, and
parser-only rejections of native completion are not blindly recalculated.
The read-only native-evidence audit above recovers unambiguous truncated WFU
names before exclusion, preserving successful calculations. The supplied job's
runtime showed `/scratch2/hram` permission denied and a successful home-cache
fallback; the restart therefore tries that known working location first.

Monitor in a separate shell (this does not submit jobs):

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
/home/hram/KinBot/.venv/bin/python -m kinbot.orbital_campaign progress \
    /home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525
tail -n 60 /home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525/controller.log
squeue -u hram -o '%.18i %.16P %.48j %.10T %.10M %.6D %R'
BASH
```

The original example grids have serial checkpoint dependencies along each
basis/channel; ten is a maximum allocation count, not a promise of ten
simultaneously runnable probes with a single basis. Orientation probes remain
independent. Do not add controllers to fill idle slots.

When the controller has stopped after review-required completion, collect:

```bash
bash <<'BASH'
set -euo pipefail
cd /home/hram/KinBot
/home/hram/KinBot/.venv/bin/python -m kinbot.orbital_campaign collect \
    /home/hram/kinbot-runs/ch3s_ch3ss_orbital_20261014_193525 \
    /home/hram/kinbot-runs/ch3s_ch3ss_tzf12_review.tar.gz
BASH
```

The archive includes excluded-basis text/export evidence at its relocated
paths. Large WFU files remain on Blodgett. Do not delete `archived_excluded`
until the old diagnostics and any useful basis comparisons have been reviewed.

## Active spaces and remaining accuracy questions

The current reconnaissance uses **CAS(13,8) for CH3S and CAS(19,12) for CH3SS**
as diagnostic heavy-atom-valence spaces. CH3SS's production input instead has
a provisional **CAS(9,7)**. Neither count is an approved orbital selection.
The supplied CH3S DZ orientation probe has four fractional occupations
1.51137, 1.51137, 1.49450, 0.51729, and four other valence occupations
1.99808, 1.99228, 1.98755, 1.98755. That supports testing a four-orbital space,
but this single geometry/basis supplies neither orbital identities nor
evidence that the nearly doubly occupied orbitals remain inert on dissociation.

| Channel | Spin-free product terms | Diagnostic roots in each spin manifold |
| --- | --- | --- |
| CH3S C–S → CH3 + S | doublet radical × S(3P) | 3 doublets and 3 quartets |
| CH3SS C–S → CH3 + S2 | doublet radical × S2(3Σ) | 1 doublet and 1 quartet |
| CH3SS S–S → CH3S + S | CH3S(2E) × S(3P) | 6 doublets and 6 quartets |

These counts follow the declared spin-free terms, not a decision to sum every
root into the inner reactive flux. Their full spin-free microstate totals are
18, 6 and 36. Resolved SO/vibronic levels replace the appropriate spin counting;
they are not extra factors to multiply into it.

Review each channel separately. Keep the breaking-bond bonding/antibonding
pair, fragment radical orbitals, and complete near-degenerate sulfur p/CH3S
E manifolds. Compare their characters and CI contributions at the bonded
geometry, interaction region, all orientation probes, and separated limit.
For CH3S, CAS(5,4), spanning the methyl radical orbital plus S 3p orbitals,
is a physically motivated **candidate to test** against CAS(13,8). For CH3SS,
test CAS(9,7) against the larger probe space independently for C–S and S–S;
do not transfer the C–S choice to S–S without evidence.

Near-0/2 natural occupations suggest candidates for removal, but state averaging
can conceal state-specific correlation. Preserve whole coupled manifolds,
inspect root character/overlaps near crossings, and compare reduced/larger
spaces using **correlated interaction energies and spin gaps**, after each
curve's asymptotic zero is removed. Add hydrogen antibonding projections as a
sensitivity check if C–H relaxation enters; the heavy-atom probe omits independent
H-centered virtual targets. `START/DIAB` aids orbital continuity; it does not
certify electronic-root identity. Molpro documents the CASCI/occupation review
and separate state-specific natural orbitals in its
[AVAS](https://www.molpro.net/manual/doku.php?id=automated_construction_of_atomic_valence_active_spaces)
and [MULTI](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi) manuals.

Before using `molpro_reconnaissance select`, reorder the **actual WFU orbitals**
into the intended inactive/active/virtual blocks. Merely changing electron and
orbital counts does not select orbital identities. Save reviewed TZ-F12 records
per channel (and per bottleneck if needed), then use the surface model's
stage `parameters.vrc_tst_state_model_file` override. Keep geometry, orbital
characters, root assignments and review reasoning with the record.

The [NH2O paper, section 2.2](https://backend.orbit.dtu.dk/ws/portalfiles/portal/337661919/nh2_o_revised_final.pdf)
uses different inner/outer state averages and sometimes different active spaces.
It treats all outer electronic contributions under an explicit rapid-ISC
assumption and selects attractive inner contributions. Those assumptions and
counts are not automatically transferable to sulfur. Its SI counting and Ne
convolution verifier passes here; this is not a reproduction of its full rates.

TZ-F12 spin splitting is consistent with the computational strategy in
[Ram et al.](https://pubs.acs.org/doi/10.1021/acs.jpca.4c07388), but that paper
does not establish convergence for CH3S/CH3SS. After orbital review, a handful
of QZ-F12 checks at relevant distances/orientations is more informative than
full high-basis recon grids. Check the CC anchor and both MR gaps separately,
including CASPT2 shift/intruder sensitivity and MRCI+Q separated-fragment error.
Retain TZ production if the resulting potential/rate changes meet the chosen
tolerance; do not infer accuracy from agreement of ROHF/CASCI orbitals alone.
An illustrative 1% activated-rate tolerance at 300 K corresponds to only
RT ln(1.01) ≈ 0.0059 kcal/mol. Capture rates need direct sensitivity checks.

For further **CBH/ANL** savings, benchmark a lower reference basis pair using
complete balanced CBH reaction energies across all participants. Total-energy
convergence alone is the wrong comparison. Do not globally switch to ANL0-F12
as a basis-only shortcut: it also changes other recipe components/geometry.
Retain higher-order, CV, relativistic, SO and ZPE contributions until their
reaction-level effect is bounded. For sulfur conventional-basis checks, also
test the native tight-d family; [Molpro distinguishes cc-pV(n+d)Z from cc-pVnZ](https://www.molpro.net/manual/doku.php?id=basis_input).
No custom F12 basis or revision has been substituted.

## Experimental atomic spin–orbit levels

The user approved retaining experimental atomic SO data. The examples explicitly
list the [NIST S I ground-term levels](https://physics.nist.gov/PhysRefData/Handbook/Tables/sulfurtable5_a.htm):
0, 396.055, 573.640 cm−1 with degeneracies 5, 3, 1. KinBot does not fetch them
at runtime. They supply fragment partition functions and the ground-level shift
relative to the spin-free term centroid, −195.7561 cm−1. They do not replace
the calculated molecular doublet–quartet gaps or correlated potential curves.
The declared spatial term degeneracy supplies three spin-free S(3P) components;
it is distinct from the three resolved J levels and their 2J+1 weights.

Spin-free CC/CASPT2/MRCI at TZ-F12 does not yield these J splittings. Replacing
them with calculations would require a separately specified SO Hamiltonian
and state-interaction calculation. The same input schema can accept such
calculated levels with their source, but the production method/basis alone
does not trigger that calculation. Molecular S2 fine structure and CH3S
spin–vibronic structure remain separate, unresolved accuracy questions.

## Software fixes in this follow-up

- Molpro FILE,2 uses short lowercase local WFU names, with an explicit `-W`
  directory for recon, production corrections, ASE forces and ROTD sampling.
  The [Molpro input manual](https://www.molpro.com/manual/doku.php?id=definition_of_molpro_input_language)
  documents 31-character filenames; the observed native build truncated at 32.
  Completed old jobs recover only an unambiguous checkpoint linked to a
  hash-verified native input. A staged predecessor cannot substitute for it.
  The shared I/O retry guard also tracks short WFU names referenced by owned
  inputs, preserving rollback of corrupted attempts.
- Recon audit records concrete validation errors and includes them in controller
  failures, instead of only saying a queued job ended without valid output.
- ANL reconciliation reads results after allocation exit, and allows 60 seconds
  for result publication before declaring a missing-record scheduler failure.
- Recon submission does not consume a live predecessor's WFU. Explicit
  running/failed native records cannot be overridden by legacy export files.
- Any uncertain submission blocks new recon jobs across the manifest, including
  entries that startup/dependency filters would otherwise skip. Scheduler calls
  are bounded; sbatch timeout retains the uncertain journal.
- WFU provenance checks stream SHA-256 instead of allocating the entire file
  in memory, including VRC record validation/staging and force checkpoints.
- CBH fallback distinguishes unavailable tiers from invalid accepted data.
  Invalid energies/states/profile labels now raise instead of silently selecting
  a cheaper tier. Legitimate missing-component fallback remains supported.
- Stopped known allocations can be explicitly archived/retried, and excluded
  basis directories can be archived without deleting their diagnostic evidence.

Validation: the full suite passed **1,321 tests and 391 subtests**. After the
final long-name execution fixture and archive-path guard, the affected recon,
campaign and I/O retry tests passed **79 tests**. All 161 non-template Python
modules parse, all 14 Bash blocks in the current/runbook/readiness documents
pass `bash -n`, and `git diff --check` passes. The actual NH2O SI topology,
electronic counting and convolution verifier passes. The suite warnings are
the pre-existing `imp` deprecation and Sella external-cache permissions.

RHF/ROHF-UCC, intentional ANL extensions, and the arithmetic mean of the two
MR gaps remain intact. A larger basis or averaging two methods is not by itself
proof of higher accuracy. Licensed Molpro/MESS and live Slurm acceptance remain
cluster checks; local regressions do not certify sulfur production chemistry.
