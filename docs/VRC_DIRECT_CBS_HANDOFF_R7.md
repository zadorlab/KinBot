# Direct CASPT2-F12 CBS handoff, r7, 2026-10-22

The user chose direct CASPT2-F12 aDZ-F12/aTZ-F12 extrapolation after the native
reference panel exposed CC/MR correspondence concerns. This replaces the
CC-reference-plus-MR-gap construction for **every existing CH3S/CH3SS VRC
channel**, including all species-specific overrides. It also installs the
reviewed CH + H2S channel from its saved orbital screen in one in-place driver
replacement. No native result is relabelled as a different calculation.

## Scientific changes

At each accepted relaxed geometry, extrapolate the target CASPT2-F12 energies
with power 3.7, then subtract the corresponding target CBS energy at the
separate 30 A asymptote. The correction compares this relative high-level
curve with the existing frozen sampling curve. No CC anchor or spin-gap term
is used. Existing channels keep their CAS, state averages, selected roots,
orbital seeds, shift/IPEA, sampling levels, scan points and dividing surfaces.
Thus their MR scan graphs remain compatible. Their already running QC can
finish; obsolete CC corrections will not enter the new potential.

CH + H2S is newly installed with CAS(5,4), four equally weighted doublets,
D1 reactive, and source natural orbitals 8,12,13,14 from each 2 A checkpoint.
The sulfur p donor has index 7 at 30 A, but the implementation packs the
2 A donor once and propagates its active subspace; it does not reselect index
8 at later geometries. The failed auxiliary quartet reference is excluded
from this new model's state average. See [the orbital evidence](CH_H2S_ORBITAL_REVIEW.md).

The direct model is an explicit accuracy tradeoff. It does not cure an
inadequate active space, root change or H transfer. The C + SH3 whole-scan
connectivity hold remains, with its bounded retries and independent work
continuation. Final rates still require those scientific reviews.

## Single handoff on Blodgett

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
run="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
proposal="$run/ch_h2s_direct_model_proposal.json"
"$HOME/KinBot/.venv/bin/python" examples/anl/tools/requested_orbital_recon.py proposal \
  "$HOME/kinbot-runs/r7_ch_h2s_orbital_recon_20261022" \
  "$run/ch3s/input.json" \
  examples/anl/ch3s_anl1_f12/orbital_reviews/ch_h2s_cs_selected.json \
  "$proposal"
"$HOME/KinBot/.venv/bin/python" -m kinbot.run_refresh "$run" \
  --vrc-direct-caspt2-cbs --vrc-reviewed-models "$proposal" --apply
"$HOME/KinBot/.venv/bin/python" -m kinbot.run_monitor "$run"
BASH
```

The proposal step validates request identity, source geometry, all eight
completed orbital-screen jobs and native WFU/output hashes without publishing
a live edit. The refresh validates all mapped and initial channels and both
correction bases before pausing controllers, then repeats the proposal/evidence
checks while controllers are quiesced. Changed inputs/evidence or unresolved
live input requests stop the migration before input edits. Case inputs,
resolved settings, every existing well input and the persistent
`.kinbot_search_options.json` channel overrides are updated together. Inputs
are archived; databases, scan records, native QC and submitted jobs remain at
their original paths. New controllers adopt saved jobs. The final audit prints
`Verified current VRC inputs` and the direct CBS contract after launch.

An interrupted refresh must be resumed with the **same run_refresh command and
proposal file**, without regenerating the proposal or changing code mid-handoff.
Once completed, a same-code repeated refresh command is a no-op. Historical
native inputs, old proposal archives and failed diagnostic recipes intentionally
retain the settings that actually generated their evidence.

## Future spin-splitting review

New spin-splitting orbital selections and changed contracts at full refresh
must include `vrc_tst_reference_review: {"path": "...", "sha256": "..."}`.
An unchanged legacy model is not retroactively certified by this rule.
Direct selections do not require a CC reference review.

The receipt requires completed validated native CC and MR energies plus
state-specific MR/ROHF density matrices at an interior point and the configured
asymptote in every reference/gap basis. It binds the actual CAS, orbital records,
state average, selected roots, correction levels and shifts. Native execution
must record the contract fingerprint established before staging. The existing
production probe builder records that fingerprint; its current DZ-only panel
is not sufficient for a two-basis acceptance. Never retrofit a fingerprint onto
old results to make them pass.

For an explicitly prepared scalar-leaf parameter file:

```bash
python -m kinbot.vrc_reference_review contract leaf.json
# Stage matching native probes with that reference_review_contract_sha256.
# checks.json lists {"manifest": "/absolute/reconnaissance.json", "job": "..."}.
python -m kinbot.vrc_reference_review leaf.json checks.json draft.json
python -m kinbot.vrc_reference_review leaf.json checks.json approved.json --findings findings.json
```

`findings.json` must contain explicit, nonempty reviews of `state_pairing`,
`cc_quality`, `degenerate_subspaces` and `orbital_continuity`. The approval
command records the user's scientific decision; density similarity alone
cannot make it. The gate reruns the numerical analysis and verifies hashes
before publication. It rejects missing, failed, changed or differently bound
checks; it does not bypass CC convergence errors or auto-select the closest
root. Attach the printed path/hash to the proposed stage or scalar input.
