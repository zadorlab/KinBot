# R8 bindings: independent closure decisions, 19 November 2026

B–W3 is one separate native calculation. It does not gate this evidence review.
The remaining work has been separated into recovered provenance, counting rules,
proposed electronic models, and unresolved physical connections.
**R9 is not ready for a full production launch.**

## Attached driver

`git pull` and preflight passed in the supplied Blodgett output. The original
`drive` command waits for the task to finish and prints its summary only at the
end. That output does not establish a hang or confirm a Slurm submission.
The updated driver prints attachment immediately, task/job changes, and an
unchanged-state heartbeat every 60 seconds. It keeps the existing owner and
submission rules.

Read the current run from a second terminal:

```bash
cd "$HOME/KinBot" || exit
KB_CLOSURE="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/validation_20261114/w3_complex_b_saddle_20261119"
.venv/bin/python -B -m kinbot.anl.dispatch status "$KB_CLOSURE" --cached
cat "$KB_CLOSURE/state.json"
squeue --me -h -o '%i|%T|%R|%j|%Z'
```

`submitted` plus a matching queue row means Slurm owns the calculation; its
state/reason distinguishes running from pending. `staged` has not been admitted;
`submitting` requires reconciliation before another submission. A missing queue
row alone is not a calculation result. The running controller keeps its loaded
code after a pull; an update is not a reason to restart a live calculation.

## Closed decisions

| Item | Existing evidence | Action |
|---|---|---|
| Three reverse-saddle mode records | Native imaginary modes −1157.3337, −1077.8991 and −1287.1405 cm⁻¹; atom-frame errors below 0.000001 Å | Bind the recovered native modes and source hashes; no frequency retry |
| `rxn_11` / reflected branch counting | Self-mirror saddle; even unstable mode, overlap 0.999999685; existing weights 0.5 + 0.5 | Preserve both explicit-population branches with total surface weight one |
| Global `rxn_4` / `rxn_36` deduplication | Cross-matching saddles are attached to distinct enantiomer populations | Do not collapse them merely because geometry and mode match; their endpoint assignments remain a separate question |
| CH2SH H-loss inner path (`h101`) | Existing two-sided B2PLYP path and endpoint spectra; duplicate seed already identified | Retain one inner saddle, reuse the completed outer electronic evidence |
| TS5 old H-exit asset reuse | All saved H-exit cases 14/15/16/29 leave a chain C–S–S fragment | Reject binding those assets to HC(=S)SH + H |
| Optional CC and completed orbital campaigns | Their failures/completions do not change the established B2PLYP connectivity conclusions | No blanket retry, no repeated reconnaissance |

The original TS5 target fragment has C–S distances **1.62484 and 1.73766 Å** and
S–S **3.03041 Å**. The checked chain fragments have one C–S at 1.55–1.61 Å,
the other at 2.97–3.36 Å, and S–S at 2.13–2.27 Å. No sulfur permutation or
rigid rotation changes these distance sets. The shared formula/name “HCSSH”
is insufficient for checkpoint reuse.

## Remaining bindings and concrete dispositions

| Item | Finding | Disposition |
|---|---|---|
| Legacy `rxn_reverse_93359f49` | Saved entrance matches complex B within 0.000194 Å RMSD; serialized edge names A. Its new reverse L2 IRC is unconverged | Reject the A assignment. A B reassignment is provisional until the reverse path is supported. Retain the independently validated paper-TS1 A route |
| Legacy `26b8db7c` and `b18dcbc8` | Native first-order mode provenance exists, but the supplied path failures remain | Do not count mode recovery as new connectivity validation. Establish retained-network necessity before requesting another path |
| `rxn_4` / `rxn_36` handed endpoints | L1 discovery histories do not independently verify the L2 handed assignments | Keep the extra routes conditional; do not call their population mapping settled |
| Mirror C–S capture | One hand uses two doublets + two quartets; the later opposite-hand review uses one doublet + one quartet | Concrete canonical proposal: later CAS(5,4), 1D+1Q, D1 flux for both hands. Preserve separate populations; confirm replacement of the older state model before applying |
| Mirror S–S capture | One hand uses CAS(5,4), the later opposite-hand review CAS(7,5), both 3D+3Q | Concrete canonical proposal: later CAS(7,5), 3D+3Q, D1 flux for both hands. Confirm the model change; no new orbital screening |
| `h101` outer choice | One-doublet CAS(3,3) return checks preserve the subspace; near energy 6.06291 kcal/mol at aTZ. Two-doublet CAS(3,3) versus CAS(5,5) near difference 5.47260 kcal/mol | Prepare the existing one-doublet CAS(3,3) outer candidate, retain the established inner TS, and keep the human selection explicit |
| `h102`: W4 → CH2SS + H | Frozen two-doublet CAS(3,3) and CAS(5,5) near energies differ by only 0.08212 kcal/mol at aTZ; the relaxed B2PLYP scan switches chemistry | Reuse the completed outer evidence. A retained direct W4 H-loss channel still needs a state-consistent inner path; neither frozen points nor the contaminated scan establish it |
| `h102` CAS(5,4) DZ return | Energy difference only 0.0003263 kcal/mol, but minimum squared active-subspace singular value 0.723112 | Reject an energy-only tracking pass for that control; keep better-tracked existing candidates |
| TS5 separated exit | W3 departure, CH(SH)2 recapture and HC(=S)SH + H are distinct pieces. The saved chain-fragment H-exit model is inapplicable | Preserve the bound inner route. An exact HC(=S)SH + H outer model is still required if that separated exit is retained; H101/H102 do not supply it |

All radial energies use each model's own r30 reference, CASPT2-F12 with IPEA
0.25 and shift 0.1. They are frozen-fragment interaction energies, not relaxed
activation barriers. The 578-point campaign finished, but its final record
explicitly leaves production-model acceptance false. The model proposals here
preserve that distinction.

## Kinetics and R9

The provenance-enriched input has **identical physical tokens** to the
conditional CH3SS model already run natively. The corresponding execution
receipt remains applicable. This revision does not fabricate new rates or
change chemical connectivity in that input.

The earlier native MESS results still contain negative individual-species
finite-pressure coefficients. Positive unreduced propagation establishes a
usable numerical route for tested time-dependent observables, not acceptance
of those reduced rate tables. The final chemical choices and the desired
initial radical energy/distribution, bath and observation time must precede
validation of the final R9 observable. The 148.142 kcal/mol ground-state radical
scout is a declared model scenario; full photon-energy deposition is not an
experimentally established energy partition.

The remaining consequential work is therefore:

1. Review the independent B–W3 native saddle result when available.
2. Resolve the retained H102 inner path and exact TS5 outer-model coverage;
   completed orbital controls and optional CC retries do not fill these gaps.
3. Apply approved electronic-model choices and select/reconcile the conditional
   legacy and handed routes using their actual impact on the selected network.
4. Regenerate the final network, then validate its requested MESS observable
   with the appropriate grain/time/reduction checks before staging production.

No additional remote calculation is submitted by this evidence package.

## Evidence and verification

- `binding_decisions.json`: per-issue status, complete proposed mirror
  contracts, basis-specific seed hashes, H-loss model energies and tracking.
- `h_loss_fragment_binding.json`: exact native-source and target geometries,
  atom mappings, invariant distances, comparison for both basis assets.
- `network_binding_evidence.json.gz`: native mode/frame/provenance and
  population-aware comparisons.
- `h_loss_model_evidence.json.gz`: complete joined point records, input/output/
  WFU hashes, state definitions, natural occupations and return checks.
- `ch3ss_provenance_bound.inp.gz`: conditional MESS input with recovered mode
  metadata; chemical tokens unchanged.
- `manifest.json`: hashes for this review package.

A separate verification pass rehashes sources, checks geometry distance
invariants, validates mode counts and physical-token identity, and recomputes
reported energy/return differences from stored records. This is a same-reviewer
verification; no additional agents or production-model approvals were used.
