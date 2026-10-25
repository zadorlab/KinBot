# Remaining CH3SS orbital requests, 2026-10-25

This batch contains the four requests from parent `791871710370000000002`.
The five previously reviewed CH3SS channels need accepted-model reconciliation,
not another orbital screen. No production model is published by this batch.

Atom numbers are one-based and follow the request order C1, S2, S3, H4, H5, H6.

| Request suffix | Physical channel and separation | Diagnostic space and roots |
| --- | --- | --- |
| `hom_sci_1_2` | Parent C1–S2 cleavage → CH + H₂SS | Heavy-atom CASCI(19,12); two equally weighted doublets support the CH orbital pair. |
| `hom_sci_2_3` | Parent S2–S3 cleavage → CHSH₂ radical + S | Heavy-atom CASCI(19,12); three doublets and three quartets support the S spatial manifold coupled to the doublet fragment. Quartet ROHF initialization. |
| `R_Addition_CSm_R_5_1_2` | C(SH₂) + SH product-complex exit; S2–S3 separation | Full-valence CASCI(19,15), two doublets; product fragments [1,2,4,6] and [3,5]. |
| `r12_insertion_R_2_1_5` | HCS + H₂S product-complex exit; C1–H4 separation | Full-valence CASCI(19,15), two doublets; product fragments [1,3,5] and [2,4,6]. |

The exit screens use normal L2 `_IRC_R_prod_high` records: row 843 for
complex `791430280000000000002`, and row 1162 for
`791230350000000000002`, as supplied on 2026-10-25. Staging reads the latest
record with the explicit name, requires normal completion, checks atom order,
the accepted network complex, fragment connectivity, charge/multiplicity and
the separation coordinate. It saves the actual row ID and geometry in
`geometry_source.json`. The original parent geometry is not used for these
exits. Internal distances within each fragment remain fixed during translation.

Each screen has independent AVAS initializations at 2, 3.5, 6 and 30 Å in
avdz-f12 and avtz-f12: 32 fixed-geometry jobs altogether. Heavy-atom targets
are C 2s/2p and both S 3s/3p; the exit screens also include all three H 1s
targets so hydrogen-containing fragment partners are available for inspection.
These broad spaces are reconnaissance, not final compact CAS choices. Their
state averages do not specify how many reactive fluxes enter production.
No additional relaxed-geometry reconnaissance is requested here.

The two S2–S3 requests have different product identities and separate evidence
directories. Their later production publication must distinguish those
reactions explicitly: the current bond-only channel selector cannot represent
both as separate models in one parent map. Do not approve one screen as evidence
for both. This batch deliberately performs no model publication.

## Live fix and launch on Blodgett

Pull the commit identified in the accompanying handoff, then:

```bash
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_PY="$PWD/.venv/bin/python"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_REVIEW="$KB_RUN/orbital_recon/remaining_20261025"
KB_PLAN="$PWD/examples/anl/orbital_review_batches/r7_remaining_20261025/screens.json"

"$KB_PY" -B -m kinbot.run_control resubmit "$KB_RUN/ch3ss/input.json" --reconcile-models
"$KB_PY" -B -m kinbot.run_control resubmit "$KB_RUN/ch3ss/input.json" --reconcile-models --apply
if [ ! -f "$KB_REVIEW/review_launch.json" ]; then
  "$KB_PY" -B -m kinbot.requested_orbital_recon batch \
    "$KB_RUN" "$KB_PLAN" "$KB_REVIEW" --shared-pool
fi
KB_LOG="$KB_REVIEW/run.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$KB_PY" -B -u -m kinbot.requested_orbital_recon run \
  "$KB_REVIEW" --report-failures >"$KB_LOG" 2>&1 </dev/null &
printf 'Review PID: %s\nLog: %s\n' "$!" "$KB_LOG"
```

The first reconciliation command previews; the second queues accepted models
through the current owning drivers. It does not alter scientific choices,
cancel QC, or restart controllers. Previously searching workers inherit at their
ready boundary. The new accepted-snapshot inheritance fix is also retained in
the source for future controller generations. A queued/applied receipt alone is
not proof of native scan completion: check the worker review records and monitor.

Shared-pool staging creates a fresh review directory under the production
policy. Each submission uses the existing admission lock, ledger and exclusions.
The review cap is three within the existing total of 23, with no separate
reservation and no policy edits. Full capacity leaves the remaining jobs
prepared until a slot becomes available. Existing production allocations stay
owned by their drivers. A second review owner is refused by the campaign lock;
rerunning after owner exit adopts existing job IDs and results.

```bash
"$KB_PY" -B -m kinbot.run_control resubmit "$KB_RUN/ch3ss/input.json" --status
"$KB_PY" -B -m kinbot.requested_orbital_recon monitor "$KB_REVIEW"
"$KB_PY" -B -m kinbot.requested_orbital_recon evidence "$KB_REVIEW"
"$KB_PY" -B -m kinbot.run_monitor "$KB_RUN"
```

Local validation covers the supplied geometries/decks, fragment-preserving
translations, all-or-nothing preflight, shared admission with retained work,
timeouts, purged terminal attempts, model reconciliation across generations,
and repeated invocation. Scheduler/native outcomes on Blodgett remain to be
verified from the resulting receipts and output.

## Interpretation of the 08:47 PDT monitor

CH3S L1/L2 exploration was complete. Its four retained VRC surfaces each had
48 accepted scan points including 30 Å, and their four RotdPy owners were
sampling. ANL geometry/correction work was still running; MESS was waiting for
required results. The operator-rejected C + SH₃ channel remains excluded.

CH3SS was still exploring: five active exploration workers and one unstarted
well. Its initial C–S surface had 40 accepted points (through 13 Å; 14 Å running),
and S–S had 19 (through 5.6 Å; 5.8 Å running). The nine review holds comprised
the five missed inherited models and these four new requests. New outer-selected
CBH/ANL work was still behind the L2 exploration gate. The snapshot had 22
allocations against a 23-job pool.

This is not yet the final combined PES/CBH network. Final kinetic pruning,
removal of orphan products, cross-PES reuse decisions and optional expensive-task
simplifications remain after both L1/L2 explorations finish, as requested.
Pruning must retain every surviving reaction's products and the transitive CBH
reference dependencies, including references shared with another target/PES.
Outer L2 path selection is not itself evidence of negligible kinetic flux.
