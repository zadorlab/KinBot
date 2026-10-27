# Native results and compact follow-up, 26 October 2026

The transferred root-spectrum archive verifies: **860 files, zero collection
issues**, SHA-256 `005278deac5b8d5cb29ec425c52c6e37c905f167a54a2fb56544c4faa9a62e4f`.
All **20 jobs completed**. Their native results are retained in
[native_character_evidence.txt](native_character_evidence.txt),
[root_spectrum_audit.json](root_spectrum_audit.json), and the local archive
`examples/anl/tools/root-spectrum-20261026T162140Z/orbital-evidence.tar.gz`.
Verification proves file integrity, not scientific acceptance.

**Next finite panel: 52 jobs, cap 3 within the production pool.** Forty-four
CASSCF jobs test explicitly identified compact spaces; eight SCF jobs identify
orbitals for two further requests discovered in the latest capture. None of
these jobs publishes a model. Eleven of the original twelve selections remain
prepared. Five requests from the next parent and two from the newest parents
are also tracked; no request has been silently dropped.

## What the fixed-orbital calculation resolved

For the C(SH2)+SH product exit, all four CAS(9,6) CASCI calculations retained
exactly the intended six-dimensional orbital span: minimum overlap singular
value is at least **0.999999999999987**. This distinguishes the earlier CASSCF
orbital substitution from seed-transport loss. At 30 Å, the two lowest doublets
form the SH-component pair in these fixed orbitals. Their splitting is
**0.239 cm−1 (DZ), 0.217 cm−1 (TZ)**. D3/D4 coincide with Q1/Q2 to the printed
precision and their CI vectors have three-open-shell character; D3 lies
23,358/24,698 cm−1 above D1. Those are fixed-orbital diagnostic gaps, **not
isolated-fragment adiabatic singlet–triplet energies**.

Near the complex, Q1 is between D2 and D3. D1 has substantial three-open-shell
configuration mixing while D2 is dominated by `222a20` in the new diagnostic
natural basis. Root numbers therefore cannot be assumed to follow the same
fragment state from near to far. The earlier two-doublet optimizations lost a
required direction; repeating them is exhausted.

The next comparison optimizes four doublets and two quartets with **1/6 weight
per state**. These are orbital/state-identity support, not six reactive flux
contributions. State averaging over interacting states is supported by the
[Molpro MULTI documentation](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#defining_the_optimized_states);
the chosen six-state window is a hypothesis derived from this native spectrum.
It does not guarantee retention or authorize changing the requested fragment
multiplicities. After optimization, examine CI vectors, active/inactive orbital
character, near/far correlation and asymptotic splittings before selecting
reactive roots. A changed state ordering would require explicit review.

Two spaces are compared from the **same geometry-specific, basis-specific
CASCI natural checkpoints**:

* **CAS(5,4), source natural 19–22:** the SH p pair plus carbon frontier pair.
  Closing source 17/18 removes two occupied spectator directions. Their
  diagnostic occupations are about 1.996/1.980 near and 1.999/1.980 far;
  occupation alone is not the justification—the AO character and state
  configurations identify the retained frontier manifold. Test their necessity
  against the larger reference after orbital relaxation.
* **CAS(9,6), source natural 17–22:** the existing six-orbital union, retaining
  the two occupied C/S valence directions as the comparison. Its orbital count
  is unchanged. The averaging support changes from the exhausted two-doublet
  trial; no all-valence space is introduced.

Both start independently at the unchanged L2 complex endpoint and 30 Å in
avdz-f12 and avtz-f12: eight jobs. These are optimized CASSCF probes, unlike
the completed fixed-orbital diagnostics. No orbital is frozen.

## New parent: minimal-space comparisons

Parent `791871560500020000002-sb6e34ded5baaf03b`; atom order C1 S2 S3 H4 H5 H6.
All indices below refer to **this channel's own 30 Å canonical SCF record
2100.2**, separately hashed for each basis. Transport reads a copy on 2100.3;
optimized natural orbitals are written to 2140.2. Packed indices are bookkeeping,
not proof that the initial character survives CASSCF.

| Reaction / fragments | Probe and source indices | Why included; state average | Jobs |
| --- | --- | --- | --- |
| C1–S2 `hom_sci_1_2`: CH=(1,6), HSSH=(2,3,4,5) | (3,3): 18,21,22; (5,4): 18,20,21,22, both bases | CH occupied sigma/frontier manifold and two carbon p directions; compare addition of occupied sulfur donor combination 20. Two equal doublets for CH orbital components. | 8 |
| S2–S3 `hom_sci_2_3`: HCSH=(1,2,4,6), SH=(3,5) | (5,4): DZ 18,19,21,22; TZ 18,19,21,23. (7,5) adds 20. | SH p pair 19/21; occupied C/S p direction 18 and the virtual with opposite C/S valence-p signs (DZ22/TZ23); compare occupied carbon-rich sigma direction 20. Two equal doublets. | 8 |
| `12_shift_S_F_3_2_4` and `12_shift_S_R_3_2_4`: HCS=(1,2,6), H2S=(3,4,5) outer exit | (3,3): 19,21,22; (5,4): 19,20,21,22, both bases | HCS pi/SOMO/pi-star; compare H2S sulfur lone pair 20. One doublet. Exact shared L2 endpoint is checked again; separate reaction bindings remain. | 8 shared |
| `r12_insertion_R_4_2_1`: HCSSH=(1,2,3,5,6)+H4 outer exit | (1,1): 21; two separate (3,2) controls: 19,21 and 20,21 | H4 1s SOMO; compare S3-rich donor 19 and C/S frontier donor 20 separately. One doublet each. | 12 |

Scission probes use 2 and 30 Å. Product exits use their actual L2 endpoint and
30 Å, generated by rigid fragment translation. C1–H4 and H5–H4 are outer
separation anchors, not parent covalent bonds being cleaved. No fragment
internal bond or geometry is relaxed by the panel.

The SCF evidence distinguishes the DZ/TZ virtual ordering for HCSH+SH; taking
the lowest virtual in both bases would miss that distinction. Printed AO
coefficients are amplitudes, not populations. Final decisions require optimized
natural character, CI/root continuity and endpoint-relative energy comparisons,
not merely successful termination or a nearly 2/0 occupation.

Packed active windows are 21 for (1,1); 20–21 for (3,2); 20–22 for (3,3);
19–22 for (5,4); 18–22 for (7,5); 17–22 for (9,6).

## Two further parents discovered in the export

Both are C1–S2 homolysis to HCSH=(1,3,5,6)+SH=(2,4), with their own parent
geometry and sulfur atom assignment:

* `791671560700020000002-sdc19d1516ced03fc_hom_sci_1_2`
* `791671560700020000002-sedd8f7729b69b19b_hom_sci_1_2`

[Captured requests and network records](new_requests_evidence.json) are checked
against the live files. Each gets SCF character at 2/30 Å in both bases—eight
jobs total. No CAS indices are transplanted from the S–S scission above, despite
the common product formula. These requests were not present in the prior
five-request staging inventory.

## Blodgett launch and collection

From your SSH session, after pulling the accompanying commit:

```bash
cd "$HOME/KinBot"
git pull --ff-only origin composite
bash examples/anl/orbital_review_batches/r7_state_support_20261026/run_blodgett.sh
```

This validates all seven new request definitions, checks endpoint sharing,
preflights all 52 decks and every source hash, stages a fresh sibling campaign,
and starts its review owner. Admission uses the current shared cap and excluded
nodes, with at most three review allocations. It preserves completed panels,
production jobs, live decisions and sampling. Repeating the command resumes
unchanged work; changed scientific settings are refused. No production model,
PES restart or cancellation is requested here.

It also writes `production_status_before.json`. The corrected collector follows
failed allocation journals to their **archived attempt stderr/stdout**. Transfer
that JSON immediately if convenient; it need not wait for the new calculations.
The preceding capture included journals but omitted these nested streams, so
it could not establish the native cause of RotdPy job 52937329.

Monitor:

```bash
"$HOME/KinBot/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/orbital_recon/state_support_20261026"
```

Collect after completion:

```bash
bash "$HOME/KinBot/examples/anl/orbital_review_batches/r7_state_support_20261026/collect_blodgett.sh"
```

Transfer the printed export directory to local `examples/anl/tools/`. It includes
native outputs, WFU/XML, status and request evidence. Previous verified exports
remain the source record; they need not be transferred again.

## Production status and remaining handoff

The [latest status audit](status_audit.md) records the 09:21 monitor and the
captured failures. CH3SS is still exploring; final PES/CBH pruning is premature.
The unresolved tasks remain explicit:

1. Review the optimized six-state comparison and the new compact candidates;
   identify compact indices for the two newest parents from their SCF evidence.
   No unreviewed root or orbital choice is automatically published.
2. Finish the RotdPy allocation-recovery audit using actual failed-attempt streams.
   In-memory retry counts and ambiguous scheduler handling remain known code
   concerns. This panel fixes evidence collection, **not** that controller logic.
3. Combine accepted choices with the eleven prepared selections; reconcile the
   five previously published channels. Preserve the operator's excluded channel.
4. Complete/test the controller-refresh handoff, load the endpoint fix, preserve
   live jobs/limits/exclusions and sampling ownership, and verify actual worker
   scan dispatch/adoption beyond the case-level receipt. The endpoint defect
   calls for VRC preparation, not earlier PES/IRC/L2/ANL/CBH reruns.
5. Retain the full existing/new active-space explanations, native hashes and
   scientific limits with the final publication. The previous ANL/MESS/live-choice
   work remains in scope; this export does not provide new final MESS validation.

No complete production handoff or successful native recovery is claimed by this
update. The local tests cover all 52 plans, source preservation, shared-pool
admission with live jobs, repeated invocation, and failed-stream capture.
# Detached execution after retirement of the source run

The source `_r7` is now retired. Use `run_detached_blodgett.sh` instead of
`run_blodgett.sh`. It stages the same 52 explicit probes beneath
`~/kinbot-runs/ch3s_ch3ss_orbital_reviews_20261027/state_support`, with a separate
three-job resource policy and inherited node exclusions. It reads the original
requests, completed outputs and hash-bound checkpoints without editing them.
The general `requested_orbital_recon --resource-root` option separates the
evidence tree from the scheduling owner; it does not copy or restart PES workers.

The launcher also captures **all** current orbital-review requests, their
installed channel definitions, accepted reaction records and both IRC endpoint
directions. Energy pruning and orbital approval remain independent: this panel
still includes the C(SH2)+SH channel even though the 150 kcal/mol audit excludes
its asymptote. Further pending requests are reconciled against the inventory;
the seven-screen panel is not a claim that every requested model is approved.

Local simulated-site tests exercise all 52 jobs, admission at the new pool,
unchanged source files, unrelated live work, and repeated invocation. Native
results and scientific CAS acceptance remain pending until the remote probes
finish and their evidence is reviewed.
