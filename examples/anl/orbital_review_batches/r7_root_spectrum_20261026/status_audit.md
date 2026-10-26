# Evidence and production status, 26 October 2026

The [08:30:20 PDT monitor](production_monitor_20261026_083020.txt) is retained
verbatim, SHA-256
`ebc5bd87f6914d67f9a3be9d72075052560ffead379b43c9a19b45ede6b11861`.
It collects observations over time; later log lines can postdate the header.
Per-worker and overview counts therefore need not be simultaneous.

| Area | Observed result | Next dependency |
| --- | --- | --- |
| Pool | Cap 23, revision 14; 22 allocations, CH3S 9 / CH3SS 13; diagnostic 0; live-control revision 1, 16 applied requests | Review admission shares this cap and current exclusions |
| CH3S exploration | Postprocessing; four wells, zero unstarted; 9 retained / 20 omitted by pathway policy | Continue retained thermochemistry and sampling |
| CH3S CBH-ANL | Across 246 displayed graphs: 1,118 complete, 29 skipped, 325 staged, 9 submitted, 284 waiting | Owners dispatch dependencies as accepted results arrive; these are task counts, not species counts |
| CH3S VRC | Four retained state scans resolved, each 48 accepted entries through 30 Å | Four RotdPy campaigns running; final flux validation remains |
| CH3S sampling | Allocations include surfaces 101 and 108; final surface files remain 0/116 per campaign | Intermediate sample work is separate from final surface convergence; inspect failure records below |
| CH3SS exploration | Seven wells marked VRC, five exploring; 4/5 exploration slots, two unstarted wells | New wells/requests still appear; final PES/CBH pruning remains premature |
| Initial CH3SS C–S | 48 accepted entries, resolved; RotdPy running, surfaces 58/59 visible | Sampling/flux validation |
| Initial CH3SS S–S | 32 accepted, 9.0 Å running after 8.5 Å | Remaining scan/asymptote and orbital/root/connectivity review |
| Reviewed CH2SHS S–S | 12 accepted, 4.4 Å running after 4.2 Å | Remaining scan and review |
| Same parent's C–S | One accepted, 2.2 Å running; allocation elapsed about 18 h 25 min | Native iterations needed to assess progress; elapsed alone is not a stall verdict |
| Earlier CH2(SH)S exits | Old controller still has the known endpoint-preparation hold | Final code refresh and VRC preparation; no earlier chemistry rerun indicated |
| New orbital selections | Eleven choices prepared; C(SH2)/SH unresolved; five additional requests received | Four state diagnostics and four distinct SCF geometry screens, then compact-space review |
| CH3SS CBH-ANL | New work waits for outer PES L2 selection | Complete exploration/selection before final combined network assessment |
| MESS | Waiting, executable ready | Accepted thermal and required rate/flux data; no completed solve in this evidence |

The operator-rejected C+SH3 channel remains excluded. This audit makes no new
reaction pruning, thermochemical simplification or allocation cancellation.
Sharing and per-species/per-task ANL controls remain available; their use must
retain other species' dependencies and reference work.

## Time observations

Seven expensive CH3S allocations are shown: five L3 geometry jobs and two
higher-excitation single points. Geometry elapsed times range from **1 d 22 h
41 min to 3 d 23 h 33 min**; CCSDTQ/TZ is **1 d 15 h 48 min**, and
CCSDTQP/DZ is **3 d 7 h 54 min**. Other previously listed allocations are absent;
that alone does not prove their tasks completed. Graph/attempt records are
needed to distinguish completion, failure and replacement.

The completed natural-start panel reports final Molpro real times of roughly
10/36 s for DZ near/far and 68/354 s for TZ near/far, excluding queue admission.
The new fixed-orbital and SCF work is different; no native duration or completion
time is asserted. CH3SS exploration has no constrained finish estimate while
new and unstarted wells remain.

## Failures that need native classification

The capture reads these attempts without accepting or retrying them:

* H2 elimination `...sdc19d1516ced03fc_h2_elim_4_5`, job 52938098, and
  `...sedd8f7729b69b19b_h2_elim_4_5`, job 52937703: optimized structures have
  **zero imaginary frequencies**, so TS acceptance correctly failed. The
  monitor does not establish whether another applicable recovery exists.
* `...s254020c87a151164_r12_insertion_R_1_2_3`, job 52936698: 750-step limit
  with large residual forces. The trajectory needs diagnosis before a retry.
* `...s8fc26f5f496eeb4e_r12_insertion_R_1_2_6`, job 52936162, and
  `...sb6e34ded5baaf03b_h2_elim_5_6_high`, job 52930179: the compact monitor
  lacks enough native evidence to classify the failure or choose recovery.
* CH+H2S RotdPy surface 115 / face 0 / sample 10, job 52937329: Slurm FAILED
  at 04:57. Other warnings concern missing scheduler records on surfaces
  113–115. Scheduler disappearance is not a native calculation diagnosis.
  Capture stderr/stdout, retry journals, campaign ownership and worker evidence.

Code inspection also identifies RotdPy allocation-retry accounting held in
memory and a scheduler-unknown branch that treats an attempt as failed. These
need a provenance-preserving, restart-aware audit against the captured attempts
before the final refresh. No fix or successful recovery for those branches is
claimed by this diagnostic update. Native Molpro solver-recovery tests are
distinct from that allocation/controller issue.

Historical CH3S failed exits remain as evidence, not automatic retry requests.
Product-minimum descent recovery is recorded separately; it does not make a
zero-imaginary structure an accepted TS.

## Verified native archives and reproducibility

| Local export | Archive SHA-256 | Verification |
| --- | --- | --- |
| `examples/anl/tools/compact-continuation-20261025T232718Z/orbital-evidence.tar.gz` | `7db5f41aee419391295616af5fd029ed02ab36ff93f2a013548f2f611e8a63f1` | 4,056 archived members verified previously |
| `examples/anl/tools/natural-union-20261026T012004Z/orbital-evidence.tar.gz` | `8621cbfed310ed09d149cba2d0b6c71cd2b8f6a0478dc5b8afb4d1208c7102fa` | 731 files verified, zero collection issues |

The newer archive contains the preceding ten-job compact-union panel and four
natural-start jobs. Its two CH2/HSS selected output hashes are now bound in the
prepared selection alongside the WFU hashes. Reproduce the numerical audit with:

```bash
PYTHONPATH="$PWD" .venv/bin/python -B \
  examples/anl/orbital_review_batches/r7_root_spectrum_20261026/audit_natural_union.py \
  /path/to/extracted/compact-continuation /path/to/extracted/natural-union
```

The audit includes native state energies, active/inactive character, CI tables,
runtime entries, same-geometry MO comparisons and file hashes. Each basis uses
its own checkpoints. Archive integrity is separate from scientific approval;
the last C(SH2)/SH production selection remains unresolved.
