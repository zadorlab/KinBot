# Production observations at 09:21 PDT, 26 October 2026

Source: `examples/anl/tools/root-spectrum-20261026T162140Z/production_monitor.txt`,
header 2026-10-26T09:21:45-07:00, plus the archive's
`production_status_after.json`. Files were captured sequentially, not atomically.
No job or decision was changed during this local audit.

| Area | Observed state | Next dependency |
| --- | --- | --- |
| Shared pool | Cap 23, revision 14; 20 allocations (CH3S 9, CH3SS 11); diagnostic 0; live-control revision 1 | Review controller shares current headroom and exclusions |
| CH3S exploration | Postprocessing; 4 wells; no unstarted wells; 9 retained, 20 omitted by existing policy | Retained thermochemistry and VRC sampling continue |
| CH3SS exploration | Exploring; 8 wells at VRC, 5 exploring; 5/5 exploration slots, no unstarted wells in this snapshot | New wells can still appear; network is not final |
| Orbital diagnostics | Root-spectrum panel 20/20 complete | New 52-job finite panel; native scientific review follows |
| Pending orbital requests | 19 captured: original 12, next 5, newest 2 | Eleven original selections prepared; one original state model unresolved; new seven staged by evidence level |
| Existing VRC | Initial scans/RotdPy and previously reviewed scissions continue; known product-endpoint preparation hold remains on old controllers | Load endpoint fix during final refresh; verify owner dispatch |
| CH3S ANL | Seven expensive allocations still listed: five L3 geometries and two higher-excitation single points | Native convergence/accepted task results determine dependent dispatch |
| MESS | No completed solve established by this evidence | Accepted thermochemistry and required rates/fluxes |

Seven retained ANL allocation elapsed times: CCSDTQ/TZ job 52929237 about
1 d 16 h 39 min; CCSDTQP/DZ 52924181 about 3 d 8 h 45 min; geometry jobs
52928790, 52924269, 52924255, 52924199 and 52918345 span about
1 d 23 h 32 min to 4 d 00 h 24 min. These are elapsed observations, not remaining
time estimates. They are the requested high-level work, not a basis for automatic
simplification. No constrained completion forecast follows from this snapshot.

The four fixed-orbital diagnostic final native real times were 5.45, 3.60,
35.35 and 18.12 seconds (DZ near/far, TZ near/far), excluding queue/admission.
CASSCF optimization in the new panel is different; those times are not a
prediction of its completion time.

## Native failure classification

| Attempt | Native evidence | Disposition |
| --- | --- | --- |
| `...sdc19d1516ced03fc_h2_elim_4_5`, 52938098 | Zero imaginary frequencies | Failed TS-order validation remains; no acceptance override |
| `...sedd8f7729b69b19b_h2_elim_4_5`, 52937703 | Zero imaginary frequencies | Failed TS-order validation remains |
| `...s8fc26f5f496eeb4e_r12_insertion_R_1_2_6`, 52936162 | −39.1781, −14.2173, −3.1246 cm−1; none below configured −50 threshold. Sella warns 11 internal DOF versus 12 Cartesian vibrational DOF. | Earlier unknown is now an order-validation failure; no threshold relaxation or automatic retry |
| `...s254020c87a151164_r12_insertion_R_1_2_3`, 52936698 | 750 Sella steps; final reported fmax about 0.619 eV/Å with repeated oscillation and small trust radius. Last Gaussian force call terminated normally. | Optimization failure, not missing Gaussian termination; repeating the same settings is not justified |
| `...sb6e34ded5baaf03b_h2_elim_5_6_high`, 52930179 | 750 Sella steps; late fmax about 0.071–0.118 eV/Å, final 0.0989. Last UB2PLYP/cc-pVTZ force call terminated normally. | Unconverged L2 optimization; no automatic accepted geometry or downgrade |
| CH+H2S RotdPy surface 115 / face 0 / sample 10, 52937329 | Guard journal: return code 1, `io_failure=false`; retained stream paths recorded, but nested streams omitted by old collector | Native cause unresolved; corrected collector now reads those exact archived streams before selecting recovery |

Seven failed allocation journals were captured across five RotdPy sample IDs.
Some IDs have multiple allocation attempts. The latest same-named launcher
streams are empty and cannot explain the historical failed allocation. A normal
scheduler exit alone cannot establish a completed sample; a missing scheduler
record alone cannot establish failure. Retry accounting must survive owner
restarts and preserve evidence. Those controller fixes remain pending the
native/provenance audit; no blanket sample restart is staged.

The failed TS attempts do not become valid reactions merely because their force
calls terminated normally. Conversely, a failed search is not evidence that the
reaction is kinetically irrelevant. Final rejection/pruning decisions require
the completed networks and the applicable search/recovery history.
