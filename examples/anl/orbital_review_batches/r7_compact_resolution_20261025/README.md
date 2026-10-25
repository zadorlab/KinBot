# Review of the completed compact panel and next minimal comparisons

The user supplied `compact_review_character.20261025T203418Z.txt` on 25 October
2026. It contains **24 completed compact CASSCF probes and 32 completed SCF
character probes**. The complete supplied text is saved in
[native_character_evidence.txt](native_character_evidence.txt); its SHA-256 and
scope are recorded in [evidence_provenance.json](evidence_provenance.json).
WFU and XML hashes are included in that report. The binary checkpoints and full
native outputs remain on Blodgett; the pasted top-six AO coefficients are not a
complete coefficient matrix or a quantitative orbital-overlap calculation.

The six affected seed jobs now have native completions. In particular, the
recomputed HCS/H2S CAS(3,3) contains the intended valence orbitals, whereas the
previous forced-OVL result contained sulfur core-like character. This supports
the same-basis restart repair for these native attempts. It does not establish
that every optimized active space is chemically consistent.

This directory stages **60 fixed-geometry CASSCF jobs, 15 candidates across ten
requests**: 52 jobs for the eight newly characterized requests and eight jobs to
resolve two specific active-space changes. Every trial has **two to five active
orbitals**. There are no new SCF-only screens, full-valence CI calculations,
CASPT2 calculations or relaxed scans in this panel. It uses the existing general
`kinbot.requested_orbital_recon` driver; no installed code changes are required.

## What the completed CASSCF probes establish

For parent `791871710370000000002`, atom order is C1 S2 S3 H4 H5 H6. Source
indices below refer to **canonical SCF record 2100.2**, with separate immutable
WFU hashes for each basis. Optimized natural indices refer to **record 2140.2**;
they must not be substituted for source indices when building a new guess.
All statements about character below are interpretations of the printed native
AO coefficients and occupations, not orbital population integrals.

| Reaction and coordinate | Source choice | Observed optimized subspace | Assessment |
| --- | --- | --- | --- |
| `hom_sci_1_2`, C1–S2 → CH + H2SS | CAS(5,4), 17,18,21,22 from 30 Å; two equally weighted doublets | Active 19–22 retain the sulfur donor and CH frontier subspace. At 30 Å: sulfur donor 19 ≈2.000; CH occupied frontier 20 ≈1.954; CH p pair 21/22 ≈0.523 each. Near occupations are ≈1.983,1.497,1.487,0.032 in both bases. | Supports this orbital selection at the tested endpoints in both bases. No additional partition trial is staged. |
| `hom_sci_2_3`, S2–S3 → CHSH2 + S | CAS(5,4), 18,20,21,22 from 30 Å; three doublets and three quartets | At 30 Å active 19–21 span all three S3 p directions, each ≈1.33333; 22 is the fragment radical, ≈1.000. Near occupations ≈1.499,1.499,1.003,0.999 retain sulfur/fragment frontier character in both bases. Far doublet/quartet roots agree within about 10^-6 Eh. | Supports the atomic-S manifold plus fragment-radical selection. No additional partition trial is staged. |
| `R_Addition_CSm_R_5_1_2`, product-side S2–S3 separation → C(SH2) + SH | Two CAS(5,4) starts, 19,20,21,22 versus 19,20,21,23, from 30 Å; two doublets | Both starts reach the same far solution. Near geometry: virtual-22 DZ retains carbon frontier/SH character, but TZ contains a sulfur s-like and correlating orbital instead. Virtual-23 gives similar DZ/TZ near solutions, but one SH p component becomes inactive and an occupied C/S2 orbital enters the active window. | Neither comparison demonstrates the same intended subspace throughout. Stage only the five-orbital union described below. |
| `r12_insertion_R_2_1_5`, product-side C1–H4 separation → HCS + H2S | CAS(3,3), 19,21,22, versus donor-added CAS(5,4), 19,20,21,22; seeds from actual L2 complex; one doublet | CAS(3,3) active 20/21/22 retain HCS pi/radical/pi-star: ≈1.8995,1.0000,0.1005 at both endpoints and in both bases. In CAS(5,4), the H2S donor is active at 30 Å but becomes inactive near the complex; an HCS occupied valence orbital replaces it. | The corrected CAS(3,3) is a credible compact candidate. The larger result does not isolate a donor contribution. Compare a consistently defined HCS-only CAS(5,4). |

For the S-atom case, roots have weights 1/3 within each multiplicity block and
blocks have equal total weight 1/2. The CH+H2SS calculation uses one doublet block
with two roots of weight 1/2. These are optimization averages; they do not define
additional reactive flux. The large near-geometry root splittings and the
far-degenerate roots must be tracked in subsequent scans, rather than assuming
that a fixed root ordinal proves state continuity.

The HCS/H2S corrected CAS(3,3) natural pi orbital has valence C/S coefficients and
occupation about 1.8995; its pi-star partner is about 0.1005. The old near-2/1/0
occupations associated with sulfur core-like orbitals are not used for selection.
The new donor-added near solution has H2S donor character in inactive orbital 18,
while active 19 has HCS in-plane valence character. Thus its extra stabilization
cannot be called an H2S donor effect.

The numerical consequences illustrate why completion and lower energy alone
are insufficient. The following are **D1 CASSCF differences, E(near)-E(30 Å),
in Eh**, at fixed geometries, without counterpoise or dynamic correlation.
They are not production binding energies or rate predictions.

| Candidate | avdz-f12 | avtz-f12 |
| --- | ---: | ---: |
| HCS/H2S CAS(3,3) | +0.00145334 | +0.00156102 |
| HCS/H2S donor-start CAS(5,4), character changes between endpoints | -0.01243432 | -0.01204189 |
| C(SH2)/SH virtual-22 start, inconsistent near character across bases | +0.03334507 | +0.03730379 |
| C(SH2)/SH virtual-23 start, near/far active-space exchange | -0.00310556 | -0.00317441 |

The focused comparisons are:

- **C(SH2)+SH: CAS(7,5), source 18,19,20,21,23 from each 30 Å SCF checkpoint.**
  This adds the occupied C/S2 valence orbital 18 that entered the smaller near
  solution, while retaining the SH pair 19/21, carbon frontier 20 and candidate
  acceptor 23. Two equally weighted doublets, actual complex and 30 Å, both bases:
  four jobs. The added orbital has a specific observed purpose; no full-valence
  expansion or claim of a converged minimum space is made.
- **HCS+H2S: CAS(5,4), source 18,19,21,22 from each L2-complex SCF checkpoint.**
  It retains HCS in-plane valence 18, pi 19, radical 21 and pi-star 22, with the
  H2S donor inactive at both endpoints. One doublet, complex and 30 Å, both bases:
  four jobs. Compare its final character with the already completed CAS(3,3).

## Eight additional requests: source character and minimal hypotheses

Parent `791671450790040000002` has the same atom ordering. All selected indices
in this table refer to each request's **own 30 Å SCF checkpoint**, separately
bound to the printed DZ and TZ hashes. Similar products do not authorize using
another request's indices or checkpoint. Labels 20–22 within a three-orbital
span may exchange order without changing the intended complete span.

| Request suffix | Products and near coordinate | Compact trials; same source indices in DZ/TZ | States and reason |
| --- | --- | --- | --- |
| `hom_sci_1_2` | Triplet CH2 + doublet HSS; C1–S2, parent-derived 2 Å | CAS(3,3): **20,21,22** | One doublet and one quartet, equal total weights. Retain the two CH2 open-shell directions and HSS radical. |
| `r12_insertion_R_2_1_5` | CH2 + HSS; H5–S2, actual product complex | CAS(3,3): **20,21,22** | Same spin coupling, independently identified endpoint orbitals. |
| `r12_insertion_R_2_1_6` | CH2 + HSS; H6–S2, actual product complex | CAS(3,3): **20,21,22** | Same spin coupling, independently identified endpoint orbitals. |
| `hom_sci_2_3` | Singlet CH2S + doublet SH; S2–S3, parent-derived 2 Å | Two CAS(5,4) trials: **18,19,21,22** and **19,20,21,22** | Two equal doublets. Retain SH p pair 19/21; compare CH2S pi 18 versus sulfur donor 20, with candidate acceptor 22. |
| `r12_insertion_R_2_3_4` | CH2S + SH; H5–H4, actual product complex | Same two CAS(5,4) index sets, from this exit's own seeds | Two equal doublets; test pi versus donor participation along this separate outer coordinate. |
| `r12_insertion_R_5_1_6` | Singlet HCSSH + H5; H6–H5, actual product complex | Two CAS(3,2) trials: **19,21** and **20,21** | One doublet. H 1s SOMO 21 plus terminal-sulfur p donor 19 or mixed C/S occupied frontier 20. |
| `r12_insertion_R_6_1_2` | HCSSH + H6; H4–H6, actual product complex | **19,21** versus **20,21**, CAS(3,2) | One doublet; each seed is tied to this endpoint and atom map. |
| `r12_insertion_R_6_1_5` | HCSSH + H6; S3–H6, actual product complex | **19,21** versus **20,21**, CAS(3,2) | One doublet; each seed is tied to this endpoint and atom map. |

HCSSH abbreviates the requested singlet `[H][C]#[S]#[S][H]` connectivity; displayed
network bond orders do not themselves assign molecular orbital character.
For HCSSH+H, CAS(1,1) would retain only the H radical. The two CAS(3,2) hypotheses
also allow one occupied molecular donor to participate; they do not establish
that donor participation is necessary or that no other orbital matters. The
lowest virtual orbitals show substantial diffuse mixing and order changes
between bases, so a common virtual ordinal is not selected automatically.
For CH2S+SH, virtual 22 remains a candidate acceptor whose optimized character
must be checked, rather than a presumed pi-star assignment from orbital energy.

Each model has a near and 30 Å target in both bases. The three CH2+HSS requests
produce 12 jobs; the two CH2S+SH requests produce 16; the three HCSSH+H requests
produce 24. Together with the two focused comparisons, this is 60 jobs.
The recipes retain the existing SO-SCI/QN/exact-diagonalization settings and
convergence thresholds. They do not alter a production state model.

The six new outer-exit geometries are the normal L2 product endpoints 1763,
1780, 1784, 1802, 1811 and 1817 recorded in the
[previous endpoint evidence](../r7_compact_followup_20261025/additional_endpoint_evidence.json).
The two earlier outer exits use rows 843 and 1162. Their internal fragment
geometries remain fixed and only the specified interfragment coordinate is
translated to 30 Å. Their near points are the actual L2 complexes, **not a
forced 2 Å homolysis geometry**. The inner TSs are separate network objects.

## Blodgett handoff

Run over the user's SSH session after pulling the accompanying commit. These
are new review workers using the existing general driver; production and
RotdPy controllers need no restart for this panel. The panel cap remains three
within the existing shared production policy, whose requested cap is 23.
The commands preserve that policy and its exclusions, including node34 if
excluded. They do not increase the cap, override exclusions or cancel jobs.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
test "$(git branch --show-current)" = composite
git pull --ff-only origin composite
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PLAN="$PWD/examples/anl/orbital_review_batches/r7_compact_resolution_20261025/screens.json"
KB_PANEL="$KB_RUN/orbital_recon/compact_resolution_20261025"

"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_RUN" "$KB_PLAN" "$KB_PANEL" --shared-pool --preflight --resume
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_RUN" "$KB_PLAN" "$KB_PANEL" --shared-pool --resume
KB_LOG="$KB_PANEL/run.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.requested_orbital_recon \
  run "$KB_PANEL" --report-failures >"$KB_LOG" 2>&1 </dev/null &
printf 'Review PID: %s\nLog: %s\n' "$!" "$KB_LOG"
BASH
```

Preflight checks all source completions, indices, electron counts, required
SOMOs, native input/output/checkpoint hashes and request bindings before staging.
A changed checkpoint stops preparation; do not replace a saved hash merely to
bypass that check. Staging retains the previous 56 jobs and evidence unchanged.
The runner admits at most three review allocations subject to the shared cap;
if production occupies all capacity, review submissions wait. A repeated stage
checks the same plan; runner locking and saved job IDs prevent duplicate owners
and duplicate submissions. No new production model is published in this step.

Read-only progress and evidence capture:

```bash
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_PANEL="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/orbital_recon/compact_resolution_20261025"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor "$KB_PANEL"
KB_REPORT="$HOME/compact_resolution_character.$(date -u +%Y%m%dT%H%M%SZ).txt"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon character "$KB_PANEL" > "$KB_REPORT"
printf 'Character report: %s\n' "$KB_REPORT"
```

Review the final active character, occupations and roots at both endpoints and
in both bases before finalizing the remaining models together. Formal
subspace/state continuity and the correlated production potential are not
certified by this endpoint character report. In particular, `complete` means a
native calculation passed completion checks, not that its orbitals were accepted.

The previously published CH3S and CH3SS choices remain documented in
[ALL_ACTIVE_SPACES.md](../r7_progress_20261024/ALL_ACTIVE_SPACES.md). This report
adds selection evidence for the remaining requests; it does not re-approve
those earlier models or prove that the five-channel reconciliation has reached
production scans. No new production monitor was supplied with this character
report, so current production progress and the live pool setting are not
independently verified here.

## Validation scope

The new regression exercises the full 60-job plan using the supplied atom maps
and occupation records with synthetic checkpoints. It verifies all-panel
preflight, immutable existing work, source-hash binding, preservation of product
geometries, minimal electron/orbital counts and SOMOs, the executor's same-basis
restart, rejection of a bad last seed before any partial panel, and repeat
staging. A simulated scheduler begins with 23 production jobs, lets one finish,
retains the other 22, admits all 60 review jobs within the shared cap, reconciles
claims and resubmits none on repeat execution. No local Molpro or remote Slurm
execution of the new panel is claimed.
