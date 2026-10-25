# Compact follow-up and eight additional orbital requests

This panel uses the **16 completed SCF character probes** in
`orbital_recon/remaining_character_20261025` and the eight additional requests
under parent `791671450790040000002`. It stages **56 fixed-geometry jobs**:
24 compact CASSCF probes and 32 SCF character probes. All use `avdz-f12` and
`avtz-f12`. The panel cap is three allocations within the existing shared
production pool; it does not increase the pool limit or change exclusions.

The recorded SCF character and WFU/XML hashes are preserved in
[scf_character_evidence.txt](scf_character_evidence.txt). The six newly supplied
product endpoints are preserved in
[additional_endpoint_evidence.json](additional_endpoint_evidence.json).
Preparation reads the live requests, accepted networks, database endpoints and
source checkpoints on Blodgett. These local records document the reasoning and
provide regression fixtures; they do not override live evidence.

## Compact candidates for parent 791871710370000000002

Atom order is C1 S2 S3 H4 H5 H6. Indices below are **source canonical SCF indices**,
record `2100.2`, separately checked in both bases. They are packed into the
compact active window once. Optimized natural-orbital indices will differ.
Every explicit seed is bound to the WFU SHA-256 printed in the supplied report;
a newly recomputed SCF solution cannot silently inherit the same index choice.

| Channel | Trial space and source indices in DZ and TZ | Chemical hypothesis | Trial state average |
| --- | --- | --- | --- |
| `hom_sci_1_2`: CH + H2SS | CAS(5,4), **17,18,21,22**, 30 Å seed | S2 donor (17, two electrons); CH occupied sigma frontier (18, two electrons) and p frontier pair (21/22, one electron). Include the donor to permit interaction with the CH subspace. | Two doublets, weights 1/2 each |
| `hom_sci_2_3`: CHSH2 + S | CAS(5,4), **18,20,21,22**, 30 Å quartet seed | All three atomic-S 3p directions (18/20/21, four electrons) plus the fragment radical (22, one electron). This includes the S spatial manifold and its spin couplings. | Three doublets and three quartets, equal weights within each block and equal total spin-block weights |
| `R_Addition_CSm_R_5_1_2`: C(SH2) + SH outer exit | Two CAS(5,4) guesses: **19,20,21,22** and **19,20,21,23**, 30 Å seeds | SH pi pair (19/21, three electrons), occupied C(SH2) frontier (20, two electrons), and a candidate acceptor. Printed amplitudes do not uniquely resolve virtual 22 versus 23; compare the optimized subspaces. | Two doublets, weights 1/2 each |
| `r12_insertion_R_2_1_5`: HCS + H2S outer exit | CAS(3,3), **19,21,22**, compared with CAS(5,4), **19,20,21,22**, actual L2-complex seeds | HCS C–S pi/pi-star pair and radical; the larger trial adds the H2S donor (20). This tests whether donor participation is needed for the approach. | One ground doublet |

Each candidate has two target geometries in both bases: parent-derived 2/30 Å
for homolyses, actual L2 complex/30 Å for exits. The latter use the previously
verified product rows 843 and 1162. These are six candidate models and 24 jobs.

The far-separation seeds give interpretable fragment orbitals. The near CH+H2SS
SCF virtual 22 is diffuse/mixed, so the near ordinal is not treated as the far
CH frontier orbital. For C(SH2)+SH, the singly occupied SCF orbital is localized
primarily on different fragments at the two geometries. That observation does
not prove physical charge transfer or state correspondence. The two acceptor
guesses deliberately test a remaining ambiguity within the same compact size.

A seed transferred between geometries is only a starting guess. Convergence must
be followed by inspection of optimized natural occupations, AO/CI character,
fragment charge/spin and state continuity. If the compact near and far solutions
are consistent, intermediate points and correlated interaction energies can
then resolve any remaining ambiguity. Absolute energies of different CAS sizes
alone do not select a model. Equal optimization weights do not specify reactive
flux, and the quartet probes do not add quartet flux to production.

These jobs optimize CASSCF orbitals only. They do not launch CASPT2, relaxed
scans, full-valence CI or automatic active-space expansion. No new production
space or kinetic root weight is accepted by this handoff.

## Eight additional requests: parent 791671450790040000002

Atom order remains C1 S2 S3 H4 H5 H6. The native parent is the supplied
`[H][S]=[S]=[C]([H])[H]` structure. The atom pair in an exit row is an
**interfragment separation coordinate**, not a claim that the corresponding
covalent bond breaks in the inner reaction.

| Reaction suffix | Products, spin and one-based fragment maps | Separation pair | Near geometry | SCF initialization |
| --- | --- | --- | --- | --- |
| `hom_sci_1_2` | CH2 triplet [1,5,6] + HSS doublet [2,3,4] | C1–S2 | Parent-derived 2 Å | Quartet ROHF |
| `hom_sci_2_3` | CH2S singlet [1,2,5,6] + SH doublet [3,4] | S2–S3 | Parent-derived 2 Å | Doublet ROHF |
| `r12_insertion_R_2_1_5` | CH2 triplet [1,5,6] + HSS doublet [2,3,4] | H5–S2 | Normal L2 reverse product row 1763 | Quartet ROHF |
| `r12_insertion_R_2_1_6` | CH2 triplet [1,5,6] + HSS doublet [2,3,4] | H6–S2 | Normal L2 reverse product row 1780 | Quartet ROHF |
| `r12_insertion_R_2_3_4` | CH2S singlet [1,2,5,6] + SH doublet [3,4] | H5–H4 | Normal L2 forward product row 1784 | Doublet ROHF |
| `r12_insertion_R_5_1_6` | HCSSH singlet [1,2,3,4,6] + H doublet [5] | H6–H5 | Normal L2 forward product row 1802 | Doublet ROHF |
| `r12_insertion_R_6_1_2` | HCSSH singlet [1,2,3,4,5] + H doublet [6] | H4–H6 | Normal L2 reverse product row 1811 | Doublet ROHF |
| `r12_insertion_R_6_1_5` | HCSSH singlet [1,2,3,4,5] + H doublet [6] | S3–H6 | Normal L2 forward product row 1817 | Doublet ROHF |

Here HCSSH abbreviates the requested singlet `[H][C]#[S]#[S][H]` connectivity;
its displayed bond orders are the network representation, not an orbital
assignment. Each request gets its near geometry and a 30 Å translation in both
bases: **32 SCF jobs**. The six exits retain the internal L2 fragment geometries.
Preparation checks the exact reaction, product-complex record, normal endpoint
status, atom order, fragment partition and fragment connectivity. Similar
products and reused product-complex IDs do not merge distinct requests or
justify transferring orbital indices.

These screens print occupied and eight virtual orbitals and retain native WFU
and XML files. Their role is to identify the radical/frontier/donor subspaces:
CH2+HSS needs the two CH2 open-shell orbitals and the HSS radical character;
CH2S+SH requires the SH frontier pair and the interacting CH2S orbitals to be
identified; HCSSH+H requires the H 1s orbital and the interacting molecular
frontier/donor character to be identified. Counts, indices and state averages
will be chosen from that evidence. No full-valence space is inferred.

## Staging and execution

Run on Blodgett after pulling the accompanying code. The user owns SSH.
The existing production and RotdPy controllers keep running; these new review
workers load the pulled code. No controller refresh or production publication
is part of this panel.

```bash
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PLAN="$PWD/examples/anl/orbital_review_batches/r7_compact_followup_20261025/screens.json"
KB_PANEL="$KB_RUN/orbital_recon/compact_followup_20261025"

"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_RUN" "$KB_PLAN" "$KB_PANEL" --shared-pool --preflight --resume
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  batch "$KB_RUN" "$KB_PLAN" "$KB_PANEL" --shared-pool --resume
KB_LOG="$KB_PANEL/run.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.requested_orbital_recon \
  run "$KB_PANEL" --report-failures >"$KB_LOG" 2>&1 </dev/null &
printf 'Review PID: %s\nLog: %s\n' "$!" "$KB_LOG"
```

The first command validates the complete panel without submitting or changing
production files. The second stages it, or validates an identical existing
panel. A late preparation failure rolls back only newly created, unsubmitted
staging. The runner resumes journaled jobs and completion records; locks prevent
a duplicate owner. It respects the shared cap/exclusions at each admission,
retains existing jobs and uses the existing bounded recovery procedures.
A full 23-job pool means new review jobs wait for shared capacity.

```bash
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor "$KB_PANEL"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon \
  character "$KB_PANEL" > "$KB_PANEL/character.txt"
```

Native output/input, WFU/XML, hashes, execution records and immutable seed
provenance remain under the panel for review and later evidence export.
`complete_review_required` is computational completion, not publication.

The five previously published models and their worker reconciliation remain
separate from these twelve requests. The reconciliation code is included in
ancestor `e89b8dbb`; the last inspected worker snapshot had not yet demonstrated
all five channels passing their ready boundary. A case-level refresh receipt
alone does not establish that their scans started. Existing published selections
are documented in [ALL_ACTIVE_SPACES.md](../r7_progress_20261024/ALL_ACTIVE_SPACES.md).

## Implementation and verification

The shared recon code accepts explicitly selected, evidence-bound SCF seeds,
checks electron counts and all singly occupied orbitals, and keeps the source
canonical record on FILE,3 while new RHF/CASSCF records are written on FILE,2.
The index exchange uses Molpro's documented zero-angle `ROTATE` after explicit
`START`; see the [MULTI manual](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi).
Neither the new code nor its recovery paths branch on these species or paths.

Local tests cover the complete mixed panel, actual supplied fragment maps,
input/source hash checks, late preparation failure and retry, repeated
invocation, separate source/destination records, and a simulated scheduler
retaining production jobs under the shared cap. Those tests do not substitute
for native Molpro convergence or a verified Blodgett launch.

## Recovery after the 25 October native results

The 12:23 PDT snapshot contained four VIRTORB failures during forced-OVL RHF
initialization, two native CASSCF completions rejected by the occupation reader,
nine completed SCF character jobs, one pending job, and forty unsubmitted jobs.
This is a timestamped snapshot, not a claim about the current queue.

The four failed initializations use far-separation SCF seeds. Their native
inputs forced spatial overlap projection despite matching bases. The supplied
CH+H2SS and CHSH2+S outputs report many occupied overlaps of zero, followed by
inconsistent virtual counts. The shared executor now uses a normal same-basis
coefficient restart from the immutable source record. The owning review driver
verifies provenance, archives these attempts, and retries each once.

The two HCS+H2S L2-geometry jobs terminated normally:

| Job | Candidate | Native energy / Eh | Native elapsed / s | Final active occupations |
| --- | --- | ---: | ---: | --- |
| 52931214 | CAS(3,3), source 19,21,22 | -834.641533617762 | 10.48 | 1.99965, 1.00000, 0.00035 |
| 52931217 | CAS(5,4), source 19,20,21,22 | -834.641785308293 | 23.40 | 1.99999, 1.99965, 1.00001, 0.00036 |

The reader previously recognized only `NATURAL ORBITALS (state averaged)`.
It now also reads the bare `NATURAL ORBITALS` heading for these single-root
calculations. This corrects computational status; it does not validate their
active spaces. Their final active window contains sulfur core-like p character
at approximately -6 Eh and a high-energy correlating orbital around +13 Eh.
The intended valence pi/pi-star orbitals are outside that final active window.
These two results are therefore **not suitable selection evidence for the
intended valence CAS**.

The initial MULTI excerpts confirm record 2100.3, occupied overlaps of 1.000000,
and the expected 20.1/19.1 exchange for CAS(3,3). They report overlaps only through
occupied orbital 21, so they do not establish preservation of the selected
virtual orbital 22. They also do not distinguish a virtual initialization
problem from subsequent active-space drift. The next controlled check repeats
only these two completed jobs with the same source indices and scientific
settings, using the corrected restart. Original native files are archived, not
discarded. The six unsubmitted HCS/H2S points also use the corrected executor.
No additional CAS, extra root, frozen core or CI initialization option is added.
Any remaining character failure requires inspecting the new native evidence.

Run the following on Blodgett after pulling the repair commit supplied in the
handoff. This resumes the existing 56-job panel; it does not stage a duplicate
panel or restart production/RotdPy controllers.

```bash
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PANEL="$KB_RUN/orbital_recon/compact_followup_20261025"
KB_LOG="$KB_PANEL/seed_recovery.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.requested_orbital_recon \
  run "$KB_PANEL" --report-failures --wait-owner 86400 \
  --recheck-seed-job 52931214 --recheck-seed-job 52931217 \
  >"$KB_LOG" 2>&1 </dev/null &
printf 'Review PID: %s\nLog: %s\n' "$!" "$KB_LOG"
```

If the previous owner is still running its independent SCF jobs, this process
waits for its lock. It then rechecks the saved panel, archives the two explicitly
selected completions, recovers only the proved VIRTORB failures, and resumes
submission. The panel cap remains three within the existing shared pool of 23.
The shared node exclusions remain effective and now appear in the orbital monitor.
Repeated invocation recognizes archived old IDs and preserves newer attempts.

Read-only follow-up:

```bash
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor "$KB_PANEL"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon character "$KB_PANEL"
```

Inspect both the repaired near probes and their far/basis counterparts before
choosing a production model. The eight additional SCF screens keep their completed
results and continue. The existing five-model production reconciliation is a
separate status check; this review repair does not certify its remote outcome.
