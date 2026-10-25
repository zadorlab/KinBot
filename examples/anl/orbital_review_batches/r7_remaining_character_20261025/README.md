# Remaining VRC selections: character probes and compact-space review

This replaces the withdrawn `(19,12)/(19,15)` launch. It prepares **16 independent
SCF single points**, not a full-valence CI: four channels, two geometries and
two bases (`avdz-f12`, `avtz-f12`). The output is the orbital-character evidence
needed to specify compact multistate probes. No final active indices, optimized
multistate roots, production weights or new surface model are approved by this
stage. The existing selections in [ALL_ACTIVE_SPACES.md](../r7_progress_20261024/ALL_ACTIVE_SPACES.md)
are unchanged. The four new models will be published together after review.

## The two S–S coordinates describe different events

For parent `791871710370000000002`, atom order is C1 S2 S3 H4 H5 H6.
`hom_sci_2_3` is S2–S3 cleavage, with CHSH2 radical + atomic S products.
`R_Addition_CSm_R_5_1_2` has an inner rearrangement TS, followed by a
C(SH2) + SH complex and its outer separation. Its H5 has moved from C1 to S3.
Using S2/S3 as the outer separation coordinate does not make that rearrangement
another homolysis.

The supplied native/ASE evidence gives:

| Structure | DB row | C1–H5 / Å | S3–H5 / Å | S2–S3 / Å |
| --- | ---: | ---: | ---: | ---: |
| L2 parent | 8 | 1.095881 | 3.406387 | 1.959128 |
| L1 IRC forward endpoint | 510 | 1.114636 | 3.443372 | 1.993825 |
| L1 IRC reverse endpoint | 515 | 4.490312 | 1.360400 | 2.810540 |
| L2 reverse product complex | 843 | 3.906176 | 1.342667 | 2.485724 |

L2 electronic energies in hartree, converted using KinBot's eV conversion:
complex −835.7837550333901; C(SH2) −437.1349527923601; SH −398.63289985041007.
Their harmonic ZPVEs are 0.028528, 0.017408 and 0.006166 hartree, respectively.
Thus the electronic binding depth is 9.979 kcal/mol and the harmonic E0 depth
is **6.870 kcal/mol**. All 12 complex harmonic frequencies are positive; the
lowest is 115.9176 cm⁻¹. The logged 29.59 kcal/mol belongs to the earlier L1
assessment. This establishes a bound L2 product-complex endpoint under the
current harmonic model. It does not certify an outer potential or rate, or
prove a separate L2 IRC was calculated.

## Finite first panel

| Recipe | Products and one-based fragment maps | Coordinate and geometries | SCF initializer |
| --- | --- | --- | --- |
| `ch_h2ss_cs` | CH doublet [1,5]; H2SS singlet [2,3,4,6] | C1–S2; parent-derived 2 and 30 Å | Doublet ROHF |
| `chsh2_s_ss` | CHSH2 doublet [1,2,4,5,6]; S triplet [3] | S2–S3; parent-derived 2 and 30 Å | Quartet ROHF |
| `csh2_sh_exit` | C(SH2) singlet [1,2,4,6]; SH doublet [3,5] | S2–S3; actual L2 complex (2.485724 Å) and 30 Å | Doublet ROHF |
| `hcs_h2s_exit` | HCS doublet [1,3,5]; H2S singlet [2,4,6] | C1–H4; actual L2 complex and 30 Å | Doublet ROHF |

The exit geometries come from normal `_IRC_R_prod_high` records 843 and 1162,
bound to their exact accepted-network reaction and product-complex IDs. Preparing
the panel checks the atom order, both fragment connectivities, charge/spin
declarations and maps. Translation preserves each fragment's internal geometry.
The production correction coordinate is an interfragment atom-pair distance;
RotdPy dividing-surface/pivot distances are separate quantities.

Each job prints all occupied and eight virtual orbitals, with AO coefficients
down to 0.05; it retains the WFU and native XML orbital exports. The prepared
Molden directives are translated by the existing runner to XML for basis-format
compatibility. Printing controls follow the [Molpro documentation](https://www.molpro.com/manual/doku.php?id=general_hints_-_frequently_asked_questions).
AO coefficients are amplitudes, not populations. An SCF solution can have the
wrong separated-fragment state or charge distribution; its integer occupations
cannot validate a CAS or root average. The quartet initializer does not select
quartet reactive flux.

`orbital_character_only` prevents these SCF energies from being interpreted as
a correlated potential or used to request radial refinement. The panel does
not launch relaxed scans, CASPT2, CC checks or automatic active-space expansion.

## Compact probes to construct from the returned indices

These are **counted hypotheses for review**, not installed selections. Actual
indices must be assigned separately in each basis from the character output,
then confirmed by optimized compact CASSCF natural orbitals and CI character.
Source indices may change or mix with geometry; do not select a fixed energy
window or transfer another molecule's ordinal indices.

| Channel | First compact hypothesis and chemical contents | Initial state-average hypothesis |
| --- | --- | --- |
| CH + H2SS | CAS(5,4): CH's three frontier orbitals/three electrons plus the interacting sulfur donor/two electrons. Check whether H2SS requires an additional donor/antibond pair; that is an explicit comparison, not assumed full valence. | Two doublets, 1/2 each, for the CH orbital pair. A four-doublet comparison is justified only if nearby states or the donor interaction mix with that pair. |
| CHSH2 + S | CAS(5,4): atomic S 3p set/four electrons plus the fragment radical orbital/one electron. Test CAS(7,6) only if an identified delocalization/bond pair is missing. | Three doublets and three quartets, equal within each spin block, representing the S spatial manifold coupled to the radical. Spin-block balancing and root identities must be inspected. |
| C(SH2) + SH exit | CAS(5,4): SH frontier pi pair/three electrons plus the singlet fragment's occupied/acceptor frontier pair/two electrons. Confirm the singlet fragment's actual orbital character, rather than treating it as CHSH2. | Two doublets, 1/2 each, for SH's orbital pair. No atomic-S quartet manifold for this singlet + doublet asymptote. |
| HCS + H2S exit | Compare the HCS radical plus its C–S pi/pi-star pair, CAS(3,3), with CAS(5,4) adding the interacting H2S donor. This tests whether the donor must participate in the noncovalent approach. | One ground doublet initially. An additional averaged doublet requires evidence of a nearby interacting state; HCS is not assigned CH's orbital degeneracy. |

The fragment-term reasoning is supported by [HCS ground-state calculations](https://pubs.rsc.org/en/content/articlelanding/1998/ft/a800926k),
[CH/SH/HCS reaction spectroscopy](https://pmc.ncbi.nlm.nih.gov/articles/PMC7502777/)
and [atomic sulfur data](https://physics.nist.gov/cgi-bin/Elements/elInfo.pl?context=frames&element=16).
Those sources establish fragment terms, not the adequacy of these proposed
active spaces for the present potentials.

For each compact candidate, inspect near/far character, then the transition
region (3.5/6 Å for homolyses; an intermediate intact-fragment separation for
exits) and the same fixed geometries in both bases. Add the smallest specific
pair or extra roots only where the native result warrants it. Compare interaction
energies relative to each candidate's own asymptote if an expansion is tested;
lower absolute energy alone does not select the larger space. Equal state-average
weights within a degenerate manifold are an orbital-optimization choice, not
reactive flux weights. Ground-parent connectivity and the resulting scan determine
which adiabatic components contribute to the kinetic model.

There are no executable compact-candidate recipes with guessed indices in this
directory. SCF evidence cannot be published through the natural-orbital selection
publisher. After reviewing these outputs, prepare the finite compact calculations,
review all four, and publish them in one `review-batch --require-all` transaction.
Production then owns the relaxed/frozen scans, high-level corrections and RotdPy;
native scan-review gates remain active.

## Code changes accompanying the panel

* An optional exact `reaction` in each surface channel distinguishes two channels
  sharing a bond. It is carried through selection, additive publication, live
  inheritance, review acknowledgement, scan dispatch and reference-probe lookup.
  A bond-only channel remains compatible if unambiguous. It cannot be mixed with
  reaction-specific channels on the same bond or silently migrated.
* Both scanned fragments must match accepted products in connectivity,
  stereochemistry, charge and spin before fragment QC. A mismatch holds that
  channel; other reviewed channels proceed. Installing orbitals cannot clear an
  unresolved fragment-identity hold.
* Product-exit panels retain the exact accepted L2 endpoint, atom maps and
  network provenance in `geometry_source.json`; the parent request geometry is
  not used as the exit geometry. The log no longer calls product exits “homolytic”.
* The existing accepted-model inheritance fix and `resubmit --reconcile-models`
  remain part of the update. The five previously published CH3SS channels retain
  their existing CAS definitions and native seeds. A case-level receipt alone
  does not verify worker inheritance.
* `batch --preflight` checks all requests/decks in temporary storage before any
  campaign is created. `--resume` verifies the saved plan, recipes and ownership
  before reusing it. `character` checks saved output/WFU/XML hashes and reports
  orbital character without querying Slurm or altering results.

## Blodgett: load fixes, reconcile the five, and run the character panel

Execute in your SSH session. The generic handoff runs the full panel preflight
before replacing controllers. It then waits for the replacement PES/worker
owners and any current input publication, requests reconciliation through those
owners, prepares the new panel, and runs it. A blocked scientific revision is
reported rather than overridden.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
test "$(git branch --show-current)" = composite
git pull --ff-only origin composite
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
KB_PANEL="$KB_RUN/orbital_recon/remaining_character_20261025"
KB_PLAN="$PWD/examples/anl/orbital_review_batches/r7_remaining_character_20261025/screens.json"
KB_LOG="$KB_RUN/vrc_reconcile_character.$(date -u +%Y%m%dT%H%M%SZ).log"
nohup "$PWD/.venv/bin/python" -B -u -m kinbot.orbital_review_handoff \
  "$KB_RUN" "$KB_PLAN" "$KB_PANEL" \
  --case-input "$KB_RUN/ch3ss/input.json" --expected-max-jobs 23 \
  --apply --wait-idle 600 --wait-owner 120 \
  >"$KB_LOG" 2>&1 </dev/null &
printf 'Handoff PID: %s\nLog: %s\n' "$!" "$KB_LOG"
BASH
```

This loads the new installed code into both production PES controllers **now**,
retains running/pending Slurm work and independent RotdPy owners, and reconciles
the five previously reviewed CH3SS channels. It does not wait for the four new
selections before releasing eligible old channels. A repeated same-code refresh
is a no-op; panel preparation is resumable and duplicate handoff/review owners
are locked out. No `--replace-pending` or cancellation option is used.

The campaign admits at most three review allocations, all **inside the existing
23-allocation configured production pool**, using its existing exclusions and
admission policy. It waits when the pool is full. The handoff checks the configured
limit is 23 and does not edit it. Existing reservation settings, if present, are
also retained; the monitor reports the effective limit. It changes none of the
five accepted CAS definitions. Readiness waits are bounded; if one expires, rerun
the same block to resume from the saved refresh/panel state.

Monitor and return the native character report:

```bash
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_PANEL="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/orbital_recon/remaining_character_20261025"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor "$KB_PANEL"
"$PWD/.venv/bin/python" -B -m kinbot.requested_orbital_recon character "$KB_PANEL" > "$KB_PANEL/character.txt"
cat "$KB_PANEL/character.txt"
```

Check the reconciliation separately from the panel:

```bash
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_RUN="$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7"
"$PWD/.venv/bin/python" -B -m kinbot.run_control resubmit "$KB_RUN/ch3ss/input.json" --status
"$PWD/.venv/bin/python" -B -m kinbot.run_monitor "$KB_RUN"
"$PWD/.venv/bin/python" -B - <<'PY'
from pathlib import Path
import json
from kinbot.input_refresh import _model_map
r=Path.home()/'kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/ch3ss'
a=json.loads((r/'.kinbot_accepted_input.json').read_text())
models=_model_map(a,r)
for name,model in models.items():
 w=r/name;p=w/(name+'.json')
 if not p.exists():
  print('WORKER',name,'input not yet written');continue
 b=json.loads(p.read_text());overrides=w/'.kinbot_search_options.json'
 if overrides.exists():b.update(json.loads(overrides.read_text()))
 print('WORKER',name,'MODEL_MATCH',_model_map(b,w).get(name)==model)
 for p in sorted((w/'vrctst/orbital_review').glob('*/request.json')):
  q=json.loads(p.read_text())
  print('REQUEST',q.get('reaction'),q.get('status'),q.get('error') if q.get('status')=='needs_review' else '')
PY
```

`MODEL_MATCH True` verifies the effective worker input; `model_accepted` verifies
its acknowledgement. The monitor must then show fragment/scan work progressing
or its specific hold/queue reason. Neither the input receipt nor `model_accepted`
is proof that a native scan has started or completed. A remaining native scan or
reference-review hold is not cleared merely by restoring the orbital model.

After the compact results are reviewed, publish the four exact-reaction model
proposals together through `review-batch --require-all --apply`. That uses the
newly loaded live-control code and needs **no further controller replacement**.
This handoff submits only the character panel, not that future publication.

No new CAS is claimed accepted, no old review is relaunched, and no network
pruning, ANL downgrade, cross-PES sharing activation or QC cancellation occurs.
Those broader decisions remain deferred until both L1/L2 explorations finish.
Local regression checks cover orchestration and failure guards; they do not
constitute successful Molpro or Blodgett execution.

Validation for this update: **497 passed, 1 skipped** across the reconnaissance,
reaction identity, publication/inheritance, full controller refresh, shared
admission, submission recovery and handoff suites. The handoff tests exercise
startup waits, blocked/concurrent revisions, full-panel preflight before refresh,
duplicate-owner refusal and repeated preparation without changing existing work.
Refresh integration tests retain independent sampling owners and live allocation
IDs. These are local tests with simulated scheduler/native execution, not a
verified remote launch or successful native retry.
