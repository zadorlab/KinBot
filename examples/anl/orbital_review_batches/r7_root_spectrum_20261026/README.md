# Review update, 26 October 2026

**Completion update:** all twenty jobs completed and the transferred archive
verified. See the [native findings and next compact panel](../r7_state_support_20261026/README.md).
The launch below is historical and need not be repeated.

Twenty diagnostic jobs are prepared: **four fixed-orbital CASCI jobs** for the
unresolved C(SH2)+SH exit and **sixteen SCF character jobs** for five newly
requested channels. The two sulfur-shift requests share exactly the same
supplied endpoint and therefore one four-job character screen. All twenty jobs
use one review controller, capped at three and admitted through the existing
production pool. With 22 allocations against cap 23 in the supplied snapshot,
only one additional allocation fits initially. Actual admission uses live counts.

The previous fourteen jobs are complete and archived; do not repeat them.
Eleven production selections remain prepared, and the combined publication is
still held. This update gathers the remaining scientific and failure evidence.
It does not publish models or restart PES controllers.

## Last channel: what the completed natural-start calculation established

The [reproducible native audit](natural_union_audit.json) reads the complete
outputs and XML from the verified export. Both near calculations read the
intended natural record 2140.2, but CASSCF again replaced its SH-p direction
with an occupied valence direction. Against the intended six-dimensional
natural 17–22 span, the minimum overlap singular value is **0.121846 DZ /
0.122409 TZ**. Replacing active 17 with inactive SH 16 raises this to
**0.990473 / 0.990240**. Near energies reproduce the previous canonical-start
branch within approximately 2e-8 Eh. Repeating that initialization is exhausted.

At 30 Å the continued D1 has dominant three-open-shell CI strings
`222aba`, `222baa`, `222aab`, while D2 is dominated by `2222a0`. Together with
the AO character this suggests competition between carbon-fragment electronic
states. It does **not** establish isolated-fragment ground-state ordering or
justify replacing the requested singlet-fragment exit with a triplet channel.

| Continued natural-start result / Eh | DZ | TZ |
| --- | --- | --- |
| Near D1 | -834.477054809332 | -834.491571890347 |
| Near D2 | -834.441989877345 | -834.457632584248 |
| 30 Å D1 | -834.481898752232 | -834.496164694357 |
| 30 Å D2 | -834.470943757414 | -834.486039847537 |

The next diagnostic uses **CAS(9,6), natural 17–22**, at each of the unchanged
near/far geometries in each basis. Each starts independently from the matching
completed CAS(7,5) natural checkpoint, opening its inactive 17. Near, that adds
the missing SH-p direction; at 30 Å it adds the carbon-p direction outside the
previous five-orbital space. These are geometry-specific identities, not an
assumption that ordinal 17 always has the same character.

Four doublets and two quartets provide a diagnostic window for the two SH
orbital components coupled to singlet or triplet carbon-fragment character.
Higher roots may still enter that window. Equal weight over the six states
defines the reported diagnostic natural density; it does not select production
averaging weights or reactive flux. Near/far CI character and the quartet
comparison must be examined before choosing any subsequent optimization.

`CAS-CI` keeps orbital subspaces fixed; `NOEXTRA` disables implicit
extra-symmetry adaptation. Natural orbitals can rotate within those fixed
subspaces after the CI solve. The extraction compares the resulting active span
with its source, so loss of a direction remains detectable. This follows the
[Molpro CAS-CI definition](https://www.molpro.net/manual/doku.php?id=the_mcscf_program_multi#cas-ci).
Shared code checks checkpoint hashes, unchanged source geometry, independent
initialization and absence of CASSCF optimization options. These seeded CASCI
diagnostics cannot be published as an optimized production selection.

## Five new requests: four distinct character screens

Parent: `791871560500020000002-sb6e34ded5baaf03b`. Atom order is
**C1, S2, S3, H4, H5, H6**. Definitions and exact coordinates are transcribed
in [new_requests_evidence.json](new_requests_evidence.json), then checked against
the live request, accepted network and normal L2 database record before launch.

| Reaction suffix | Coordinate and fragment maps | Character to identify |
| --- | --- | --- |
| `hom_sci_1_2` | C1–S2 cleavage; CH=(1,6), HSSH=(2,3,4,5); 2 and 30 Å | CH sigma/frontier-pair directions and the participating sulfur donor/bonding directions |
| `hom_sci_2_3` | S2–S3 cleavage; HCSH=(1,2,4,6), SH=(3,5); 2 and 30 Å | SH p pair and the HCSH frontier directions participating in S–S association |
| `12_shift_S_F_3_2_4` | HCS=(1,2,6), H2S=(3,4,5); actual L2 endpoint and 30 Å along C1–H4 separation | HCS radical/pi directions and H2S donor character; outer fragment separation, not parent-bond cleavage |
| `12_shift_S_R_3_2_4` | Same supplied endpoint and fragment maps, distinct reaction record | Shares the preceding screen only after exact live geometry/map checks |
| `r12_insertion_R_4_2_1` | HCSSH=(1,2,3,5,6), H4; actual L2 endpoint and 30 Å along H5–H4 separation | H-atom SOMO and any fragment frontier mixing; no H–H bond is being dissociated |

Every distinct geometry is computed in avdz-f12 and avtz-f12. Only rigid
fragment translation generates the far points; internal geometry is preserved.
The three non-homolysis requests use their **product-side L2 geometries**.
SCF-only jobs print occupied and virtual orbitals and export complete XML/WFU.
They have no CAS, CI root average or correlated calculation. Native character
will determine compact candidate indices; no indices are borrowed from another
well. Matching chemistry can guide the candidate choice, but cannot establish
checkpoint identity or optimized orbital quality.

The F/R shift endpoint records are 1318/1321 in the supplied database; the
H-atom exit is record 1295. Live row IDs are recorded anew. If either shared
endpoint changes, the check refuses sharing. Separate reaction identities and
subsequent flux accounting remain intact; this is not a duplicate-path verdict.

## Blodgett commands

From the user's SSH session, pull the accompanying commit and run:

```bash
cd "$HOME/KinBot"
git pull --ff-only origin composite
bash examples/anl/orbital_review_batches/r7_root_spectrum_20261026/run_blodgett.sh
```

The script validates the five requests and endpoint sharing, preflights all
twenty decks, stages the campaign with `--resume`, captures current native
failure/request evidence, and starts only its review controller. Existing
allocations, limits, exclusions, operator decisions and sampling are preserved.
Repeated invocation uses the existing campaign and lock; completed work is not
submitted again. A changed scientific recipe is refused on resume.

Monitor without changing state:

```bash
"$HOME/KinBot/.venv/bin/python" -B -m kinbot.requested_orbital_recon monitor \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/orbital_recon/root_spectrum_20261026"
```

The launch also writes `production_status_before.json` inside that campaign.
It reads pending request definitions, both L2 endpoint directions, bounded TS
and RotdPy failure tails, controls and controller metadata. SQLite is opened
read-only. Missing or partial files are reported without altering calculations.
The capture is per-file, not an atomic whole-run snapshot. This file can be
transferred immediately while the diagnostics run.

After completion, collect the text, native files, WFU/XML and current status:

```bash
bash "$HOME/KinBot/examples/anl/orbital_review_batches/r7_root_spectrum_20261026/collect_blodgett.sh"
```

Transfer the printed export directory into local `examples/anl/tools/`.
The exporter verifies archive integrity; it does not approve the science.
The older compact-continuation and natural-union archives are already verified
locally and need not be exported again. Their hashes are in [the audit](status_audit.md).

## Publication and controller handoff still required

The original twelve-request set has eleven prepared selections and the last
state model unresolved. The five new requests also need compact selections
after this character stage. The final combined plan must account for all of
them, including the shared endpoint's separate reaction binding, and reconcile
the five previously published channels. The existing held plan cannot partially
publish. No additional channel is rejected by this update.

The active-space explanations remain in [ALL_ACTIVE_SPACES.md](../r7_progress_20261024/ALL_ACTIVE_SPACES.md)
for existing/earlier selections and the [compact-selection audit](../r7_compact_union_20261025/README.md#disposition-of-the-other-ten-requests)
plus [CH2/HSS result](../r7_compact_union_20261025/union_audit.md) for the eleven
new prepared choices. Their source indices, packed natural indices, averaged
states, reactive flux and validation limits remain distinct.

The final refresh must load endpoint/publication code into owning PES/KB
controllers and retain live allocations and RotdPy campaigns. The audited
endpoint defect calls for VRC preparation retries; it indicates no earlier
PES search, IRC, L2, ANL or CBH rerun. Remote dispatch/adoption must be checked
after refresh, beyond the receipt. The [ANL/MESS thermal and live-control
contract](../../../../docs/anl_live_controls.md) remains unchanged; no new native
thermal or final MESS validation is claimed here.

Local validation: **350 tests passed, one platform-specific test skipped**
(197 orbital/evidence tests and 153 orchestration/recovery tests). Tests cover
the exact twenty-job plan, endpoint identity/sharing, rigid
fragment motion, hash/geometry refusals, publication guards, shared-cap
admission, retention of 22 live jobs and repeated invocation. Refresh/model
inheritance/vdW/RotdPy orchestration checks also pass. These are simulated-site
software checks; the new native diagnostics and remote handoff have not run.
