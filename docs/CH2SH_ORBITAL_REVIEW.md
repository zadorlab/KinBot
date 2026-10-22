# CH2SH C–S orbital selection and live refresh, 2026-10-21

The reviewed r7 request is `471030810040000000002_hom_sci_1_2`: neutral
doublet CH2SH, atom 1 C and atom 2 S, dissociating to triplet CH2 and doublet
SH. This differs from the initial CH3S -> CH3 + S channel.

**Select CAS(5,4). Do not resume the 22-job validation panel.** The user
requested orbital assignment using the completed evidence, with production
KinBot responsible for the subsequent potential and spin-splitting energies.
The optional CC calculation was unnecessary for assigning these orbitals.
Four candidate calculations have native results; eighteen remain unsubmitted.
This handoff neither submits them nor treats the entire panel as validated.

## Chemical selection and evidence

The compact space contains the two CH2 frontier orbitals with two electrons
and the two SH pi orbitals with three electrons. The forming C–S bond and
antibond are combinations within that four-orbital span. Adding a separate
sigma/sigma-star pair would count those degrees of freedom twice.

The eight completed fixed-geometry CASCI(13,11) probes cover 2, 3.5, 6 and
30 angstrom in avdz-f12 and avtz-f12. Their natural-orbital AO coefficients
identify source orbitals 11–14 as the intended frontier space in both bases.
At 30 angstrom, 11/12 are principally sulfur pi with occupations about
1.4997 each; 13/14 are carbon-centered frontier orbitals, about 1.206/0.794.
Both carbon orbitals remain active when their order/character exchanges
between near and far geometries. The deeper sulfur 3s and C–H bonds are
inactive. Orbitals 10/15 are the additional S–H sigma/sigma-star pair tested
in CAS(7,6).

The completed optimized-orbital comparison supplied by the user is:

| Calculation | Occupations of the compact four orbitals | Added S–H sigma / sigma-star |
|---|---|---|
| CASSCF(5,4), aDZ, 2 Å | 1.70808, 1.48606, 1.02341, 0.78245 | Not active |
| CASSCF(7,6), aDZ, 2 Å | 1.71600, 1.47410, 1.02155, 0.78185 | 1.97595 / 0.03055 |
| CASSCF(7,6), aDZ, 3.5 Å | 1.51514, 1.50002, 1.18344, 0.80096 | 1.97625 / 0.02419 |
| CASSCF(7,6), aDZ, 30 Å | 1.49978, 1.49978, 1.20253, 0.79747 | 1.97629 / 0.02414 |

Compact CASSCF AO coefficients retain the two sulfur-dominated and two
carbon/bonding frontier combinations. The expansion leaves the added pair
near 2/0 without qualitatively changing the frontier occupations. This
supports the selected orbital space. It does not establish complete
two-basis correlated-potential agreement, angular accuracy or a rate-error
bound; those are not claimed as completed reconnaissance.

The recipe uses the original complete CASCI screen's **separate DZ and TZ**
2 Å WFU records, `2140.2,type=natural`, indices 11–14. Production CASSCF
optimizes this selected space (`CLOSED,10; OCC,14`, 25 total electrons).
It does not transplant a DZ checkpoint into TZ or require a new TZ probe.
Keep `r7_ch2sh_cs_orbital_recon`: its hashed checkpoints remain runtime seeds.

State averaging retains four doublets and two quartets from the diagnostic
screen. Weights are 1/4 per doublet and 1/2 per quartet before normalization,
so each spin block has half the total orbital-optimization weight. The extra
doublets allow inspection of other configurations, including possible
singlet-CH2-derived character; no such character is assigned solely by energy
rank. Only doublet root 1 is a reactive VRC surface. Quartet root 1 supplies
the MR reference side of its spin gap. The averaging weights are neither
Boltzmann populations nor reactive-flux weights.

## Parser correction

The compact aDZ native calculation terminated normally. Its optional CC
check printed both F12a and F12b summaries, while the reconnaissance parser
incorrectly required a single unscaled triples contribution in the whole
output. The parser now pairs base, triples and total energies within each
CC program and verifies the recorded selector against **F12b**. It still
rejects F12a selection, scaled triples, missing/ambiguous summaries and
inconsistent totals. This changes no native energy or production method.

The native F12b total was `-437.333094348795` Eh; the saved marker was
`-437.33309435` Eh. The tests retain the user's printed base/total values;
their synthetic triples lines are explicitly labeled as total-minus-base
format fixtures. No native re-run was performed for this repair.

## Apply the reviewed selection to the running driver

`review` is a generic recipe-driven operation. It checks request identity,
geometry provenance, complete basis-matched seed records and WFU/output
hashes; preserves existing channels; validates the new model; and requires
a live PES owner. `--apply` archives the input, atomically adds the model to
the same input file, and uses the existing journaled resubmission path.
Concurrent helper invocations cannot edit the input together. Reapplication
reconciles the same revision, including after the request is marked accepted.

```bash
bash <<'BASH'
set -euo pipefail
cd "$HOME/KinBot"
git pull --ff-only origin composite
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"
"$PWD/.venv/bin/python" "$PWD/examples/anl/tools/requested_orbital_recon.py" review \
  "$HOME/kinbot-runs/r7_ch2sh_cs_orbital_recon" \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/ch3s/input.json" \
  "$PWD/examples/anl/ch3s_anl1_f12/orbital_reviews/ch2sh_cs_selected.json" \
  --apply
BASH
```

Omit `--apply` for a read-only selection preview. This operation adds only
the reviewed CH2SH channel under `vrc_tst_surface_models`. It preserves the
current scan/sampling grids, bases, extrapolation, flux target, scratch and
resource settings. No production driver is launched or restarted and no
reconnaissance calculation is submitted. Existing QC allocations remain
owned by the running drivers. The new VRC work starts when the owner and
worker consume the revision at their normal safe polling boundaries.

```bash
"$HOME/KinBot/.venv/bin/python" -m kinbot.run_control resubmit \
  "$HOME/kinbot-runs/ch3s_ch3ss_caspt2_production_20261021_r7/ch3s/input.json" --status
```

Expect `queued`, then `published`/`applied`, with `drivers_started: 0`.
`applied` records acceptance of the input revision, not completion of VRC.
If submission is blocked, the edited input remains a recoverable proposal
and the previous input is saved under `ch3s/orbital_reviews/live_refresh`;
the live owner's accepted snapshot is unchanged until publication.

The [complete channel handoff](VRC_CHANNEL_HANDOFF_R7.md) explains how this
selection relates to the three earlier reviewed channels. The historical
screen/candidate commands remain available for deliberate new studies but
are not steps in this production continuation.
