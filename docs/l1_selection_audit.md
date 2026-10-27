# L1 energy selection and read-only inventory

`python -B -m kinbot.l1_audit CASE --limit 150 --vdw-depth 0.5`
prints JSON without changing the case, its SQLite databases, decisions, or
scheduler state. Redirect stdout to a new report outside the run. Use
`python -B -m kinbot.l1_audit REPORT.json --summarize` for the compact inventory.
Optional repeated `--review-plan SCREENS.json` arguments match requested orbital
reconnaissance recipes to the reported channels. This does not launch a panel.

The audit compares **L1 electronic energy plus harmonic ZPE** with the initial
PES well at the same recorded L1 theory. The 150 kcal/mol ceiling is absolute
relative to that initial well, not relative to each newly explored well.
Missing ZPE, mismatched stoichiometry/theory/charge, and ambiguous IRC endpoints
remain `unknown`; they cannot justify deleting work. Special legacy MP2/BLS
barrier calculations need a matched reference. A family name alone does not
prove a different theory: profiled inputs can resolve scan/BLS method and basis
to the same L1 model. Explicit different overrides remain unclassified.

For neutral unimolecular searches, the two optimized L1 IRC endpoints are
assigned with the same stereochemical identity test as the owning IRC driver.
The product-side fragment identities must agree with the accepted product map.
Endpoint assignment does not require endpoint ZPE. Sella/FC IRC templates save
optimized endpoint geometry and electronic energy without a Hessian; ordinary
single-product reactions can still use the independently completed product-well
E+ZPE after the endpoint identity is established.
The binding depth is the sum of isolated fragment E+ZPE minus complex E+ZPE.
A complex is retained only when its depth is **strictly greater** than the
requested threshold. `vdW_detection` already defaults to 0.5 kcal/mol; the audit
also reports each worker's saved value and live search override. Reapplying the
same threshold may therefore remove no additional complexes.

Schema 2 also reports two depth diagnostics, **without applying them to the
selection**:

- `electronic_only`: isolated-fragment electronic energies minus the endpoint
  electronic energy, omitting ZPE on both sides.
- `legacy_zero_zpe_diagnostic`: the calculation made by the existing detector
  when `get_qc_zpe` substitutes zero for an absent endpoint ZPE, while retaining
  fragment ZPE. This can inflate apparent binding and is not a matched 0 K
  binding energy. It is omitted for an invalid/nonfinite ZPE.

`threshold_only_replay` compares the saved `vdW_detection` threshold with the
requested threshold using that legacy formula (or matched E+ZPE when complete).
It answers what changing the threshold alone would do with the same inputs;
it does not certify a physical depth or reconstruct a historical attempt.
The main `vdw`, VRC, cascade, and work-selection decisions remain `unknown`
where matched ZPE is missing. No missing ZPE is silently supplied from L2,
from separated fragments, or by zero. Choosing electronic-only filtering or
calculating the missing L1 endpoint Hessians is a separate model/workflow choice
before a production handoff. This audit does not change the existing detector.

A shallow complex loses its separate complex-to-fragments VRC channel. Its
preceding tight-TS reaction can remain. A homolytic channel is not a vdW complex
and is not rejected by the depth threshold. When depth is unknown, the audit
does not use a high fragment asymptote to erase a potentially lower complex.

The report distinguishes:

- `exclude`: an above-limit well, a directly screened reaction/VRC channel, or
  an existing operator rejection;
- `retain`: the required energy checks pass;
- `unknown`: insufficient comparable evidence; keep for review;
- `disconnected_by_screen_in_saved_graph`: a downstream-well removal candidate,
  after allowing alternate reversible single-product paths and unknown edges.

Saved-graph disconnection is not proof that an unfinished exploration cannot
find another path. Bimolecular products are not treated as independent hops to
individual wells. A component already disconnected in the source is not counted
as a consequence of the new cutoffs.

Current ANL and VRC graph ledgers are counted once, outside attempt archives.
The report separates tasks under high-energy workers, tasks dedicated to
excluded VRC channels, and disconnected-worker candidates. These are recorded
statuses, not validated restartability or CPU savings. Shared species/CBH tasks
may remain necessary elsewhere. Sampling artifacts and native checkpoint
integrity require the separate seed/recovery audit.

Every potentially retained VRC channel carries its existing review request and
matching model declarations. `model_accepted` is a historical receipt, not a new
approval for a changed geometry, electronic state, or seed contract. A missing
request likewise does not establish approval. Unmatched panel recipes remain
unknown and must not be silently dropped.

## Optional search-time ceiling

`"l1_pes_energy_ceiling": 150` enables the additional L1 screen for new reaction
searches. Its default is `null`, preserving the existing workflow. The driver
screens completed L1 saddles before IRC and completed L1 products/retained
complexes before their L2 work. Existing barrier rules continue to apply.
It records decisions in `l1_screening.json` and preserves native evidence for
excluded searches. A new driver retains those terminal decisions. Accepted
reaction replay is not silently re-screened; applying a different ceiling to
previously accepted networks requires an explicit selection/migration.

This feature does not itself migrate a run, trim its files, validate checkpoints,
or restart controllers. The audit is the evidence input to that later handoff.

Validation includes SQLite read-only byte comparisons, real endpoint identity
assignment, missing/mismatched evidence, alternate/reverse paths, depth and
energy decisions, panel mapping, and rejection persistence through the owning
reaction driver. These are local simulated calculations, not native remote
restarts.
