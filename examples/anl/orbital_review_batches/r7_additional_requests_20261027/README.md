# Four additional requests captured on 2026-10-26 at 20:38 PDT

The source inventory contains eight accepted review records and 23 requests
still needing review. Four of the latter, from the `s8fc26f5f496eeb4e` worker,
did not have recipes in the prior panels. These 16 independent SCF jobs inspect
their actual atom assignments at near and 30 Å separations, in both DZ and TZ.
They select no CAS and approve no production model.

* C–S cleavage: CH2 (triplet, atoms 1/4/5) + HSS (doublet, 2/3/6).
* S–S cleavage: CH2SH (doublet, 1/2/4/5/6) + S (triplet, 3).
* Hydrogen exit: HCSSH (singlet, 1/2/3/4/6) + H (doublet, 5), from the
  forward L2 product record (captured row 426). The L1 barrier exceeds 150
  kcal/mol; the requested orbital review still proceeds independently.
* CH2 + HSS exit after insertion: the reverse L2 product record (row 762),
  with CH2 atoms 1/4/5 and HSS atoms 2/3/6. This is not the forward endpoint.

Quartet ROHF character is inspected for the triplet-plus-doublet fragment
pairs; the H exit uses doublet ROHF. This initialization does not choose the
production spin surface, state average, or statistical weight. Follow-up
compact CAS comparisons must include every dissociating bond/radical direction,
check near/far state identity and inactive-space leakage, and justify any added
donor or degenerate orbital using native evidence. Similar connectivity to a
previous stereoisomer does not justify copying its orbital indices.

`run_blodgett.sh` stages beneath the existing detached review pool. Both active
panels share its total three-job cap. No source run files are changed. Local
tests cover all 16 generated decks, forward/reverse endpoint mapping, fragment
translation, preflight and repeated staging; native execution remains separate.
