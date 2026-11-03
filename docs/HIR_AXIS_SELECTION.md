# Hindered-rotor axis selection

The ordinary rule is unchanged: a non-ring axis must be single-bonded in every
generated resonance structure. Mixed-order C-C, C-O, C-N and other non-sulfur
axes remain excluded. A double or triple bond in every form remains excluded,
including C=S and S=S. Merely being adjacent to sulfur is not an exception.

The new exception applies only when **an axis atom is sulfur**, the bond is
present in every generated form, and at least one form assigns it order one.
KinBot's expanded-valence sulfur assignments can otherwise veto a candidate
solely because another form assigns a higher bond order. This admits a scan;
it does not assume free rotation, a small barrier, or mirror coverage. The
existing native-point validation, potential fit, selected-Hessian projection
and optical-domain evidence checks still decide whether HIR is usable.
Other expanded-valence elements are not automatically enabled by this change.

The central bond must also be a graph bridge: removing it must disconnect its
two tops. This fixes the old ring-atom test, which incorrectly excluded a bond
joining two rings. Actual ring edges remain excluded. In a transition state,
the graph includes forming and breaking bonds, so a reaction-center loop is
not treated as an independent one-dimensional rotor. Nonterminal, noncollinear
outer reference atoms are still required. Existing eligible axes and reference
choices are retained in their prior order; newly eligible axes are appended.
Internal rotor symmetry uses the same bond-order eligibility rule.

## Native rxn20 diagnosis

The captured R8 calculation had `rotor_scan=1` and `nrotation=12`; its exported
TS was RRHO-only. Reconstructing the TS exactly as the reaction generator does,
from L1 saddle row 498, reactant row 2 and atom-mapped IRC product row 574,
explains the omission. Its generated TS resonance matrices assign orders
`[3,3]` to C1-S2 and `[3,1]` to S2-S3. The previous all-single requirement
therefore selected neither axis. Those integers are graph assignments, not
computed electronic bond orders.

With the new rule, **S2-S3** is selected, with dihedral H6-S2-S3-H4 (one-based
indices). C1-S2 remains excluded: the TS union has the C1-S2-H6 reaction loop.
The selected L2 saddle (row 693) generates twelve 30-degree scan inputs without
changing distances within either top. The fixture records the jobs, database
rows, coordinates and database hash in
`tests/reference/rotor_reaction_endpoints.json`.

This is evidence of correct setup, **not a completed native S-S potential**.
Neither the mirror factor nor the existing two MESS connections is removed
merely because the axis is now detected. Only a usable measured scan can
establish that its represented domain covers a mirror. The native mode-parity
and duplicate-surface checks remain separate.

## Updating a retained run

An ordinary controller refresh retains the accepted-search replay guard; it
does not authorize new missing characterization. Use the existing network
reassessment with `--vdw-detection` set to the run's **unchanged** current value.
It reconstructs accepted reactions from their saved TS/IRC calculations and
allows their missing conformer/HIR characterization. It does not restart
terminal rejected searches unless a separate rejected-search option is supplied.
Live calculations remain owned and are adopted by the refreshed controllers.

Previously computed scans are checked against their exact defining atoms and
geometry. A changed scan gets a separate recovery name rather than overwriting
old evidence. Successful and running compatible scans are reused; a failed
replacement is not resubmitted indefinitely. New conformational sampling can
discover a lower geometry, which follows the existing optimization/refinement
workflow rather than overwriting an accepted geometry without calculation.

Regression coverage includes native TS/IRC reconstruction under atom
permutations, the real Gaussian/Sella HIR input writer and repeat adoption,
non-sulfur exclusions, sulfur periodicity, graph bridges, unchanged old indices,
and complete controller reassessment with live QC/sampling preservation and
closed failed searches. These are local tests, not a verified Blodgett retry.
