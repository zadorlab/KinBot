# PR #100: saved THF conformers in final MESS inputs

The direct writer on master gives both members the parent's ground energy.
The second member's actual offset is written only in a comment, which MESS
does not use. The PR moves that offset into `ZeroEnergy`. The existing PES
route already applies the offset for this well; it retains the same values.

| Final input | Master `6a456c2` | PR #100 with review changes |
|---|---|---|
| Direct, member 1 | 0.00 kcal/mol | 0.00 kcal/mol |
| Direct, member 2 | 0.00 kcal/mol | 0.13 kcal/mol |
| PES, member 1 | 0.00 kcal/mol | 0.00 kcal/mol |
| PES, member 2 | 0.13 kcal/mol | 0.13 kcal/mol |

The unrounded second-member offset is **0.1348964193180358 kcal/mol**,
calculated from each member's electronic energy plus its own ZPE.
Both final routes retain both geometries and all 33 real frequencies per
member. These are saved isolated tetrahydrofuran calculations from
`THF_OH_new`, using FairChem `uma-m-1p1` / `omol`. The fixture preserves
source row IDs 32 and 26, source job names, coordinates, raw frequencies,
energies, ZPEs and units in `tests/reference/thf_mc_rrho.json`.

Relevant direct `mess.inp` change:

```diff
-          ZeroEnergy[kcal/mol]          0.0 ! 0.1348964193180358
+          ZeroEnergy[kcal/mol]          0.13 !
```

In the PES file the numerical `ZeroEnergy` values remain 0.00 and 0.13;
the differences are removal of offset comments and whitespace preservation.
The automated test exercises `MESS.write_input` and `pes.create_mess_input`,
including the final files, rather than testing only intermediate strings.

Run from a checkout with its dependencies installed:

```sh
PYTHONPATH=. python tests/test_real_mess_conformers.py
PYTHONPATH=. python tests/test_real_mess_conformers.py --write-example /tmp/thf-mess-review
```

The second command saves `direct_mess.inp` and `pes_mess.inp`. To compare
another revision, run this same script with `PYTHONPATH` pointing to that
revision's checkout and a different output directory.

This is a replay of real calculated molecular properties through the writers.
It is not a new conformer search, a reaction-network validation or a MESS
rate calculation. TS-specific tunneling offsets remain covered by the
separate final-output regressions in `tests/test_mess_conformers.py`.
