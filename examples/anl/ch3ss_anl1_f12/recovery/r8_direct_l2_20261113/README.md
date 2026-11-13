# Direct L2 checks of the paper and scouting disagreements

The 25 cases use B2PLYP-D3(BJ)/cc-pVTZ throughout, including constrained
guess preparation, saddle optimization, native Gaussian frequencies, IRCs
and IRC endpoint optimization. The ordinary KinBot QC owners handle jobs,
SCF recovery, checkpointing and the existing shared resource pool.

## Scope

| Cases | Work and purpose |
|---|---|
| Six paper seeds | Optimize the exact SI geometries of TS1–TS6 at L2 and follow both IRC directions |
| Twelve saved KinBot L2 saddles | Reuse their hash-verified native geometry, energy and frequency records, then calculate both L2 IRCs |
| One saved reverse TS1 geometry | Reoptimize and obtain native frequencies before L2 IRCs, because its full frequency output is absent from the supplied handoff |
| Six H-loss guesses | Repeat the three parent-angle, bond-length scales 1.3, 1.6 and 1.9 directly at L2 for CH2SH and CH3SS |

The native comparisons cover both CH3S isomerization saddles and its retained
H-loss saddle, both distinct CH3SS-to-CH2SSH saddles, all three distinct
CH2SSH-to-HSCH2S saddles associated with the TS6 disagreement, three alternative
isomerizations found while targeting TS5, and the existing TS2 connection.
Independent mirror rediscoveries are represented once in this panel.

The SI source is `Theoretical calculations/CH2S+SH/MESS_input.txt` from
DOI 10.1039/d4cp04873c. `paper_geometries.json` retains the coordinate and
endpoint labels with the SI file checksum. The six paper structures are
starting geometries for L2 calculations. Their energies and connectivity
remain observations to test. Every source native input/output is hash-bound
in `campaign.json`, from the November 13 export.

## Scientific completion

Each successful saddle receives both L2 IRCs. The result records the actual
endpoint identities and geometries even when they differ from the intended
connection. Native frequencies and displacement vectors are retained.
`irc_endpoints_available` means that both IRCs converged and both endpoint
optimizations published normal results. It still requires reviewing the
recorded connection, mode and paper correspondence before PES integration.
Partial IRC paths receive `irc_incomplete` and retain their endpoint evidence.
Terminal unsuccessful saddle searches remain terminal on restart.

The Sella saddle and endpoint tolerance is the scouting value of
0.0005 eV/angstrom. Each IRC has a 300-step bound and Sella `dx=0.05`.
The target-level checks apply no barrier-height pruning. The common parameter
reader requires a barrier threshold, which this diagnostic controller does
not use. The calculations run without conformer/HIR expansion, ANL, VRC or
MESS execution. They share the existing R8 node cap and use eight processors
and exclusive allocations. Completed R8 calculation files remain source
evidence. The new jobs reside under `direct_l2_checks_20261113`.

## Run and inspect

From the updated composite checkout on Blodgett run

```bash
bash examples/anl/ch3ss_anl1_f12/recovery/r8_direct_l2_20261113/run_blodgett.sh
```

The controller is detached from the login shell. Its exclusive ownership lock
prevents two controllers from submitting the same case concurrently. Repeating
the launcher adopts completed and active work through the saved QC routes and
reaction-search checkpoints.

Read the current saved status with

```bash
cd "$HOME/KinBot"
export PYTHONPATH="$PWD"
KB_CHECK="$HOME/kinbot-runs/ch3s_ch3ss_l1_l2_scout_r8/direct_l2_checks_20261113"
"$PWD/.venv/bin/python" -B -m kinbot.saddle_recheck status "$KB_CHECK"
cat "$KB_CHECK/controller.json"
```

Export the full text-only result records with

```bash
"$PWD/.venv/bin/python" -B -m kinbot.saddle_recheck status "$KB_CHECK" --details \
  > "$HOME/kinbot-evidence/r8-direct-l2-results.txt"
```

The controller invokes the existing QC and IRC writers. It does not infer
global absence of a pathway from a finite unsuccessful search. Further PES
integration and kinetic calculations follow the native evidence review.
