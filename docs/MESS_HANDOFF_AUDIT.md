# Complete MESS assembly evidence

When exploration is finished but final assembly fails, a monitor tail is not a
complete network regression. Capture both PES cases together for local replay:

```bash
"$HOME/KinBot/.venv/bin/python" -B -m kinbot.run_seed_audit "$KB_RUN" \
  --mess-handoff --output-parent "$HOME/kinbot-evidence"
```

This reads the source and writes a new evidence directory and `seed-evidence.tar.gz`
outside it. It does not execute a saved worker, submit or cancel jobs, refresh
controllers, or change acceptance decisions. The capture includes:

- Complete worker MESS fragments, summaries and the ordered `chemids` lists.
- Full databases and SQLite sidecars, preserving selected observations and
  conformer records, plus shallow input, search and decision JSON files.
- PES/worker controller logs, `me` files and XYZ structures.

It does not traverse QC scratch, orbital campaigns or sampling trees. It has a
2 GiB total default limit and a 256 MiB per-file limit; omissions, symlinks,
changing files or a changing inventory produce an explicit incomplete result
and a nonzero CLI exit. Every copied source has a content hash, and all sources
are rechecked at the end. This is not an atomic snapshot of a live run: capture
the completed, inactive handoff rather than moving calculations underneath it.

The capture status means that the selected files were copied consistently, not
that all network models are scientifically accepted or all required artifacts
exist. Keep missing artifacts as diagnostic evidence. Replay on a disposable
copy of the captured run, using its original inputs and all workers, with QC
submission disabled. Preserve source databases, accepted decisions and native
files; repair final network serialization only when the route evidence supports
it. Check final populations, barrier connections, route counting, energy
references, conformers, rotors and PST exits across both complete networks.
Then validate the generated inputs with the native MESS executable and retain
its input hashes, output and exit status. Successful Python assembly alone is
not solver acceptance, and solver acceptance alone is not a branching analysis.

`tests/test_run_seed_audit.py` exercises native-data worker artifacts through
capture, archive extraction and two repeated real PES assemblies, checking
source immutability and all archived hashes. Other tests cover SQLite WAL,
capacity limits, symlink refusal and newly appearing artifacts. These validate
the evidence handoff; they do not resolve a particular production network's
parallel-route ambiguity.
