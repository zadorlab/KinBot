# Plot saved PESViewer inputs on a headless host

Use a separate PESViewer Python environment. The old conda-forge PESViewer
1.1.0 build imports RDKit's removed `cairoCanvas` module and fails before
reading an input. Upstream commit
[`ee0e18fcdb6bc8de6fddaa3e4726a0143c785a73`](https://github.com/zadorlab/PESViewer/commit/ee0e18fcdb6bc8de6fddaa3e4726a0143c785a73)
removes the obsolete drawing APIs, but retains the version number 1.1.0.
Install that revision in the viewer environment, preserving its dependencies:

```bash
"$HOME/.local/envs/pesviewer/bin/python" -m pip install --no-deps --force-reinstall \
  "https://github.com/zadorlab/PESViewer/archive/ee0e18fcdb6bc8de6fddaa3e4726a0143c785a73.tar.gz"
```

From the KinBot checkout, use that same interpreter to export a saved run:

```bash
"$HOME/.local/envs/pesviewer/bin/python" -B -m kinbot.pesviewer_export \
  /path/to/retired-run --pattern '*/*/pesviewer.inp' \
  --output-parent "$HOME/kinbot-evidence"
```

The pattern above selects `case/worker/pesviewer.inp`. A single input file can
also be supplied; the default directory pattern searches recursively. Each
invocation creates a new output directory outside the source. Inputs are copied
byte for byte and their hashes recorded. No driver, scheduler, database or
chemistry calculation is invoked. The source and previous exports stay intact.

The runner uses PESViewer's parser, layout and static renderer, switches from
its hardcoded GUI backend to Agg, and omits molecular depiction and interactive
graph generation. Display bounds include submerged transition states and
support well-only inputs. Each worker gets:

- `pes.png` and `pes.pdf`, with W/P/TS labels and saved energies;
- `species_key.csv`, mapping labels to original names and plotted energies;
- `channels.csv`, mapping saved reactions to those labels;
- `viewer.log`, `render.json` and the copied `pesviewer.inp`.

The export root's `report.json` records failures, input hashes, package versions
and viewer source metadata. The command exits nonzero on failure and prints
the log tail. Inputs need to be trusted: upstream's reader evaluates expressions.

These are separate worker diagrams with the energy references already encoded
in their inputs. They do not combine or align different workers into a global
PES and do not recalculate or complete missing scientific results.

Local validation used native PESViewer at the pinned revision, NumPy 2.4.6,
Matplotlib 3.11.2, NetworkX 3.6.1 and RDKit 2026.3.6. Saved CH3S/CH3SS inputs
rendered successfully. Optional native tests cover a lone well, a barrierless
channel, a submerged TS, equal-height TSs, energy preservation and source
preservation. Run them with the viewer interpreter:

```bash
python -m pytest -q tests/test_pesviewer_export.py
```
