"""Render trusted, saved PESViewer inputs without running a KinBot driver.

Run this module with the Python interpreter from the PESViewer environment.
PESViewer's GUI entry point and molecular depiction generation are bypassed;
its input reader, layout and energy-diagram renderer are used directly.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def render() -> dict:
    """Render one copied input in a fresh process and output directory."""
    import math
    from importlib.metadata import distribution, version
    import pesviewer.pesviewer as viewer

    # Upstream explicitly selects TkAgg on import, overriding MPLBACKEND.
    # Switch before creating any figures; this also works without tkinter.
    viewer.matplotlib.use('Agg', force=True)
    viewer.options['save_from_command_line'] = 1
    viewer.read_input('pesviewer.inp')
    identifier = viewer.options['id']
    if not identifier or identifier in ('.', '..') or Path(identifier).name != identifier:
        raise ValueError('PESViewer id must be a plain filename component')
    viewer.options.update(save=1, show_images=0, graph_plot=0, use_xyz=0)
    for state in viewer.wells + viewer.bimolecs + viewer.tss + viewer.barrierlesss:
        viewer.linesd[state] = []
    viewer.position()
    viewer.generate_lines()

    def bounds():
        # Upstream get_sizes() fails with wells but no TS/products and can
        # clip submerged TSs. These are plotting bounds, not energy changes.
        states = viewer.wells + viewer.bimolecs + viewer.tss
        if not states or any(not math.isfinite(s.y) for s in states):
            raise ValueError('Missing or nonfinite plotted energies')
        viewer.xlow = min(s.x for s in states)
        viewer.xhigh = max(s.x for s in states)
        viewer.ylow = min(s.y for s in states)
        viewer.yhigh = max(s.y for s in states)
        margin = max(float(viewer.options['margin']), 0.05)
        viewer.xmargin = margin * max(viewer.xhigh - viewer.xlow, 1.)
        viewer.ymargin = margin * max(viewer.yhigh - viewer.ylow, 1.)

    viewer.get_sizes = bounds
    viewer.plot()
    units = viewer.options['display_units']
    rows, labels = [], {}
    for prefix, states in [('W', viewer.wells), ('P', viewer.bimolecs), ('TS', viewer.tss)]:
        for i, state in enumerate(states, 1):
            label = f'{prefix}{i}'
            labels[state] = label
            text = f'{label}: {state.y:.2f}'
            if state in viewer.textd:
                viewer.textd[state].set_text(text)
            else:
                viewer.plt.gca().annotate(text, (state.x, state.y),
                                        xytext=(0, -12), textcoords='offset points',
                                        ha='center', fontsize=8)
            rows.append([label, state.name, state.y, units])
    with Path('species_key.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['label', 'name', 'plotted_energy', 'units'])
        writer.writerows(rows)
    with Path('channels.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['name', 'reactant', 'product', 'transition_state'])
        for channel in viewer.tss + viewer.barrierlesss:
            writer.writerow([channel.name, labels[channel.reactant],
                             labels[channel.product], labels.get(channel, '')])
    viewer.plt.savefig('pes.png', bbox_inches='tight')
    viewer.plt.savefig('pes.pdf', bbox_inches='tight')
    viewer.plt.close('all')
    result = {
        'wells': len(viewer.wells), 'products': len(viewer.bimolecs),
        'transition_states': len(viewer.tss), 'barrierless': len(viewer.barrierlesss),
        'backend': viewer.matplotlib.get_backend(), 'units': units,
        'python': sys.executable,
        'packages': {name: version(name) for name in
                     ('pesviewer', 'matplotlib', 'numpy', 'rdkit', 'networkx')},
        'viewer_source': distribution('pesviewer').read_text('direct_url.json'),
        'display_overrides': {'save': 1, 'show_images': 0, 'graph_plot': 0, 'use_xyz': 0},
    }
    Path('render.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def export(source: Path, output_parent: Path, *, pattern='**/pesviewer.inp', timeout=180) -> dict:
    source = source.resolve(strict=True)
    base = source.parent if source.is_file() else source
    output_parent = output_parent.resolve()
    if output_parent.is_relative_to(base):
        raise ValueError('Output must be outside the source run directory')
    inputs = [source] if source.is_file() else sorted(source.glob(pattern))
    inputs = [path for path in inputs if path.is_file()]
    if not inputs:
        raise ValueError(f'No saved inputs matched {pattern!r} in {source}')
    output_parent.mkdir(parents=True, exist_ok=True)
    destination = Path(tempfile.mkdtemp(prefix='pesviewer-', dir=output_parent))
    print('OUTPUT_DIRECTORY', destination, flush=True)
    report = {'source': str(source), 'directory': str(destination), 'renders': []}
    env = dict(os.environ, MPLBACKEND='Agg', MPLCONFIGDIR=str(destination / '.matplotlib'))
    for path in inputs:
        relative = path.relative_to(base)
        folder = destination / relative.parent
        folder.mkdir(parents=True, exist_ok=True)
        content = path.read_bytes()
        copied = folder / 'pesviewer.inp'
        copied.write_bytes(content)
        log = folder / 'viewer.log'
        item = {'source': str(path), 'input_sha256': hashlib.sha256(content).hexdigest(),
                'directory': str(folder), 'log': str(log)}
        try:
            with log.open('w') as handle:
                result = subprocess.run([sys.executable, '-B', str(Path(__file__).resolve()),
                                         '--render'], cwd=folder, env=env, stdout=handle,
                                        stderr=subprocess.STDOUT, timeout=timeout)
            item['returncode'] = result.returncode
            item['status'] = 'complete' if result.returncode == 0 else 'failed'
            if result.returncode == 0:
                for name in ('pes.png', 'pes.pdf', 'species_key.csv', 'channels.csv', 'render.json'):
                    data = (folder / name).read_bytes()
                    if not data:
                        raise ValueError(f'Empty output: {name}')
                item['render'] = json.loads((folder / 'render.json').read_text())
        except (OSError, ValueError, subprocess.TimeoutExpired) as error:
            item.update(status='failed', error=str(error))
        report['renders'].append(item)
        (destination / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
        print(item['status'].upper(), relative.parent, 'LOG', log, flush=True)
        if item['status'] == 'failed':
            print(item.get('error', '\n'.join(log.read_text(errors='replace').splitlines()[-15:])),
                  flush=True)
            # An import/dependency failure affects every input: report it once.
            if 'ModuleNotFoundError' in log.read_text(errors='replace'):
                break
    report['successful'] = sum(item['status'] == 'complete' for item in report['renders'])
    report['failed'] = sum(item['status'] != 'complete' for item in report['renders'])
    report['not_attempted'] = len(inputs) - len(report['renders'])
    (destination / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print('FINISHED', json.dumps({k: report[k] for k in ('successful', 'failed', 'not_attempted')}),
          flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', type=Path, nargs='?')
    parser.add_argument('--output-parent', type=Path)
    parser.add_argument('--pattern', default='**/pesviewer.inp')
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--render', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.render:
        print(json.dumps(render()))
        return
    if args.source is None or args.output_parent is None:
        parser.error('source and --output-parent are required')
    result = export(args.source, args.output_parent, pattern=args.pattern, timeout=args.timeout)
    if result['failed'] or result['not_attempted']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
