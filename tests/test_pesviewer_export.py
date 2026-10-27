"""Optional native viewer tests: run with the separate PESViewer interpreter."""
import csv
import importlib.util
from pathlib import Path

import pytest

from kinbot.pesviewer_export import export


def input_text(wells='A 0', products='', ts='', barrierless=''):
    return f'''\
> <id> example
> <options>
units kcal/mol
save 0
show_images 1
graph_plot 1
> <wells>
{wells}
> <bimolec>
{products}
> <ts>
{ts}
> <barrierless>
{barrierless}
'''


def test_reject_output_inside_source(tmp_path):
    source = tmp_path / 'source'
    source.mkdir()
    (source / 'pesviewer.inp').write_text(input_text())
    with pytest.raises(ValueError, match='outside'):
        export(source, source / 'plots')
    assert sorted(p.name for p in source.iterdir()) == ['pesviewer.inp']


@pytest.mark.skipif(importlib.util.find_spec('pesviewer') is None,
                    reason='Requires the separate PESViewer environment')
@pytest.mark.parametrize('content, expected_energies, channel_count', [
    (input_text(), [0.], 0),
    (input_text(products='B 3', barrierless='dissociation A B'), [0., 3.], 1),
    (input_text(wells='A 0\nB -5', ts='T -2 A B'), [0., -5., -2.], 1),
    (input_text(wells='A 0\nB 0', ts='T1 2 A B\nT2 2 A B'), [0., 0., 2., 2.], 2),
])
def test_native_headless_export_preserves_source_and_energies(
        tmp_path, monkeypatch, content, expected_energies, channel_count):
    monkeypatch.delenv('DISPLAY', raising=False)
    source = tmp_path / 'source' / 'case' / 'well'
    source.mkdir(parents=True)
    (source / 'pesviewer.inp').write_text(content)
    report = export(tmp_path / 'source', tmp_path / 'plots', pattern='*/*/pesviewer.inp')
    assert report['successful'] == 1, report
    assert report['failed'] == report['not_attempted'] == 0
    folder = Path(report['renders'][0]['directory'])
    assert (folder / 'pes.png').read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    assert (folder / 'pes.pdf').read_bytes().startswith(b'%PDF')
    with (folder / 'species_key.csv').open() as handle:
        rows = list(csv.DictReader(handle))
    assert [float(row['plotted_energy']) for row in rows] == expected_energies
    assert {row['units'] for row in rows} == {'kcal/mol'}
    with (folder / 'channels.csv').open() as handle:
        assert len(list(csv.DictReader(handle))) == channel_count
    assert report['renders'][0]['render']['backend'].lower() == 'agg'
    assert sorted(p.name for p in source.iterdir()) == ['pesviewer.inp']
    assert (source / 'pesviewer.inp').read_text() == content
    assert (folder / 'pesviewer.inp').read_text() == content
