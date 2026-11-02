"""Mirror completion must apply the same hom_sci/saddle precedence as PES."""
import hashlib
import json
from pathlib import Path
import re

import pytest

from kinbot import pes
from kinbot.mess_mirrors import (annotate_capture_placeholder,
    complete_mirror_channels, _PATH)
from kinbot.mess_networks import _HEADER, split_model
from kinbot.parameters import Parameters
from tests.test_mess_pes_artifacts import at

FIXTURE = Path(__file__).parent / 'reference/mess_mirror_pst_handoff.json'


def stage_mirror_capture(tmp_path, legacy_exit=False):
    from kinbot.species_routing import mess_filename
    row = json.loads(FIXTURE.read_text())
    for name, item in row['files'].items():
        assert hashlib.sha256(item['text'].encode()).hexdigest() == item['sha256']
        target = tmp_path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(item['text'])
    if legacy_exit:
        reaction = row['barrierless'][0]
        product = row['products'][0]
        source = tmp_path/reaction[0]/mess_filename(product, 0, reaction=reaction[1])
        source.rename(source.with_name(mess_filename(product, 0)))
    source = tmp_path / 'input.json'
    source.write_text(json.dumps(row['parameters']))
    par = Parameters(str(source), show_warnings=False).par
    return row, par


def assemble_mirror_capture(tmp_path, row, par, reverse=False):
    wells = row['wells'][::-1] if reverse else row['wells']
    with at(tmp_path):
        pes.create_mess_input(par, wells, row['products'], row['reactions'],
            row['barrierless'], row['complexes'], row['well_energies'],
            row['product_energies'], row['parents'], row['mass'], False)
    networks = json.loads((tmp_path/'me/mess_networks.json').read_text())
    # Sharing a bimolecular exit does not connect these two mirror wells.
    assert len(networks) == 2
    return {n['stem']: (tmp_path/'me'/(n['stem']+'.inp')).read_text() for n in networks}


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('legacy_exit', [False, True])
def test_native_mirror_pst_handoff_assembles_without_duplicate_flux(
        tmp_path, reverse, legacy_exit):
    row, par = stage_mirror_capture(tmp_path, legacy_exit)
    before = {p: p.read_bytes() for p in tmp_path.rglob('*.mess')}
    previous = None
    for _ in range(2):
        models = assemble_mirror_capture(tmp_path, row, par, reverse)
        text = ''.join(models.values())
        blocks = [b for model in models.values() for b in split_model(model)[1]]
        barriers = [b for b in blocks if _HEADER.search(b)[1] == 'Barrier']
        assert len(barriers) == 2
        assert all(_PATH.findall(b) == ['ordinary'] for b in barriers)
        assert all('Core PhaseSpaceTheory' not in b for b in barriers)
        assert 'superseded by calculated saddle' in text
        assert sum('derived by global reflection' in b for b in barriers) == 1
        assert len({tuple(sorted(_HEADER.search(b)[2].split()[1:3]))
                    for b in barriers}) == 2
        assert re.findall(r'ImaginaryFrequency\[1/cm\]\s+(\S+)', text) == ['560.35'] * 2
        assert re.findall(r'ZeroEnergy\[kcal/mol\]\s+(\S+)', ''.join(barriers)) == ['77.44'] * 2
        for model in models.values():
            assert complete_mirror_channels(model).split() == model.split()
        assert all(p.read_bytes() == contents for p, contents in before.items())
        if previous is not None:
            assert models == previous
        previous = models


def test_placeholder_tag_requires_known_homolysis_and_exclusively_pst():
    block = ('  Barrier B W P\n RRHO\n Core PhaseSpaceTheory\n End\n End\n')
    marked = annotate_capture_placeholder(block, 'well_hom_sci_1_2')
    assert '! kinbot_hom_sci_placeholder well_hom_sci_1_2' in marked
    assert annotate_capture_placeholder(marked, 'well_hom_sci_1_2') == marked
    assert annotate_capture_placeholder(block, 'well_other_1_2') == block
    sampled = block.replace('Core PhaseSpaceTheory', 'Core Rotd')
    assert annotate_capture_placeholder(sampled, 'well_hom_sci_1_2') == sampled
    mixed = block.replace(' End\n End', ' End\n Core Rotd\n End\n End')
    assert annotate_capture_placeholder(mixed, 'well_hom_sci_1_2') == mixed
    with pytest.raises(ValueError, match='conflicting'):
        annotate_capture_placeholder(marked, 'well_hom_sci_2_3')


def test_unclassified_parallel_routes_are_still_not_silently_summed():
    from tests.test_mess_mirror_completion import mirror_path_fixture
    text = mirror_path_fixture().replace('! kinbot_stereopath path2\n', '')
    with pytest.raises(ValueError, match='distinct classified'):
        complete_mirror_channels(text)


def test_known_placeholder_does_not_remove_distinct_classified_paths():
    from tests.test_mess_mirror_completion import mirror_path_fixture
    text = mirror_path_fixture()
    placeholder = annotate_capture_placeholder(
        ' Barrier fallback a r\n RRHO\n Core PhaseSpaceTheory\n End\n End\n',
        'well_hom_sci_1_2')
    text = text.replace('End ! end kinetics', placeholder + 'End ! end kinetics')
    final = complete_mirror_channels(text)
    assert final.count('Union ! 2 stereochemical pathways') == 2
    assert final.count('! kinbot_stereopath path1') == 2
    assert final.count('! kinbot_stereopath path2') == 2
    assert 'Core PhaseSpaceTheory' not in final
    assert complete_mirror_channels(final).split() == final.split()
