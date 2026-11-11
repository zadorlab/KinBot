from types import SimpleNamespace

import pytest
import numpy as np
from ase import Atoms
from ase.db import connect

from kinbot import constants
from kinbot.execution import file_sha256
from kinbot.pes import get_energy
from kinbot.pes_product_identity import ProductIdentities
from kinbot.run_format import ensure_current_run
from kinbot.stationary_pt import StationaryPoint
from kinbot.stereo_routing import guard_well_job


# Two native L2 geometries whose legacy graph sums collide. The full graphs
# differ in hydrogen attachment and C/S connectivity.
GEOMETRIES = [
    [[1.19216394, .07787878, 1.23704391], [2.38722655, 1.11562936, 1.18716748],
     [-.21189131, -.11531496, 2.16661499], [-.65802893, -1.23300258, 1.56529131]],
    [[.74676084, .22469945, -.03541007], [.78891676, .90113910, 1.38257377],
     [.20687402, -.07400391, 2.97679915], [1.02626123, .61298374, -1.00344290]],
]


@pytest.fixture
def collision(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ensure_current_run(create=True)
    workers = ['first', 'second']
    names = []
    for i, (worker, geometry) in enumerate(zip(workers, GEOMETRIES)):
        ensure_current_run(worker, create=True)
        species = StationaryPoint('fragment', 0, 2, atom=['C', 'S', 'S', 'H'], geom=np.array(geometry))
        species.characterize()
        name = str(species.chemid)
        names.append(name)
        db = connect(worker + '/kinbot.db')
        guard_well_job(SimpleNamespace(db=db), species, species.geom, name + '_well')
        db.write(Atoms(species.atom, positions=species.geom), name=name + '_well_high',
                 data=dict(status='normal', energy=(-100 + i) / constants.EVtoHARTREE,
                           zpe=.01, charge=0, multiplicity=2))
    assert names[0] == names[1]
    return workers, names[0]


def test_catalog_and_energy_selection_preserve_native_files(collision):
    workers, name = collision
    before = {w: file_sha256(w + '/kinbot.db') for w in workers}
    for order in (workers, workers[::-1], workers):
        catalog = ProductIdentities(order)
        assert len(catalog.aliases) == 2
        for i, worker in enumerate(workers):
            alias = catalog.product(worker, name)
            assert alias != name and catalog.original(alias) == name
            assert get_energy(order, alias, 0, 1, product_aliases=catalog.aliases) == pytest.approx((-100+i, .01))
            source = f'! kinbot_population 1_{name} 1_{name}\n{{fr_name_{name}}}\n'
            annotated = catalog.annotate(source, worker)
            assert f'1_{alias}' in annotated
            assert '{fr_name_' + name + '}' in annotated
            assert catalog.original_key(f'1_{alias}') == f'1_{name}'
    assert before == {w: file_sha256(w + '/kinbot.db') for w in workers}


def test_direct_reader_refuses_collision_and_single_identity_is_unchanged(collision):
    workers, name = collision
    with pytest.raises(ValueError, match='different complete stereoisomer identities'):
        get_energy(workers, name, 0, 1)
    catalog = ProductIdentities(workers[:1])
    assert not catalog.aliases
    assert catalog.product(workers[0], name) == name
    assert catalog.original_key(name + '_1') == '1_' + name


def test_alias_cannot_adopt_wrong_geometry_or_unverified_reference(collision):
    workers, name = collision
    catalog = ProductIdentities(workers)
    alias = catalog.product(workers[0], name)
    db = connect(workers[0] + '/kinbot.db')
    db.write(Atoms('CSSH', positions=GEOMETRIES[1]), name=name + '_well_high',
             data=dict(status='normal', energy=-200., zpe=.01))
    with pytest.raises(ValueError, match='Unable to find an energy'):
        get_energy(workers, alias, 0, 1, product_aliases=catalog.aliases)
