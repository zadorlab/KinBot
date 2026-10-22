from copy import deepcopy

import numpy as np
import pytest
from ase import Atoms
from ase.build import molecule
from ase.db import connect

from kinbot import constants
from kinbot.anl.automatic import _model_record
from kinbot.anl.dispatch import _geometry_hash
from kinbot.anl import pes_ts_reuse as reuse, reference_frequencies as recovery
from kinbot.execution import atomic_json


@pytest.fixture
def reference(tmp_path, monkeypatch):
    atoms = molecule('CH3')
    native = [514.4159, 1427.6941, 1427.7039, 3139.5836, 3320.8406, 3320.8435]
    with connect(str(tmp_path/'kinbot.db')) as db:
        row = db.write(atoms, name='methyl_well_high', data=dict(status='normal',
            energy=-39.811242216002/constants.EVtoHARTREE, zpe=.02996, frequencies=native))
    model = dict(name='methyl', smiles='[CH3]', charge=0, multiplicity=2,
        symbols=list(atoms.symbols), positions=atoms.positions.tolist(), symmetry_factor=6.,
        hindered_rotor_source='none-no-rotatable-bonds', hindered_rotor_blocks='',
        frequencies_cm_inverse=[180.0851, 514.4160, 793.5796, 2058.8831, 2745.2693, 3265.5741],
        source=dict(source_job='methyl_well_high', source_row_id=row,
            electronic_energy_hartree=-39.811242216002, harmonic_zpe_hartree=.02996,
            geometry_sha256=_geometry_hash(atoms)))
    path = tmp_path/'thermochemistry_model.json'
    atomic_json(path, dict(schema=1, species={'[CH3]': model}))
    entry = dict(smiles='[CH3]', charge=0, multiplicity=2)
    with monkeypatch.context() as m:
        m.setattr(recovery, 'recover', lambda r:r)
        old = _model_record(path, entry, {'calculator':'gaussian', 'basis':'cc-pVTZ'})
    return old, native, path, entry


def test_saved_native_row_repairs_export_without_mutating_source(reference):
    old, native, path, entry = reference
    contents = path.read_bytes()
    new = _model_record(path, entry, old['l2_theory'])
    assert new['frequencies_cm_inverse'] == native
    assert new['thermochemistry_model']['frequencies_cm_inverse'] == native
    assert new['molecule'] == old['molecule']  # Existing QC workflow provenance stays valid.
    assert path.read_bytes() == contents
    assert recovery.verified_migration(old, new)
    assert new['frequency_recovery']['source_row_id'] == 1


def test_registry_accepts_only_rederived_frequency_migration(reference, tmp_path):
    old, native, path, entry = reference
    new = recovery.recover(old)
    item = dict(record=old, graph=str(tmp_path/'graph'))
    registry = dict(records={'ref':reuse.signature(old)}, catalog={'ref':deepcopy(item)})
    current = {'ref':dict(item,record=new)}
    actual = reuse.history(tmp_path, registry, current)
    assert actual['ref']['record']['frequencies_cm_inverse'] == native
    assert registry['records']['ref'] == reuse.signature(new)
    assert registry['frequency_migrations']['ref']['previous_record'] == old
    for field in ('energy_hartree','harmonic_zpe_hartree'):
        bad = deepcopy(new); bad[field] += .001
        with pytest.raises(ValueError, match='provenance'):
            reuse.history(tmp_path, dict(records={'ref':reuse.signature(old)},
                catalog={'ref':item}), {'ref':dict(item,record=bad)})


def test_wrong_selected_row_is_not_replaced_by_latest_matching_job(reference):
    old, native, path, entry = reference
    with connect(str(path.parent/'kinbot.db')) as db:
        row = db.get(id=1)
        db.write(row.toatoms(), name=row.name, data=dict(row.data))
        db.update(1, data=dict(row.data, zpe=.5))
    with pytest.raises(ValueError, match='selected L2 row'):
        recovery.recover(old)


def test_rotor_reduced_spectrum_not_replaced_by_full_harmonic_list(reference):
    old, native, path, entry = reference
    old['thermochemistry_model'].update(hindered_rotor_source='kinbot_l2',
        hindered_rotor_blocks='Rotor Hindered\nEnd\n', frequencies_cm_inverse=[999.]*5)
    new = recovery.recover(old)
    assert new['frequencies_cm_inverse'] == native
    assert new['thermochemistry_model']['frequencies_cm_inverse'] == [999.]*5
    assert not new['frequency_recovery']['thermal_spectrum_recovered']
    assert new['frequency_recovery']['thermal_spectrum_requires_reexport']


def test_recovered_harmonic_spectrum_enables_strict_owner_reconciliation(reference, tmp_path, monkeypatch):
    from kinbot.anl import ts_reuse
    old, native, path, entry = reference
    canonical = deepcopy(old)
    canonical['frequencies_cm_inverse'] = native
    def item(record, name):
        return dict(record=record, kind='stable', entry=entry,
            graph=str(tmp_path/name), base=name, name='stable:[CH3]|0|2',
            plan=dict(resources={},l2_theory=old['l2_theory'],requested_method='ANL1-F12'))
    a, b = item(canonical,'a'), item(old,'b')
    limits=ts_reuse.options({})
    assert reuse.compare(a,b,limits) is None
    registry=dict(schema=1,owners=['a','b'],aliases={},evidence={},catalog={'a':a,'b':b},
        records={'a':reuse.signature(canonical),'b':reuse.signature(old)})
    atomic_json(tmp_path/reuse.NAME,registry)
    b=dict(b,record=recovery.recover(old))
    monkeypatch.setattr(reuse,'catalog',lambda *a:{'a':item(canonical,'a'),'b':b})
    mapping,_=reuse.discover(tmp_path,limits,reconcile_owners=True)
    assert mapping['aliases']=={'b':'a'}
    assert mapping['frequency_migrations']['b']['previous_record']==old
    assert mapping['evidence']['b']['measurements']['frequency_cm_inverse']==0.


@pytest.mark.parametrize('name,native,valid', [
    ('H', [], True), ('He', [], True), ('H2', [4400.], True),
    ('CH3', [514.4159, 1427.6941, 1427.7039, 3139.5836, 3320.8406, 3320.8435], True),
    ('H', [0.], False), ('H2', [], False),
])
def test_qc_writer_numpy_spectra_survive_database_roundtrip(tmp_path, name, native, valid):
    atoms = Atoms(name) if name in ('H', 'He') else molecule(name)
    smiles = {'H':'[H]', 'He':'[He]', 'H2':'[H][H]', 'CH3':'[CH3]'}[name]
    mult = 2 if name in ('H', 'CH3') else 1
    energy, zpe = -1., 0. if len(atoms) == 1 else .01
    job = name + '_well_high'
    with connect(str(tmp_path/'kinbot.db')) as db:
        row_id = db.write(atoms, name=job, data=dict(status='normal',
            energy=energy/constants.EVtoHARTREE, zpe=zpe,
            frequencies=np.array(native)))
        assert isinstance(db.get(id=row_id).data['frequencies'], np.ndarray)
    model = dict(name=name, smiles=smiles, charge=0, multiplicity=mult,
        symbols=list(atoms.symbols), positions=atoms.positions.tolist(), symmetry_factor=1.,
        frequencies_cm_inverse=[1.] * len(native), hindered_rotor_blocks='',
        hindered_rotor_source='none-no-rotatable-bonds',
        source=dict(source_job=job, source_row_id=row_id,
            electronic_energy_hartree=energy, harmonic_zpe_hartree=zpe,
            geometry_sha256=_geometry_hash(atoms)))
    path = tmp_path/'thermochemistry_model.json'
    atomic_json(path, dict(schema=1, species={smiles:model}))
    original = path.read_bytes()
    if not valid:
        with pytest.raises(ValueError, match='invalid minimum frequencies'):
            _model_record(path, dict(smiles=smiles, charge=0, multiplicity=mult),
                          {'calculator':'gaussian'})
        return
    record = _model_record(path, dict(smiles=smiles, charge=0, multiplicity=mult),
                           {'calculator':'gaussian'})
    assert record['frequencies_cm_inverse'] == native
    assert record['thermochemistry_model']['frequencies_cm_inverse'] == native
    assert path.read_bytes() == original
    with connect(str(tmp_path/'kinbot.db')) as db:
        np.testing.assert_array_equal(db.get(id=row_id).data['frequencies'], native)


@pytest.mark.parametrize('native', [None, np.array([]), np.array([514.]),
    np.array([[1., 2., 3., 4., 5., 6.]]), np.array([1., 2., 3., 4., 5., np.nan]),
    np.array([1., 2., 3., 4., 5., np.inf]), np.array([-100., 2., 3., 4., 5., 6.]),
    ['1', '2', '3', '4', '5', '6'], np.array([True]*6)])
def test_malformed_native_spectrum_still_blocks_recovery(reference, native):
    old, _, path, _ = reference
    with connect(str(path.parent/'kinbot.db')) as db:
        db.update(1, data={'frequencies':native})
    with pytest.raises(ValueError, match='selected L2'):
        recovery.recover(old)
