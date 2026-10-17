import pytest
from kinbot.reference_identity import reference_smiles_match


@pytest.mark.parametrize('left,right,match', [
    ('[S]S', 'S=[SH]', True),  # native HSS reference from Blodgett
    ('C[S]', '[CH3][S]', True),
    ('[S]S', '[S][S]', False),  # losing a hydrogen
    ('CCO', 'COC', False),     # isomer, not an alternative representation
    ('CC=O', 'C=CO', False),   # proton transfer
    ('[O][O]', 'O=O', False),  # different unpaired-electron count
    ('[NH4+]', 'N', False),
    ('[2H]S[S]', '[H]S[S]', False),
    ('F/C=C/F', 'F/C=C\\F', False),
])
def test_reference_identity_checks(left, right, match):
    assert reference_smiles_match(left, right) is match


def test_export_preserves_requested_reference_with_perceived_provenance(tmp_path, monkeypatch):
    import json
    import numpy as np
    from types import SimpleNamespace
    from kinbot.mess import MESS
    from kinbot.anl.automatic import _model_record
    species = SimpleNamespace(smiles='S=[SH]', charge=0, mult=2, natom=3,
        atom=['S', 'S', 'H'], geom=np.array([[0., 0., 0.], [1.98, 0., 0.], [2.3, 1.3, 0.]]),
        energy=-795., zpe=.01, sigma_ext=1, nopt=1, dihed=[])
    writer = object.__new__(MESS)
    writer.par = {'smiles': '[S]S', 'thermochemistry_model_output': str(tmp_path / 'model.json')}
    monkeypatch.setattr(writer, '_stable_frequencies', lambda s: [300., 700., 2500.])
    monkeypatch.setattr(writer, 'make_rotors', lambda *a: '')
    monkeypatch.setattr('kinbot.mess.routing_name', lambda s: 'reference')
    model_path = writer.write_thermochemistry_model(species, tmp_path / 'model.json')
    model = json.loads(model_path.read_text())['species']['[S]S']
    assert model['source']['perceived_smiles'] == 'S=[SH]'
    assert model['source']['requested_smiles'] == '[S]S'
    assert _model_record(model_path, {'smiles': '[S]S', 'charge': 0, 'multiplicity': 2}, {})['smiles'] == '[S]S'
    writer.par['smiles'] = '[S][S]'
    with pytest.raises(ValueError, match='connectivity differs'):
        writer.thermochemistry_model(species)
