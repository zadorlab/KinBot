"""Changed rotor definitions must not relabel previously calculated scans."""
from pathlib import Path
from types import SimpleNamespace
import pytest
from kinbot.qc import QuantumChemistry


@pytest.mark.parametrize('backend', ['sella','gauss','qchem'])
def test_changed_scan_definition_is_rejected_before_overwrite(tmp_path,backend):
    job=str(tmp_path/'scan')
    text={'sella':'base_0_fix = [idx-1 for idx in [1,2,3,4]]\n',
          'gauss':"kwargs = {'addsec':'1 2 3 4 F'}\n",
          'qchem':"kwargs = {'addsec':'\\n$opt\\nCONSTRAINT\\ntors 1 2 3 4 30.0\\nENDCONSTRAINT\\n$end'}\n"}[backend]
    path=Path(job+'.py');path.write_text(text)
    QuantumChemistry._check_hir_definition(job,[1,2,3,4],True)
    QuantumChemistry._check_hir_definition(job,[4,3,2,1],True)
    with pytest.raises(ValueError,match='Run fresh scans'):
        QuantumChemistry._check_hir_definition(job,[5,2,3,4],True)
    assert path.read_text()==text
    with pytest.raises(ValueError,match='no saved scan input'):
        QuantumChemistry._check_hir_definition(str(tmp_path/'missing'),[1,2,3,4],True)



def test_qc_hir_checks_definition_before_rewriting_or_reusing(tmp_path, monkeypatch):
    import json
    from unittest.mock import Mock
    from ase.build import molecule
    from kinbot.parameters import Parameters
    from kinbot.stationary_pt import StationaryPoint
    import kinbot.qc as module
    monkeypatch.chdir(tmp_path)
    Path('hir').mkdir()
    job = 'hir/ts_hir_0_00'
    text = 'base_0_fix = [idx-1 for idx in [1,2,3,4]]\n'
    script = Path(job+'.py')
    script.write_text(text)
    Path('input.json').write_text(json.dumps(dict(barrier_threshold=100, smiles='CO',
        rotor_scan=0, conformer_search=0, queuing='local')))
    qc = QuantumChemistry(Parameters('input.json', show_warnings=False).par)
    qc.check_qc = Mock(return_value='not found')
    qc.submit_qc = Mock()
    monkeypatch.setattr(module, 'routed_qc_job', lambda *args: job)
    atoms = molecule('CH3OH')
    point = StationaryPoint('ts', 0, 1, atom=atoms.get_chemical_symbols(), geom=atoms.positions)
    point.characterize()
    point.wellorts = 1
    replacement = qc.qc_hir(point, point.geom, 0, 0, [[5,1,2,4]], False)
    assert replacement.startswith(job + '_recovery_')
    assert script.read_text() == text
    assert Path(replacement + '.py').exists()
    qc.check_qc.return_value = 'error'
    qc.submit_qc.reset_mock()
    assert qc.qc_hir(point, point.geom, 0, 0, [[5,1,2,4]], False) == replacement
    qc.submit_qc.assert_not_called()  # A failed replacement is terminal.
