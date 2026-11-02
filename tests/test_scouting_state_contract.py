"""Cheap scouting changes theory/cost, not an already specified fragment spin."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from kinbot.electronic_states import apply_electronic_state_model, mess_electronic_levels
from kinbot.parameters import Parameters
from kinbot.stationary_pt import StationaryPoint
from kinbot.species_routing import routing_name
from kinbot.theory import automatic_cbh_anl_requested


def test_scout_preserves_declared_diatomic_state_through_job_and_mess_identity():
    root=Path(__file__).resolve().parents[1]/'examples/anl'
    production=json.loads((root/'ch3ss_anl1_f12/ch3ss.json').read_text())
    par=Parameters(str(root/'ch3ss_anl1_f12/scouting_r8/input.json'),show_warnings=False).par
    declared=next(x for x in production['electronic_state_models'] if x['smiles']=='[S][S]')
    assert par['electronic_state_models']==production['electronic_state_models']
    # The ordinary closed-shell Lewis representation initially suggests spin 1.
    # Apply the explicit contract at the same point used for discovered products.
    fragment=StationaryPoint('product',0,1,smiles='S=S')
    fragment.characterize()
    old=routing_name(fragment)
    apply_electronic_state_model(fragment,par)
    assert fragment.mult==3 and routing_name(fragment)!=old
    assert str(fragment.chemid).endswith('3')
    assert mess_electronic_levels(fragment,par)==(1,'0 3')
    assert not automatic_cbh_anl_requested(par) and not par['L3_calc']
    assert not par['rotdpy_run'] and not par['vrc_tst_scan'] and not par['vrc_tst_noscan']


def test_both_scouts_remain_input_only_pst_with_existing_cutoffs():
    root=Path(__file__).resolve().parents[1]/'examples/anl'
    for case in ('ch3s','ch3ss'):
        par=Parameters(str(root/(case+'_anl1_f12/scouting_r8/input.json')),show_warnings=False).par
        assert par['me']==1 and par['run_me']==0 and par['rotdpy_mess_mode']=='off'
        assert par['l1_pes_energy_ceiling']==150 and par['vdW_detection']==.5
        assert par['l1_vdw_depth_check'] and par['correct_submerged']==1
        production=json.loads((root/(case+'_anl1_f12')/(case+'.json')).read_text())
        assert par['electronic_state_models']==production['electronic_state_models']
        sulfur=StationaryPoint('S',0,3,smiles='[S]')
        sulfur.characterize()
        assert mess_electronic_levels(sulfur,par)==(3,'0 5\n396.055 3\n573.64 1')


@pytest.mark.parametrize('case', ['ch3s', 'ch3ss'])
def test_all_flagged_degenerate_fragments_have_explicit_models(case):
    from kinbot.electronic_states import audit_electronic_states, spin_orbit_energy_shift
    from kinbot import constants
    root = Path(__file__).resolve().parents[1]/'examples/anl'
    par = json.loads((root/(case+'_anl1_f12/scouting_r8/input.json')).read_text())
    states = [dict(smiles=s, charge=0, multiplicity=m)
              for s, m in [('[S]', 3), ('[SH]', 2), ('[CH3][S]', 2), ('[CH]', 2)]]
    report = audit_electronic_states(states, par)
    assert not report['requires_review']
    assert all(row['status'] == 'configured_correction' for row in report['states'])
    ch = SimpleNamespace(smiles='[H][C]', charge=0, mult=2)
    assert mess_electronic_levels(ch, par) == (2, '0 2\n28.14675 2')
    assert spin_orbit_energy_shift(ch, par) == pytest.approx(
        -28.14675/2/constants.AUtoCM)
    # Fragment orbital multiplicity must not multiply a reactive doublet.
    assert mess_electronic_levels(ch, par, ground_only=True) == (1, '0 2')


def test_cached_molecular_fragment_adopts_ch_ladder_without_qc_changes():
    from kinbot.mess_cached_species import refresh_cached_population_models
    root = Path(__file__).resolve().parents[1]/'examples/anl'
    par = json.loads((root/'ch3s_anl1_f12/scouting_r8/input.json').read_text())
    text = '''Reactant w
Well w ! 42 ! C
  RRHO
  End
Bimolecular p ! 130130000000000000002_1
  Fragment f ! 130130000000000000002 ! [H][C]
    RRHO
      Geometry[angstrom] 2
        C 0 0 0
        H 0 0 1.12
      Frequencies[1/cm] 1
        2874.5
      ElectronicLevels[1/cm] 1
        0 2
      ZeroEnergy[kcal/mol] 0
    End
  GroundEnergy[kcal/mol] 120
End
End ! end kinetics
'''
    updated = refresh_cached_population_models(text, par)
    assert '28.14675 2\n' in updated
    assert 'H 0 0 1.12' in updated and '2874.5' in updated
    assert 'GroundEnergy[kcal/mol] 120' in updated
    assert refresh_cached_population_models(updated, par) == updated
