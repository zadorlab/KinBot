"""Scalar reference shifts agree between worker and PES-wide handoffs."""
from types import SimpleNamespace
from pathlib import Path
import numpy as np
import pytest
from ase import Atoms
from ase.db import connect
from kinbot import constants
from kinbot.energy import species_zero_k_hartree, ZeroKEnergy
from kinbot.electronic_states import apply_electronic_state_model
from kinbot.stationary_pt import StationaryPoint
from kinbot.pes import get_energy
from kinbot.mess_cached_species import refresh_cached_population_models


def parameters():
    return dict(electronic_state_models=[dict(smiles='[S]', charge=0,
        ground_multiplicity=3, source='declared test term',
        spin_orbit_correction_hartree=-.001,
        levels=[dict(energy_cm_inverse=0., multiplicity=3, degeneracy=5),
                dict(energy_cm_inverse=400., multiplicity=3, degeneracy=3)])])


def test_worker_energy_shift_is_not_mutated_or_double_applied():
    p=SimpleNamespace(name='atom', smiles='[S]', charge=0, mult=3,
                      energy=-398., zpe=0.)
    apply_electronic_state_model(p,parameters())
    for _ in range(2):
        assert species_zero_k_hartree(p) == pytest.approx(-398.001)
        assert p.energy == -398.
    p.final_zero_k_energy=ZeroKEnergy('[S]',-398.002,'ANL','native-evidence',0,3)
    assert species_zero_k_hartree(p) == -398.002
    del p.final_zero_k_energy
    p.wellorts=1
    assert species_zero_k_hartree(p) == -398.


def test_pes_cached_atomic_levels_and_energy_use_computed_spin(tmp_path,monkeypatch):
    from kinbot.run_format import ensure_current_run
    monkeypatch.chdir(tmp_path)
    ensure_current_run(create=True)
    point=StationaryPoint('atom',0,3,atom=['S'],geom=np.zeros((1,3)))
    point.characterize();key=str(point.chemid)
    ensure_current_run(key,create=True)
    db=connect(str(Path(key)/'kinbot.db'))
    db.write(Atoms('S'),name=key+'_well',data=dict(status='normal',
        energy=-398./constants.EVtoHARTREE,zpe=0.,frequencies=[]))
    before=(Path(key)/'kinbot.db').read_bytes()
    shifts={};states={}
    for _ in range(2):
        assert get_energy([key],key,0,0,electronic_state_parameters=parameters(),
            correction_cache=shifts,state_cache=states)==pytest.approx((-398.001,0.))
    assert shifts[key]==-.001 and states[key]['multiplicity']==3
    assert (Path(key)/'kinbot.db').read_bytes()==before
    text=f'''Reactant w
Well w ! 42 ! C
  RRHO
  End
Bimolecular p ! {key}_1
  Fragment f ! {key}
    Atom
      Name S
      ElectronicLevels[1/cm] 1
        0 3
    End
  GroundEnergy[kcal/mol] 1
End
End ! end kinetics
'''
    updated=refresh_cached_population_models(text,parameters(),selected_states=states)
    assert '0 5\n' in updated and '400 3\n' in updated
    assert refresh_cached_population_models(updated,parameters(),selected_states=states)==updated
    # A symbol alone cannot distinguish an ion from a neutral fragment.
    assert refresh_cached_population_models(text,parameters())==text
    par=parameters();par['electronic_state_models'][0]['ground_multiplicity']=1
    for level in par['electronic_state_models'][0]['levels']:level['multiplicity']=1
    with pytest.raises(ValueError,match='different multiplicity'):
        get_energy([key],key,0,0,electronic_state_parameters=par)
    with pytest.raises(ValueError,match='cached QC spin'):
        refresh_cached_population_models(text,par,selected_states=states)
