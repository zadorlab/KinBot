"""Conformers in one rate model must match its completed IRC products."""
import numpy as np

from kinbot import bimol
from kinbot.bimol_rate_models import ensemble
from kinbot.conformer_records import ConformerRecord
from kinbot.reaction_path import path_geometry_allowed
from test_bimol_stereo_integration import donor_and_hydrogen, abstraction


def test_rate_ensemble_filters_product_configuration_allowed_during_refinement():
    species = donor_and_hydrogen('CC[C@H](O)C')
    hydrogens = [atom for atom, symbol in enumerate(species.atom)
                 if symbol == 'H' and species.bond[1, atom]]
    assert len(hydrogens) == 2
    reaction = abstraction(species, donor=1, hydrogen=hydrogens[0])
    ts = bimol.build_bimolecular_ts(reaction, reaction.qc)
    bimol.update_ts_topology_from_irc(reaction)
    swapped = ts.geom.copy()
    swapped[hydrogens] = swapped[hydrogens[::-1]]
    assert path_geometry_allowed(ts, swapped)
    records = []
    for index, geometry in enumerate((ts.geom, swapped)):
        records.append(ConformerRecord(f'conformer_{index}', index,
            f'calculation_{index}', 'valid', geometry=tuple(map(tuple, geometry)),
            electronic_energy_hartree=ts.energy + index * .0001,
            zpe_hartree=ts.zpe,
            zero_energy_hartree=ts.energy + ts.zpe + index * .0001,
            frequencies_cm1=tuple(ts.freq), retained=True))
    ts.conformer_inventory = tuple(records)

    selected = ensemble([ts], 'specified')

    assert selected.bimol_reactant_only is False
    assert selected.conformer_index == [0]
    np.testing.assert_allclose(selected.conformer_geom[0], ts.geom)
    assert selected.conformer_inventory[0].retained
    assert not selected.conformer_inventory[1].retained
    assert selected.conformer_inventory[1].exclusion_reason == 'different stereochemical pathway'
    assert ts.bimol_reactant_only is True
    assert ts.conformer_inventory == tuple(records)
