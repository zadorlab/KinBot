"""Native integrator/parser smoke tests, not sulfur kinetics or Slurm tests.

Set KINBOT_TEST_MESS to a real MESS executable to exercise its two-TS parser.
The analytic zero-potential two-atom model checks RotdPy's phase-space measure.
"""
import os
from pathlib import Path
import subprocess

import numpy as np
import pytest

from kinbot.vrc_multistate import render_mess
from kinbot.rotdpy_mess import format_rotd_rrho


def free_pair_table(radius):
    from rotd_py.fragment.monoatomic import Monoatomic
    from rotd_py.sample.multi_sample import MultiSample
    from rotd_py.system import Surface
    from rotd_py.flux.flux import MultiFlux
    from rotd_py.flux.fluxbase import FluxBase
    from rotd_py.analysis import get_Ne_from_e
    from rotd_py import rotd_math as units

    class FreePair(MultiSample):
        def get_energies(self, calculator, face_id=0, flux_id=0):
            return np.zeros(1)

    fragments = [Monoatomic('H', positions=[[0, 0, 0]]) for _ in range(2)]
    surface = Surface(pivotpoints={'0': [[0, 0, 0]], '1': [[0, 0, 0]]},
                      distances=np.array([[radius]]))
    surface.surf_id = 0
    sample = FreePair('analytic_free_pair', fragments, surface)
    energies = np.linspace(0., 60000., 601)  # Kelvin, as in real ROTD inputs
    flux = MultiFlux(FluxBase(np.array([300., 600.]), energies, np.array([0.]),
        'MICROCANONICAL', dict(flux_rel_err=1., pot_smp_min=32, pot_smp_max=32,
                              tot_smp_min=32, tot_smp_max=32, smp_len=32)),
                     sample=sample, calculator={}).flux_array[0]
    flux.run(32, 0)
    normalized = flux.normalize()
    states = get_Ne_from_e(normalized.e_sum[:, 0], fragments)
    reduced_mass = fragments[0].get_total_mass()/2
    exact = 2*reduced_mass*radius**2*energies*units.Kelv
    assert flux.acct_smp() == 32 and flux.fail_smp() == 0
    assert states == pytest.approx(exact, rel=1e-12, abs=1e-12)
    return np.column_stack([energies*units.Kelv*units.Hartree2cm, states])


def test_actual_rotd_microcanonical_integral_has_analytic_normalization():
    first = free_pair_table(4.)
    second = free_pair_table(6.)
    assert second[1:, 1]/first[1:, 1] == pytest.approx(2.25)


def test_actual_messpf_counts_declared_atomic_levels_once(tmp_path):
    executable = os.environ.get('KINBOT_TEST_MESSPF')
    if not executable:
        pytest.skip('Set KINBOT_TEST_MESSPF to run the native atomic partition check')
    from kinbot.mess import MESS
    from kinbot.anl.thermochemistry import parse_messpf_output
    from ase import units
    temperatures = [298.15, 600., 1500.]
    levels = [{'energy_cm_inverse': e, 'degeneracy': g}
              for e, g in [(0., 5), (396.055, 3), (573.64, 1)]]
    outputs = []
    for name, selected in [('fine', levels), ('ground', levels[:1])]:
        inp = tmp_path/(name+'.inp')
        inp.write_text(MESS._atomic_partition_function_input('sulfur', 'S', selected, temperatures))
        result = subprocess.run([executable, inp.name], cwd=tmp_path,
                                text=True, capture_output=True, timeout=60)
        (tmp_path/(name+'.stdout')).write_text(result.stdout+'\n'+result.stderr)
        assert result.returncode == 0, result.stderr
        outputs.append(parse_messpf_output(tmp_path/(name+'.dat')))
    for fine, ground, temp in zip(*outputs, temperatures):
        expected = sum(level['degeneracy']*np.exp(-level['energy_cm_inverse']*units.invcm/(units.kB*temp))
                       for level in levels)/5
        assert np.exp(fine.ln_q-ground.ln_q) == pytest.approx(expected, rel=1e-4)


@pytest.mark.parametrize('smiles,symbols,positions,frequencies,gap', [
    ('[SH]', ['S', 'H'], [[0., 0., 0.], [1.34, 0., 0.]], [2600.], 376.835),
    ('C[S]', ['C', 'S', 'H', 'H', 'H'],
     [[0., 0., 0.], [0., 0., 1.8], [1., 0., -.3],
      [-.5, .866, -.3], [-.5, -.866, -.3]], [500.] * 9, 255.),
])
def test_actual_messpf_molecular_levels_do_not_repeat_so_energy_shift(
        tmp_path, smiles, symbols, positions, frequencies, gap):
    import json
    from copy import deepcopy
    from ase import units
    from kinbot import kb_path
    from kinbot.mess import MESS
    from kinbot.anl.thermochemistry import parse_messpf_output
    executable = os.environ.get('KINBOT_TEST_MESSPF')
    if not executable:
        pytest.skip('Set KINBOT_TEST_MESSPF to run native molecular partition checks')
    writer = MESS.__new__(MESS)
    writer.par = json.loads((Path(__file__).parents[1] /
        'examples/anl/ch3s_anl1_f12/ch3s.json').read_text())
    writer.par['TemperatureList'] = [298.15, 600., 1500.]
    writer.par.update(freq_uq_ref=100., freq_uq_max_exp=4.)
    writer.pftpl = Path(kb_path, 'tpl/mess_pf.tpl').read_text()
    model = dict(smiles=smiles, name='test', charge=0, multiplicity=2,
                 symbols=symbols, positions=positions, symmetry_factor=1.,
                 frequencies_cm_inverse=frequencies,
                 electronic_levels=[dict(energy_cm_inverse=0., degeneracy=2)],
                 source={'record_sha256': 'unchanged scalar QC'})
    original = deepcopy(model)
    outputs = []
    for name in ('fine', 'ground'):
        if name == 'ground':
            for item in writer.par['electronic_state_models']:
                item['levels'] = item['levels'][:1]
                item.pop('spin_free_reference', None)
        inp = writer.write_partition_function_model(model, tmp_path / (name + '.inp'))
        result = subprocess.run([executable, inp.name], cwd=tmp_path,
                                text=True, capture_output=True, timeout=60)
        assert result.returncode == 0, result.stderr
        outputs.append(parse_messpf_output(tmp_path / (name + '.dat')))
    for fine, ground in zip(*outputs):
        expected = 1 + np.exp(-gap * units.invcm / (units.kB * fine.temperature_k))
        assert np.exp(fine.ln_q - ground.ln_q) == pytest.approx(expected, rel=1e-4)
    assert model == original


@pytest.mark.parametrize('combination', ['single', 'statistical', 'minimum'])
def test_actual_mess_accepts_kinbot_rotd_and_existing_two_ts_blocks(tmp_path, combination):
    executable = os.environ.get('KINBOT_TEST_MESS')
    if not executable:
        pytest.skip('Set KINBOT_TEST_MESS to run the native MESS smoke test')
    source = free_pair_table(4.)
    np.savetxt(tmp_path/'Ne.dat', source)
    template = (Path(__file__).parents[1]/'kinbot/tpl/mess_rotd_rrho.tpl').read_text()
    def rrho(stage):
        return format_rotd_rrho(template, stoich='H2', file='Ne.dat', symmetry=1,
            nfreq=0, freq='', hinderedrotor='', nelec=1, electronic_levels='          0 1',
            zeroenergy=0., provenance_sha256='analytic-test')
    stages = ['single'] if combination == 'single' else ['inner', 'outer']
    barrier = render_mess({'combination': combination,
                          'staged_bottlenecks': [{'name': n} for n in stages]}, rrho)
    inp = '''TemperatureList[K] 300 600
PressureList[bar] 1
EnergyStepOverTemperature 0.5
ExcessEnergyOverTemperature 20
ModelEnergyLimit[kcal/mol] 100
CalculationMethod direct
ChemicalEigenvalueMax 0.2
MicroRateOutput micro.out
MicroEnerMin[kcal/mol] 0.1
MicroEnerMax[kcal/mol] 10
MicroEnerStep[kcal/mol] 0.5
Model
 EnergyRelaxation
  Exponential
   Factor[1/cm] 100
   Power 0.85
   ExponentCutoff 15
  End
 CollisionFrequency
  LennardJones
   Epsilons[1/cm] 100 100
   Sigmas[angstrom] 3 3
   Masses[amu] 2 28
  End
 Well W
  Species
   RRHO
    Geometry[angstrom] 2
     H 0 0 0
     H 0 0 1.0
    Core RigidRotor
     SymmetryFactor 2
    End
    Frequencies[1/cm] 1
     1000
    ZeroEnergy[kcal/mol] -5
    ElectronicLevels[1/cm] 1
     0 1
   End
 End
 Bimolecular P
  Fragment H1
   Atom
    Name H
    ElectronicLevels[1/cm] 1
     0 2
   End
  Fragment H2
   Atom
    Name H
    ElectronicLevels[1/cm] 1
     0 2
   End
  GroundEnergy[kcal/mol] 0
 End
 Barrier B W P
'''+barrier+'End\n'
    (tmp_path/'smoke.inp').write_text(inp)
    result = subprocess.run([executable, 'smoke.inp'], cwd=tmp_path,
                            text=True, capture_output=True, timeout=60)
    (tmp_path/'stdout.txt').write_text(result.stdout+'\n'+result.stderr)
    assert result.returncode == 0, result.stderr[-4000:]
    output = (tmp_path/'smoke.out').read_text()
    assert 'High Pressure Rate Coefficients' in output
    assert 'nan' not in output.lower()
    assert (tmp_path/'micro.out').stat().st_size > 0
    # Check the actual native combination, not just acceptance of keywords.
    # Equal bottlenecks give N/2 statistically and N for the minimum rule.
    import re
    from ase import units
    lines = (tmp_path/'micro.out').read_text().splitlines()
    values = []
    for line in lines[2:]:
        if not re.fullmatch(r'\s*[\d.Ee+-]+\s+[\d.Ee+-]+\s*', line):
            break
        values.append([float(x) for x in line.split()])
    values = np.asarray(values)
    assert len(values) >= 10
    expected = np.interp(values[:, 0]*units.kcal/units.mol/units.invcm, source[:, 0], source[:, 1])
    if combination == 'statistical':
        expected /= 2
    assert values[:, 1] == pytest.approx(expected, rel=2e-5)
