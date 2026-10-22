"""Run MESS's statistical two-transition-state model without explicit wells.

MESS's Variational input accepts one inner RRHO and one outer RRHO, not
conformer Unions. Ask MESS for the conformer sums and supply those states
through its existing Rotd table reader. MESS performs
the two-surface combination, tunneling, thermal integration and translation.
Pathways sharing the chosen outer surface share its capture flux once.
"""
import argparse
import json
import math
from pathlib import Path
import re
import shlex
import subprocess

import numpy as np

from kinbot import constants


def variational(block, *, outer=None):
    """Move Eckart to the Variational level, as required by MESS."""
    match = re.search(r'^\s*Tunneling\s+Eckart\b.*?^\s*End[^\n]*\n', block, re.M | re.S)
    tunnel = match.group() if match else ''
    if match:
        block = block[:match.start()] + block[match.end():]
    text = 'Variational\n'
    text += block + '\n' + tunnel
    if outer is not None:
        text += '\n2TSMethod statistical\n' + re.sub(r'^\s*RRHO\s*$', 'OuterRRHO', outer,
                                                     count=1, flags=re.M)
    return text + '\nEnd ! Variational\n'


def table_model(filename, mass, ground=None):
    """The table already includes all rotations, vibrations and degeneracy."""
    ground_line = '' if ground is None else f'GroundEnergy[1/cm] {ground:.12g}\n'
    energy = ('ElectronicEnergy[1/cm] 0' if ground is None
              else f'ZeroEnergy[1/cm] {ground:.12g}')
    return (f'RRHO\nMass[amu] {mass:.12g}\n+++\nCore Rotd\nFile {filename}\n'
            f'Units 1/cm\nNosTolerance 1e-300\n{ground_line}End\n'
            f'ElectronicLevels[1/cm] 1\n0 1\n{energy}\nEnd ! RRHO\n')


def header(stem, temperatures, limit):
    return ('TemperatureList[K] ' + ' '.join(map(str, temperatures)) + '\n'
            f'ModelEnergyLimit[kcal/mol] {limit}\nAtomDistanceMin[angstrom] 0.5\n'
            f'LogOutput {stem}.log\nRateOutput {stem}.out\n')


def state_input(stem, models, mass, grid, limit):
    """Use MESS's model-only mode to print states, without solving a network.

    MESS requires a well to attach barriers. The core-less reference below
    deliberately selects its model-only mode. Its states are never used in
    any rate. No RC/PC model or collisional calculation is involved.
    """
    cm_to_kcal = constants.AUtoKCAL / constants.AUtoCM
    text = (f'AtomDistanceMin[angstrom] 0.5\nTemperatureList[K] 300\nPressureList[torr] 760\n'
            f'EnergyStepOverTemperature 0.2\nExcessEnergyOverTemperature 30\n'
            f'ModelEnergyLimit[kcal/mol] {limit}\nLogOutput {stem}.log\n'
            f'MicroRateOutput {stem}.micro\n'
            f'MicroEnerMin[kcal/mol] {grid[0] * cm_to_kcal:.14g}\n'
            f'MicroEnerMax[kcal/mol] {grid[1] * cm_to_kcal:.14g}\n'
            f'MicroEnerStep[kcal/mol] {grid[2] * cm_to_kcal:.14g}\n'
            'Model\nEnergyRelaxation\nExponential\nFactor[1/cm] 200\nPower 0.85\nExponentCutoff 15\nEnd\n'
            f'CollisionFrequency\nLennardJones\nEpsilons[1/cm] 10 10\nSigmas[angstrom] 3 3\nMasses[amu] 4 {mass}\nEnd\n'
            '! This is a state-count request, not a physical intermediate.\n'
            '! Omitting its core instructs MESS to skip rate calculations.\n'
            f'Well state_count_only\nSpecies\nRRHO\nMass[amu] {mass}\n'
            f'ZeroEnergy[1/cm] {grid[0] - 100}\nEnd\nEnd\n')
    for i, model in enumerate(models):
        text += (f'Bimolecular sink_{i}\nDummy\n'
                 f'Barrier s{i} state_count_only sink_{i}\n' + model + '\n')
    return text + 'End ! Model\n'


def read_states(filename, grid, count):
    section = Path(filename).read_text().split('microcanonical number of states:', 1)[1]
    rows = []
    for line in section.splitlines()[2:]:
        if not line.strip():
            if rows:
                break
            continue
        rows.append([float(w) for w in line.split()])
    table = np.asarray(rows)
    if table.ndim != 2 or table.shape[1] != count + 1 or not np.all(np.isfinite(table)):
        raise ValueError(f'{filename}: invalid MESS state counts.')
    # MESS prints energies with six significant digits. Reconstruct the
    # requested regular grid, rather than using rounded energies in sums.
    energies = grid[0] + np.arange(len(table)) * grid[2]
    if (len(table) < 4 or abs(energies[-1] - grid[1]) > grid[2] * 1.1
            or not np.allclose(table[:, 0] * constants.AUtoCM / constants.AUtoKCAL,
                               energies, rtol=1e-5, atol=1e-4)):
        raise ValueError(f'{filename}: MESS energy grid differs from the request.')
    return energies, table[:, 1:]


def write_model(directory, stem, branches, outer, mass, lower_cm1, limit, temperatures, note, outer_threshold,
                energy_reference):
    from kinbot.bimol_mess import union
    grid = [math.floor(lower_cm1) - 2,
            math.floor(limit * constants.AUtoCM / constants.AUtoKCAL) - 2, 1]
    if grid[1] <= grid[0] or not temperatures or min(temperatures) <= 0:
        raise ValueError('Invalid energy or temperature range for the two-transition-state model.')
    prefix = ''.join('! ' + line + '\n' for line in note.splitlines())
    models = [union(outer)] + [union(b['blocks']) for b in branches]
    (directory / f'{stem}_states.inp').write_text(prefix + state_input(stem + '_states', models, mass, grid, limit))
    text = prefix + state_input(stem + '_combined', [variational(
        table_model(stem + '_inner.dat', mass), outer=table_model(stem + '_outer.dat', mass))], mass, grid, limit)
    (directory / f'{stem}.inp').write_text(text)
    (directory / f'{stem}_states.json').write_text(json.dumps(dict(stem=stem, grid=grid, mass=mass,
        temperatures=temperatures, limit=limit, note=note, outer_threshold_cm1=outer_threshold,
        energy_reference_hartree=energy_reference,
        branches=[{k: v for k, v in b.items() if k != 'blocks'} for b in branches]), indent=2) + '\n')


def save_table(filename, energies, states):
    if not np.all(np.isfinite(states)) or np.any(states < 0):
        raise ValueError(f'{filename}: MESS returned nonfinite or negative state counts.')
    positive = np.flatnonzero(states > 0)
    if len(positive) < 3 or positive[0] == 0:
        raise ValueError(f'{filename}: energy grid must include the threshold and at least three positive values.')
    first = positive[0] - 1
    if np.any(states[first + 1:] <= 0):
        raise ValueError(f'{filename}: state counts have an unexpected gap above threshold.')
    np.savetxt(filename, np.column_stack((energies[first:], states[first:])), fmt='%.12g')
    return float(energies[first])


def read_rates(filename, temperatures):
    rows = []
    for line in Path(filename).read_text().splitlines():
        words = line.split()
        if len(words) != 3:
            continue
        try:
            rows.append([float(w) for w in words])
        except ValueError:
            continue
    values = np.asarray(rows)
    if (values.shape != (len(temperatures), 3) or not np.all(np.isfinite(values))
            or np.any(values[:, 1:] < 0) or not np.allclose(values[:, 0], temperatures)):
        raise ValueError(f'{filename}: incomplete or invalid forward/reverse MESS rates.')
    return values


def run_model(input_file, command, state_command):
    """Run on the same local host or queue worker as the ordinary MESS job."""
    directory = Path(input_file).resolve().parent
    stem = Path(input_file).stem
    config = json.loads((directory / f'{stem}_states.json').read_text())

    def run(filename, executable=command):
        output = directory / (Path(filename).stem + '.stdout')
        with output.open('w') as handle:
            subprocess.run(executable + [filename], cwd=directory, stdout=handle,
                           stderr=subprocess.STDOUT, check=True)

    run(stem + '_states.inp', state_command)
    energies, counts = read_states(directory / (stem + '_states.micro'), config['grid'], len(config['branches']) + 1)
    # At an exact threshold MESS may print a tiny negative interpolation
    # residual (e.g. -1e-9 states at E=0). The channel is closed there.
    counts[energies <= config['outer_threshold_cm1'], 0] = 0.
    # A reaction cannot reach a separated endpoint below its ground energy,
    # even when tunneling and the selected outer surface extend below it.
    inner = []
    for i, branch in enumerate(config['branches']):
        values = counts[:, i + 1].copy()
        values[energies <= branch['threshold_cm1']] = 0.
        if np.any(values < 0):
            raise ValueError('MESS returned negative inner states above the separated-endpoint threshold.')
        inner.append(values)
    inner_total = np.sum(inner, axis=0)
    mass = config['mass']
    models = {}
    for side, states in (('inner', inner_total), ('outer', counts[:, 0])):
        filename = f'{stem}_{side}.dat'
        ground = save_table(directory / filename, energies, states)
        models[side] = table_model(filename, mass, ground)
    prefix = ''.join('! ' + line + '\n' for line in config['note'].splitlines())
    text = prefix + state_input(stem + '_combined', [variational(models['inner'], outer=models['outer'])],
                                mass, config['grid'], config['limit'])
    (directory / f'{stem}.inp').write_text(text)
    run(stem + '.inp', state_command)
    combined_grid, effective = read_states(directory / (stem + '_combined.micro'), config['grid'], 1)
    if not np.array_equal(combined_grid, energies):
        raise ValueError('Combined MESS barrier changed the energy grid.')
    threshold = max(config['outer_threshold_cm1'], min(b['threshold_cm1'] for b in config['branches']))
    effective[energies <= threshold, 0] = 0.
    results = []
    for i, (branch, states) in enumerate(zip(config['branches'], inner)):
        # Parallel pathways share the one chosen capture surface. This splits
        # the MESS-computed total flux; it does not apply capture again.
        share = np.divide(states, inner_total, out=np.zeros_like(states), where=inner_total > 0)
        filename = f'{stem}_channel_{i:04d}'
        ground = save_table(directory / (filename + '.dat'), energies, effective[:, 0] * share)
        text = prefix + header(filename, config['temperatures'], config['limit']) + branch['fragments']
        text += f'Barrier {filename}\n' + variational(table_model(filename + '.dat', mass, ground))
        (directory / (filename + '.inp')).write_text(text)
        run(filename + '.inp')
        results.append(read_rates(directory / (filename + '.out'), config['temperatures']))
    output = ['Rate Units: cm^3 molecule^-1 s^-1',
              '! Each channel uses its own separated R and P partition functions.']
    output += [f'! channel_{i:04d}: R={b["endpoints"]["R"]} P={b["endpoints"]["P"]}'
               for i, b in enumerate(config['branches'])]
    output.append('T[K] ' + ' '.join(f'channel_{i:04d}_forward channel_{i:04d}_reverse' for i in range(len(results))))
    for row, temperature in enumerate(config['temperatures']):
        output.append(str(temperature) + ' ' + ' '.join(f'{value:.9g}' for result in results for value in result[row, 1:]))
    (directory / (stem + '.out')).write_text('\n'.join(output) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mess-command', default='mess-abstraction')
    parser.add_argument('--state-command', default='mess')
    parser.add_argument('input')
    args = parser.parse_args()
    run_model(args.input, shlex.split(args.mess_command), shlex.split(args.state_command))


if __name__ == '__main__':
    main()
