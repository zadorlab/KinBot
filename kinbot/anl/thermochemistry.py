"""MESS partition functions and CBH formation enthalpies at temperature."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import numpy as np

from kinbot.energy import formation_enthalpy_0k_kj_mol


R_KJ_MOL_K = 0.00831446261815324
R_J_MOL_K = R_KJ_MOL_K * 1000.
CAL_TO_J = 4.184


def _count_hindered_rotor_blocks(text: str) -> int:
    """Count MESS ``Rotor Hindered`` blocks in a rendered species model."""
    return len(re.findall(
        r'^\s*Rotor\s+Hindered(?:\s|$)', text,
        flags=re.IGNORECASE | re.MULTILINE))

@dataclass(frozen=True)
class PartitionFunctionPoint:
    temperature_k: float
    ln_q: float
    d_ln_q_d_temperature: float
    d2_ln_q_d_temperature2: float
    entropy_cal_mol_k: float
    heat_capacity_cal_mol_k: float

    @property
    def heat_content_0k_kj_mol(self) -> float:
        """Ideal-gas H(T)-H(0), including the pV=RT term."""
        return R_KJ_MOL_K * (
            self.temperature_k ** 2 * self.d_ln_q_d_temperature
            + self.temperature_k)


def parse_messpf_output(source) -> tuple[PartitionFunctionPoint, ...]:
    """Parse the single-species table written by the official ``messpf``."""
    if hasattr(source, 'read'):
        text = source.read()
    elif isinstance(source, Path):
        text = source.read_text()
    elif isinstance(source, str) and '\n' not in source:
        text = Path(source).read_text()
    else:
        text = str(source)
    lines = text.splitlines()
    start = next((index for index, line in enumerate(lines)
                  if 'Z_0' in line and 'Z_1' in line and 'Z_2' in line), None)
    if start is None:
        raise ValueError('MESSPF output does not contain the partition-function table.')
    points = []
    number = re.compile(r'^[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+\-]?\d+)?$')
    for line in lines[start + 1:]:
        words = line.split()
        if len(words) != 6 or not all(number.match(word) for word in words):
            if points:
                break
            continue
        values = [float(word.replace('D', 'E').replace('d', 'e')) for word in words]
        if not all(math.isfinite(value) for value in values):
            raise ValueError('MESSPF output contains a nonfinite thermochemical value.')
        points.append(PartitionFunctionPoint(*values))
    if not points:
        raise ValueError('MESSPF partition-function table has no data rows.')
    return tuple(points)


def point_at_temperature(points, temperature_k=298.15, tolerance_k=0.02):
    """Select an explicitly computed temperature; do not silently interpolate."""
    matches = [point for point in points
               if abs(point.temperature_k - temperature_k) <= tolerance_k]
    if len(matches) != 1:
        raise ValueError(f'MESSPF output needs exactly one row at {temperature_k} K.')
    return matches[0]


def cbh_formation_enthalpy_at_temperature(formation, heat_contents):
    """Propagate a CBH formation tree to one common temperature.

    ``heat_contents`` maps every species in the tree to its calculated
    ideal-gas ``H(T)-H(0)``.  CBH-0 terminates in pinned ATcT values at that
    temperature.  Higher rungs recurse through lower CBH results, so no
    elemental reference-state increment enters the calculation.
    """
    from kinbot.anl.atct import ATcTRecord, _canonical_smiles
    from kinbot.anl.cbh import FormationEnthalpy

    if not isinstance(formation, FormationEnthalpy):
        raise TypeError('Expected a CBH formation record.')
    canonical = {}
    for smiles, value in heat_contents.items():
        key = _canonical_smiles(smiles)
        if key is None or key in canonical:
            raise ValueError('CBH heat-content identities are invalid or repeated.')
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f'{smiles}: nonfinite heat content.')
        canonical[key] = number

    def heat(smiles):
        try:
            return canonical[_canonical_smiles(smiles)]
        except KeyError as error:
            raise ValueError(
                f'{smiles}: missing MESSPF heat content for the CBH reaction.') from error

    def solve(node):
        if not node.stoichiometry:
            reference = node.references.get(node.target_smiles)
            if (len(node.references) != 1
                    or not isinstance(reference, ATcTRecord)
                    or reference.formation_298k_kj_mol is None
                    or not math.isfinite(reference.formation_298k_kj_mol)):
                raise ValueError('Terminal CBH-0 reference lacks ATcT Hf(298 K).')
            value = float(reference.formation_298k_kj_mol)
            return value, {
                'target_smiles': node.target_smiles, 'rung': 0,
                'kind': 'atct_terminal', 'atct_id': reference.atct_id,
                'formation_temperature_kj_mol': value,
            }
        if node.stoichiometry.get(node.target_smiles) != -1:
            raise ValueError('CBH target coefficient must be -1.')
        if set(node.references) != set(node.stoichiometry) - {node.target_smiles}:
            raise ValueError('CBH references disagree with the stored reaction.')
        reference_values = {}
        reference_trace = {}
        for smiles, reference in node.references.items():
            if isinstance(reference, ATcTRecord):
                if node.rung != 0:
                    raise ValueError('Only CBH-0 may use a direct ATcT reference.')
                value = reference.formation_298k_kj_mol
                if value is None or not math.isfinite(value):
                    raise ValueError(f'{smiles}: ATcT reference lacks Hf(298 K).')
                reference_values[smiles] = float(value)
                reference_trace[smiles] = {
                    'kind': 'atct', 'atct_id': reference.atct_id,
                    'formation_temperature_kj_mol': float(value)}
            elif isinstance(reference, FormationEnthalpy):
                value, trace = solve(reference)
                reference_values[smiles] = value
                reference_trace[smiles] = trace
            else:
                raise ValueError(f'{smiles}: invalid CBH formation reference.')
        thermal_reaction = math.fsum(
            coefficient * heat(smiles)
            for smiles, coefficient in node.stoichiometry.items())
        reaction_temperature = math.fsum((
            node.reaction_energy_0k_kj_mol, thermal_reaction))
        reference_sum = math.fsum(
            node.stoichiometry[smiles] * value
            for smiles, value in reference_values.items())
        result = reference_sum - reaction_temperature
        if not math.isfinite(result):
            raise ValueError('CBH finite-temperature formation result is nonfinite.')
        return result, {
            'target_smiles': node.target_smiles, 'rung': node.rung,
            'kind': 'cbh_reaction',
            'stoichiometry': dict(node.stoichiometry),
            'reaction_energy_0k_kj_mol': node.reaction_energy_0k_kj_mol,
            'thermal_reaction_correction_kj_mol': thermal_reaction,
            'reaction_enthalpy_temperature_kj_mol': reaction_temperature,
            'reference_formation_sum_temperature_kj_mol': reference_sum,
            'formation_temperature_kj_mol': result,
            'references': reference_trace,
        }

    return solve(formation)


def nasa7_heat_capacity_j_mol_k(coefficients, temperature_k):
    """Evaluate the NASA-7 constant-pressure heat capacity."""
    a1, a2, a3, a4, a5 = coefficients[:5]
    temperature = float(temperature_k)
    return R_J_MOL_K * (a1 + a2 * temperature + a3 * temperature ** 2
                        + a4 * temperature ** 3 + a5 * temperature ** 4)


def nasa7_enthalpy_kj_mol(coefficients, temperature_k):
    """Evaluate the NASA-7 molar enthalpy, including its absolute anchor."""
    a1, a2, a3, a4, a5, a6 = coefficients[:6]
    temperature = float(temperature_k)
    h_over_rt = (a1 + a2 * temperature / 2. + a3 * temperature ** 2 / 3.
                 + a4 * temperature ** 3 / 4.
                 + a5 * temperature ** 4 / 5. + a6 / temperature)
    return R_KJ_MOL_K * temperature * h_over_rt


def nasa7_entropy_j_mol_k(coefficients, temperature_k):
    """Evaluate the NASA-7 standard molar entropy."""
    a1, a2, a3, a4, a5, _, a7 = coefficients
    temperature = float(temperature_k)
    s_over_r = (a1 * math.log(temperature) + a2 * temperature
                + a3 * temperature ** 2 / 2. + a4 * temperature ** 3 / 3.
                + a5 * temperature ** 4 / 4. + a7)
    return R_J_MOL_K * s_over_r


def _nasa_cp_coefficients(points):
    """Fit Cp/R in T/1000 coordinates, then return standard NASA powers."""
    if len(points) < 5:
        raise ValueError('A NASA-7 Cp range needs at least five temperatures.')
    temperatures = np.asarray([point.temperature_k for point in points], dtype=float)
    cp_over_r = np.asarray(
        [point.heat_capacity_cal_mol_k * CAL_TO_J / R_J_MOL_K
         for point in points], dtype=float)
    x = temperatures / 1000.
    scaled, _, rank, _ = np.linalg.lstsq(
        np.column_stack([x ** power for power in range(5)]), cp_over_r,
        rcond=None)
    if rank != 5:
        raise ValueError('MESSPF temperature grid is singular for NASA-7 fitting.')
    return [float(value / 1000. ** power)
            for power, value in enumerate(scaled)]


def _nasa_anchor(cp_coefficients, temperature_k, enthalpy_kj_mol,
                 entropy_j_mol_k):
    a1, a2, a3, a4, a5 = cp_coefficients
    temperature = float(temperature_k)
    h_without_constant = R_KJ_MOL_K * temperature * (
        a1 + a2 * temperature / 2. + a3 * temperature ** 2 / 3.
        + a4 * temperature ** 3 / 4. + a5 * temperature ** 4 / 5.)
    a6 = (enthalpy_kj_mol - h_without_constant) / R_KJ_MOL_K
    s_without_constant = R_J_MOL_K * (
        a1 * math.log(temperature) + a2 * temperature
        + a3 * temperature ** 2 / 2. + a4 * temperature ** 3 / 3.
        + a5 * temperature ** 4 / 4.)
    a7 = (entropy_j_mol_k - s_without_constant) / R_J_MOL_K
    return [*cp_coefficients, a6, a7]


def fit_nasa7(points, formation_298_kj_mol, *, reference_temperature_k=298.15,
              transition_temperature_k=1000.):
    """Fit an H/S-continuous two-range NASA-7 model to a MESSPF table.

    Cp is fitted independently in each range.  The low-temperature integration
    constants reproduce Hf(298.15) and MESSPF S(298.15); the high-temperature
    constants enforce H and S continuity at the range boundary.
    """
    points = tuple(points)
    reference = point_at_temperature(points, reference_temperature_k)
    if not math.isfinite(formation_298_kj_mol):
        raise ValueError('NASA-7 formation-enthalpy anchor is nonfinite.')
    low_points = tuple(point for point in points
                       if point.temperature_k <= transition_temperature_k)
    high_points = tuple(point for point in points
                        if point.temperature_k >= transition_temperature_k)
    low_cp = _nasa_cp_coefficients(low_points)
    high_cp = _nasa_cp_coefficients(high_points)
    reference_entropy = reference.entropy_cal_mol_k * CAL_TO_J
    low = _nasa_anchor(low_cp, reference_temperature_k,
                       formation_298_kj_mol, reference_entropy)
    boundary_h = nasa7_enthalpy_kj_mol(low, transition_temperature_k)
    boundary_s = nasa7_entropy_j_mol_k(low, transition_temperature_k)
    high = _nasa_anchor(high_cp, transition_temperature_k,
                        boundary_h, boundary_s)
    cp_discontinuity = (nasa7_heat_capacity_j_mol_k(
        high, transition_temperature_k) - nasa7_heat_capacity_j_mol_k(
            low, transition_temperature_k))

    def rmse(data, coefficients):
        residuals = [nasa7_heat_capacity_j_mol_k(coefficients,
                                                 point.temperature_k)
                     - point.heat_capacity_cal_mol_k * CAL_TO_J
                     for point in data]
        return math.sqrt(math.fsum(value * value for value in residuals)
                         / len(residuals))

    return {
        'format': 'NASA7',
        'temperature_low_k': min(point.temperature_k for point in points),
        'temperature_mid_k': float(transition_temperature_k),
        'temperature_high_k': max(point.temperature_k for point in points),
        'low_coefficients': low,
        'high_coefficients': high,
        'enthalpy_anchor_temperature_k': float(reference_temperature_k),
        'enthalpy_anchor_kj_mol': float(formation_298_kj_mol),
        'entropy_anchor_j_mol_k': reference_entropy,
        'cp_discontinuity_at_mid_j_mol_k': cp_discontinuity,
        'low_cp_rmse_j_mol_k': rmse(low_points, low),
        'high_cp_rmse_j_mol_k': rmse(high_points, high),
    }


def messpf_thermochemistry_record(species, output, heat_contents,
                                  temperature_k=298.15):
    """Build a CBH-reaction thermochemistry and NASA/PAC99 record."""
    from dataclasses import asdict as dataclass_dict

    points = parse_messpf_output(output)
    point = point_at_temperature(points, temperature_k)
    formation_0k = formation_enthalpy_0k_kj_mol(species)
    formation = species.formation_enthalpy_0k
    formation_temperature, cbh_trace = cbh_formation_enthalpy_at_temperature(
        formation, heat_contents)
    nasa7 = fit_nasa7(points, formation_temperature,
                      reference_temperature_k=temperature_k)
    return {
        'schema': 1,
        'species': {'name': species.name, 'chemid': str(species.chemid),
                    'smiles': species.smiles, 'atoms': dict(Counter(species.atom))},
        'temperature_k': point.temperature_k,
        'formation_0k_kj_mol': formation_0k,
        'species_heat_content_0_to_t_kj_mol': point.heat_content_0k_kj_mol,
        'formation_298k_kj_mol': formation_temperature,
        'formation_temperature_method': 'cbh_reaction',
        'cbh_formation': dataclass_dict(formation),
        'cbh_temperature_trace': cbh_trace,
        'messpf': asdict(point),
        'messpf_temperature_grid': [asdict(item) for item in points],
        'nasa7': nasa7,
    }


def write_thermochemistry_record(species, output, destination, heat_contents,
                                 temperature_k=298.15):
    record = messpf_thermochemistry_record(
        species, output, heat_contents, temperature_k)
    Path(destination).write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    return record


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def audit_messpf_run(run_dir):
    """Verify native MESSPF outputs and the derived thermochemistry record.

    This audit is intentionally independent of the live KinBot objects used
    to write the inputs.  It binds every input and output to the execution
    record, reparses the native table, and checks the Hf(0)->Hf(298.15) and
    NASA-7 anchor identities stored by :class:`kinbot.mess.MESS`.  The
    finite-temperature formation enthalpy is independently reconstructed
    through the stored CBH reaction tree; elemental reference-state thermal
    increments are neither accepted nor used.
    """
    root = Path(run_dir).resolve()
    me = root if root.name == 'me' else root / 'me'
    directory = me / 'partition_functions'
    manifest_path = directory / 'manifest.json'
    execution_path = directory / 'execution.json'
    thermo_path = me / 'thermochemistry_298.json'
    for path in (manifest_path, execution_path, thermo_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    try:
        manifest = json.loads(manifest_path.read_text())
        execution = json.loads(execution_path.read_text())
        thermo = json.loads(thermo_path.read_text())
    except json.JSONDecodeError as error:
        raise ValueError('MESSPF audit record is not valid JSON.') from error
    if (manifest.get('schema') != 1 or execution.get('schema') != 1
            or execution.get('program') != 'messpf'
            or thermo.get('schema') != 1):
        raise ValueError('MESSPF audit records use an unsupported schema.')
    manifest_species = manifest.get('species')
    executions = execution.get('species')
    thermo_species = thermo.get('species')
    if (not isinstance(manifest_species, dict) or not manifest_species
            or set(manifest_species) != set(executions or {})
            or not isinstance(thermo_species, dict)):
        raise ValueError('MESSPF species coverage differs across its records.')

    network_ids = {
        ident for ident, entry in manifest_species.items()
        if entry.get('role', 'network') == 'network'}
    if set(thermo_species) != network_ids:
        raise ValueError(
            'Derived thermochemistry must cover exactly the network species.')
    model_library = manifest.get('reference_model_library')
    if model_library is not None:
        if (not isinstance(model_library, dict)
                or set(model_library) != {'path', 'sha256'}):
            raise ValueError('Invalid CBH reference-model provenance.')
        model_path = Path(model_library['path']).expanduser().resolve()
        if (not model_path.is_file()
                or _sha256(model_path) != model_library['sha256']):
            raise ValueError('CBH reference-model library hash differs.')

    checks = {}
    parsed_outputs = {}
    heat_contents = {}
    for ident, declared in sorted(manifest_species.items()):
        native = executions[ident]
        if native.get('status') != 'complete' or native.get('returncode') != 0:
            raise ValueError(f'{ident}: MESSPF execution is not complete.')
        input_path = directory / native.get('input', '')
        output_path = directory / native.get('output', '')
        for path, key in ((input_path, 'input_sha256'),
                          (output_path, 'output_sha256')):
            if (not path.is_file() or native.get(key) != _sha256(path)):
                raise ValueError(f'{ident}: MESSPF {key} verification failed.')
        input_text = input_path.read_text()
        hindered_rotor_blocks = _count_hindered_rotor_blocks(input_text)
        points = parse_messpf_output(output_path)
        point = point_at_temperature(points)
        smiles = declared.get('smiles')
        if not isinstance(smiles, str) or not smiles:
            raise ValueError(f'{ident}: MESSPF SMILES identity is missing.')
        from kinbot.anl.atct import _canonical_smiles
        canonical = _canonical_smiles(smiles)
        if canonical is None:
            raise ValueError(f'{ident}: invalid MESSPF SMILES identity.')
        previous = heat_contents.get(canonical)
        if (previous is not None
                and not math.isclose(previous, point.heat_content_0k_kj_mol,
                                     rel_tol=0., abs_tol=1e-10)):
            raise ValueError(f'{smiles}: conflicting MESSPF heat contents.')
        heat_contents[canonical] = point.heat_content_0k_kj_mol
        parsed_outputs[ident] = (points, point, hindered_rotor_blocks,
                                 native, declared)

    from kinbot.anl.cbh import formation_from_dict
    for ident in sorted(network_ids):
        points, point, hindered_rotor_blocks, native, declared = \
            parsed_outputs[ident]
        record = thermo_species[ident]
        stored_point = record.get('messpf', {})
        for key, value in asdict(point).items():
            try:
                stored = float(stored_point[key])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(
                    f'{ident}: stored MESSPF point is incomplete.') from error
            if not math.isclose(stored, value, rel_tol=0., abs_tol=1e-10):
                raise ValueError(
                    f'{ident}: stored MESSPF point differs from native output.')
        try:
            formation_0k = float(record['formation_0k_kj_mol'])
            formation_temperature = float(record['formation_298k_kj_mol'])
            heat_content = float(record[
                'species_heat_content_0_to_t_kj_mol'])
            formation = formation_from_dict(record['cbh_formation'])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                f'{ident}: thermochemistry identity is incomplete.') from error
        if (record.get('formation_temperature_method') != 'cbh_reaction'
                or not math.isclose(
                    formation_0k, formation.formation_0k_kj_mol,
                    rel_tol=0., abs_tol=1e-9)
                or not math.isclose(
                    heat_content, point.heat_content_0k_kj_mol,
                    rel_tol=0., abs_tol=1e-10)):
            raise ValueError(f'{ident}: stored CBH thermochemistry differs.')
        expected_temperature, expected_trace = \
            cbh_formation_enthalpy_at_temperature(
                formation, heat_contents)
        if not math.isclose(formation_temperature, expected_temperature,
                            rel_tol=0., abs_tol=1e-9):
            raise ValueError(
                f'{ident}: Hf(298.15 K) differs from the CBH reaction.')
        if record.get('cbh_temperature_trace') != expected_trace:
            raise ValueError(
                f'{ident}: stored finite-temperature CBH trace differs.')
        nasa = record.get('nasa7', {})
        low = nasa.get('low_coefficients')
        high = nasa.get('high_coefficients')
        middle = float(nasa.get('temperature_mid_k', 0.))
        if (nasa.get('format') != 'NASA7'
                or not isinstance(low, list) or len(low) != 7
                or not isinstance(high, list) or len(high) != 7
                or middle <= 0.):
            raise ValueError(f'{ident}: NASA-7 fit is incomplete.')
        if not math.isclose(
                nasa7_enthalpy_kj_mol(low, point.temperature_k),
                formation_temperature, rel_tol=0., abs_tol=1e-8):
            raise ValueError(f'{ident}: NASA-7 enthalpy anchor differs.')
        if not math.isclose(
                nasa7_entropy_j_mol_k(low, point.temperature_k),
                point.entropy_cal_mol_k * CAL_TO_J,
                rel_tol=0., abs_tol=1e-8):
            raise ValueError(f'{ident}: NASA-7 entropy anchor differs.')
        if (not math.isclose(nasa7_enthalpy_kj_mol(low, middle),
                             nasa7_enthalpy_kj_mol(high, middle),
                             rel_tol=0., abs_tol=1e-7)
                or not math.isclose(nasa7_entropy_j_mol_k(low, middle),
                                    nasa7_entropy_j_mol_k(high, middle),
                                    rel_tol=0., abs_tol=1e-7)):
            raise ValueError(f'{ident}: NASA-7 ranges are not H/S continuous.')
        checks[ident] = {
            'name': declared.get('name'),
            'input_sha256': native['input_sha256'],
            'output_sha256': native['output_sha256'],
            'temperature_points': len(points),
            'formation_0k_kj_mol': formation_0k,
            'formation_298k_kj_mol': formation_temperature,
            'heat_content_0_to_298_15_kj_mol': heat_content,
            'hindered_rotor_blocks': hindered_rotor_blocks,
            'nasa7_low_cp_rmse_j_mol_k': nasa.get('low_cp_rmse_j_mol_k'),
            'nasa7_high_cp_rmse_j_mol_k': nasa.get('high_cp_rmse_j_mol_k'),
            'cbh_rung': formation.rung,
        }
    reference_checks = {}
    for ident in sorted(set(manifest_species) - network_ids):
        points, point, hindered_rotor_blocks, native, declared = \
            parsed_outputs[ident]
        reference_checks[ident] = {
            'name': declared.get('name'), 'smiles': declared.get('smiles'),
            'input_sha256': native['input_sha256'],
            'output_sha256': native['output_sha256'],
            'temperature_points': len(points),
            'heat_content_0_to_298_15_kj_mol':
                point.heat_content_0k_kj_mol,
            'hindered_rotor_blocks': hindered_rotor_blocks,
        }
    return {
        'schema': 1, 'status': 'passed', 'program': 'messpf',
        'directory': str(directory), 'species': checks,
        'cbh_reference_species': reference_checks,
        'reference_model_library': model_library,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Audit native MESSPF thermochemistry and NASA fits')
    parser.add_argument('run_dir', type=Path)
    args = parser.parse_args(argv)
    print(json.dumps(audit_messpf_run(args.run_dir), indent=2,
                     sort_keys=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
