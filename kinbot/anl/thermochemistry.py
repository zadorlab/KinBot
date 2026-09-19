"""MESS partition functions and the 0 K to 298.15 K formation conversion."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import re
from typing import Mapping

import numpy as np

from kinbot.energy import formation_enthalpy_0k_kj_mol


R_KJ_MOL_K = 0.00831446261815324
R_J_MOL_K = R_KJ_MOL_K * 1000.
CAL_TO_J = 4.184

# H(298.15 K) - H(0 K), per mole of atoms in the elemental reference state.
# Values are from the pinned NIST-JANAF reference-state tables. Molecular
# reference values are divided by two. Carbon is the graphite reference state.
ELEMENT_REFERENCE_HEAT_CONTENT_298_KJ_MOL_ATOM = {
    'H': 8.467 / 2.,
    'C': 1.051,
    'N': 8.670 / 2.,
    'O': 8.683 / 2.,
    'F': 8.825 / 2.,
    'Cl': 9.181 / 2.,
}
ELEMENT_REFERENCE_SOURCE = {
    'H': 'https://janaf.nist.gov/tables/H-050.html',
    'C': 'https://janaf.nist.gov/tables/C-002.html',
    'N': 'https://janaf.nist.gov/tables/N-023.html',
    'O': 'https://janaf.nist.gov/tables/O-029.html',
    'F': 'https://janaf.nist.gov/tables/F-054.html',
    'Cl': 'https://janaf.nist.gov/tables/Cl-073.html',
}


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


def formation_enthalpy_298_kj_mol(formation_0k_kj_mol: float,
                                  species_heat_content_kj_mol: float,
                                  atoms,
                                  elemental_heat_content: Mapping[str, float] | None = None):
    """Convert Hf(0 K) to Hf(298.15 K) using elemental reference states."""
    elemental = (ELEMENT_REFERENCE_HEAT_CONTENT_298_KJ_MOL_ATOM
                 if elemental_heat_content is None else elemental_heat_content)
    counts = Counter(atoms)
    missing = sorted(set(counts) - set(elemental))
    if missing:
        raise ValueError('No 298.15 K elemental reference heat content for: '
                         + ', '.join(missing))
    values = (formation_0k_kj_mol, species_heat_content_kj_mol,
              *(count * elemental[element] for element, count in counts.items()))
    if not all(math.isfinite(value) for value in values):
        raise ValueError('Formation-enthalpy conversion contains a nonfinite value.')
    element_increment = math.fsum(
        count * elemental[element] for element, count in counts.items())
    return formation_0k_kj_mol + species_heat_content_kj_mol - element_increment


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


def messpf_thermochemistry_record(species, output, temperature_k=298.15):
    """Build the auditable datum needed to anchor a NASA/PAC99 fit."""
    points = parse_messpf_output(output)
    point = point_at_temperature(points, temperature_k)
    formation_0k = formation_enthalpy_0k_kj_mol(species)
    formation_298 = formation_enthalpy_298_kj_mol(
        formation_0k, point.heat_content_0k_kj_mol, species.atom)
    nasa7 = fit_nasa7(points, formation_298,
                      reference_temperature_k=temperature_k)
    return {
        'schema': 1,
        'species': {'name': species.name, 'chemid': str(species.chemid),
                    'smiles': species.smiles, 'atoms': dict(Counter(species.atom))},
        'temperature_k': point.temperature_k,
        'formation_0k_kj_mol': formation_0k,
        'species_heat_content_0_to_t_kj_mol': point.heat_content_0k_kj_mol,
        'element_reference_heat_content_0_to_298_15_kj_mol_atom':
            dict(ELEMENT_REFERENCE_HEAT_CONTENT_298_KJ_MOL_ATOM),
        'element_reference_sources': dict(ELEMENT_REFERENCE_SOURCE),
        'formation_298k_kj_mol': formation_298,
        'messpf': asdict(point),
        'messpf_temperature_grid': [asdict(item) for item in points],
        'nasa7': nasa7,
    }


def write_thermochemistry_record(species, output, destination, temperature_k=298.15):
    record = messpf_thermochemistry_record(species, output, temperature_k)
    Path(destination).write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    return record
