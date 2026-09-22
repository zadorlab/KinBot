"""Immutable conformer observations alongside the legacy parallel arrays.

Indices identify calculations, not positions in a filtered list. Energies in
the conformer search are E + ZPE; electronic energy is never guessed from it.
"""
from dataclasses import asdict, dataclass, replace
import logging
import subprocess
import numpy as np

from kinbot.calculation import array_fingerprint


def hessian_record(qc, job, geometry, atoms, *, row=None):
    """Read an existing Hessian; a supplied row uses only its stored matrix."""
    try:
        stored_only = row is not None
        if row is None:
            from kinbot.species_routing import resolve_job
            rows = list(qc.db.select(name=resolve_job(qc.db, job)))
            if not rows:
                return {}
            row = rows[-1]
        if (list(row.symbols) != list(atoms)
                or not np.allclose(row.positions, geometry, rtol=0., atol=1.e-8)):
            return {}
        hess = np.asarray(row.data.get('hess', []) if stored_only
                          else qc.read_qc_hess(job, len(atoms)), dtype=float)
        if hess.shape != (3 * len(atoms), 3 * len(atoms)) or not np.all(np.isfinite(hess)):
            return {}
        # KinBot stores calc_vibrations Cartesian Hessians in database rows.
        # Native Q-Chem's weighted matrix comes from its output parser instead;
        # the current backend setting must not relabel a stored Sella matrix.
        weighted = False if stored_only else qc.hessian_is_massweighted()
        if not isinstance(weighted, (bool, np.bool_)):
            return {}
        reference = dict(source_job=job, source_row_id=row.id,
            hessian_source_job=row.name, atoms=list(map(str, atoms)),
            geometry_sha256=array_fingerprint(geometry),
            hessian_sha256=array_fingerprint(hess), hessian_massweighted=bool(weighted),
            hessian_unit='hartree / (bohr^2 * amu)' if weighted else 'hartree / bohr^2')
        return dict(hessian=tuple(tuple(map(float, r)) for r in hess),
                    hessian_reference=reference)
    except (AttributeError, TypeError, ValueError, KeyError, OSError,
            NotImplementedError, subprocess.SubprocessError) as error:
        logging.getLogger('KinBot').debug('No stored optical Hessian for %s: %s', job, error)
        return {}


@dataclass(frozen=True)
class ConformerRecord:
    member_id: str
    index: int
    source_job: str | None
    status: str
    geometry: tuple | None = None
    zero_energy_hartree: float | None = None
    electronic_energy_hartree: float | None = None
    zpe_hartree: float | None = None
    frequencies_cm1: tuple | None = None
    retained: bool = False
    exclusion_reason: str | None = None
    sigma_ext: int | None = None
    mirror_states: int | None = None
    remaining_optical_weight: float | None = None
    population_ratio: float | None = None
    stereo_identity: str | None = None
    optical_population: str | None = None
    mirror_coverage: str | None = None
    mirror_partner_ids: tuple = ()
    duplicate_of: int | None = None
    attempted_source_job: str | None = None
    optical_evidence: dict | None = None
    hessian: tuple | None = None
    hessian_reference: dict | None = None

    def as_dict(self):
        # Counting needs the matrix in memory, not a second copy of every
        # Hessian in network/observation JSON. Its source and diagnostics remain.
        result = asdict(replace(self, hessian=None))
        result.pop('hessian')
        result['hessian_available'] = self.hessian is not None
        result.update(geometry_unit='angstrom', representation='RRHO',
                      atom_mapping='same order as species.atom')
        return result


def inventory(species, geometries, energies, frequencies, valid, sources=None):
    """Capture every result before duplicate or population filtering."""
    if valid and all(status != 0 for status in valid) and not len(geometries):
        # The legacy all-failed search returns no geometry/property arrays.
        geometries = energies = frequencies = [None] * len(valid)
    lengths = {len(geometries), len(energies), len(frequencies), len(valid)}
    if len(lengths) != 1:
        raise ValueError('Conformer geometry/property/status arrays have different lengths.')
    sources = sources or {}
    records = []
    prefix = str(getattr(species, 'name', 'conformer'))
    for index, status in enumerate(valid):
        source = sources.get(index, {})
        record = ConformerRecord(f'{prefix}:conformer:{index}', index,
                                 source.get('source_job'), 'failed')
        if status == 0:
            geom = np.asarray(geometries[index], dtype=float)
            freq = frequencies[index]
            if (geom.shape != (len(species.atom), 3)
                    or not np.all(np.isfinite(geom))
                    or not np.isfinite(energies[index])):
                raise ValueError(f'Invalid properties for successful conformer {index}.')
            record = replace(record, status='valid',
                geometry=tuple(tuple(float(x) for x in atom) for atom in geom),
                zero_energy_hartree=float(energies[index]),
                electronic_energy_hartree=source.get('electronic_energy_hartree'),
                zpe_hartree=source.get('zpe_hartree'),
                hessian=source.get('hessian'),
                hessian_reference=source.get('hessian_reference'),
                frequencies_cm1=tuple(float(x) for x in freq) if freq is not None else ())
        records.append(record)
    return tuple(records)


def retain(species, records, indices):
    """Keep the unfiltered inventory and an index-keyed retained view."""
    indices = set(indices)
    species.conformer_inventory = tuple(
        replace(record, retained=record.index in indices,
                exclusion_reason=(None if record.index in indices else
                                  record.exclusion_reason or
                                  ('not retained' if record.status == 'valid' else record.status)))
        for record in records)
    species.conformer_records = {record.index: record
                                 for record in species.conformer_inventory
                                 if record.retained}


def update_member(species, index, *, source_job=None, accepted=True, hessian_data=None):
    """Refresh a retained L2 member by its original calculation index."""
    records = getattr(species, 'conformer_records', {})
    if index not in records:
        return
    record = records[index]
    if accepted:
        offset = species.conformer_index.index(index)
        record = replace(record, source_job=source_job,
            geometry=tuple(tuple(float(x) for x in atom)
                           for atom in species.conformer_geom[offset]),
            electronic_energy_hartree=float(species.conformer_energy[offset]),
            zero_energy_hartree=float(species.conformer_zeroenergy[offset]),
            zpe_hartree=float(species.conformer_zeroenergy[offset]
                              - species.conformer_energy[offset]),
            frequencies_cm1=tuple(map(float, species.conformer_freq[offset])),
            hessian=(hessian_data or {}).get('hessian'),
            hessian_reference=(hessian_data or {}).get('hessian_reference'),
            sigma_ext=None, mirror_states=None, remaining_optical_weight=None,
            optical_evidence=None)
    else:
        record = replace(record, status='failed', retained=False,
                         attempted_source_job=source_job,
                         exclusion_reason='L2 calculation failed')
    records[index] = record
    species.conformer_inventory = tuple(
        record if member.index == index else member
        for member in species.conformer_inventory)
