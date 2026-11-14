"""Accepted basis-specific minima carried through the L1 -> L2 lineage.

These are immutable geometry inputs, not synthetic completed ANL tasks. Native
VPT2 work still runs and is parsed normally at each basis's own minimum.
"""
from copy import copy, deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms


BASES = ('cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ')


def canonical_basis(basis):
    return next((b for b in BASES if b.casefold() == str(basis).casefold()), None)


def basis_pair(basis):
    basis = canonical_basis(basis)
    if basis is None:
        raise ValueError('Anharmonic extrapolation requires a matching TZ/QZ basis family.')
    return list(BASES[2:] if basis.startswith('aug-') else BASES[:2])


def digest(record):
    return hashlib.sha256(json.dumps(record, sort_keys=True,
        separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def seal(record):
    record.pop('record_sha256', None)
    record['record_sha256'] = digest(record)
    return record


def geometry_hash(symbols, positions):
    from kinbot.anl.dispatch import _geometry_hash
    return _geometry_hash(Atoms(symbols=symbols, positions=positions))


def validate(record, molecule):
    """Prove atom order/state and immutable geometry before staging or assembly."""
    payload = {k: v for k, v in record.items() if k != 'record_sha256'}
    if record.get('record_sha256') != digest(payload):
        raise ValueError('Accepted basis geometry record hash differs.')
    if (record.get('schema') != 1 or record.get('status') != 'normal'
            or record.get('stationary_point') != 'minimum'
            or record.get('symbols') != molecule['symbols']
            or record.get('charge') != molecule.get('charge', 0)
            or record.get('multiplicity') != molecule.get('multiplicity', 1)):
        raise ValueError('Accepted basis geometry atom order/state differs.')
    positions = np.asarray(record['positions'])
    if positions.shape != (len(molecule['symbols']), 3) or not np.isfinite(positions).all():
        raise ValueError('Invalid accepted basis geometry.')
    if record['geometry_sha256'] != geometry_hash(record['symbols'], record['positions']):
        raise ValueError('Accepted basis geometry hash differs.')
    if record['masses'] != Atoms(record['symbols']).get_masses().tolist():
        raise ValueError('Explicit isotopes require an isotope-aware ANL input; cannot reuse as default masses.')
    if not record.get('source_job') or record.get('source_row_id') is None:
        raise ValueError('Accepted basis geometry lacks a source calculation.')
    theory = record.get('theory', {})
    if (theory.get('calculator') != 'gaussian' or theory.get('method') != 'B2PLYP'
            or theory.get('dispersion') != 'GD3BJ' or canonical_basis(theory.get('basis')) is None):
        raise ValueError('Accepted basis geometry has an unsupported theory.')
    return record


def remember(species, qc, row):
    """Called only when the owning optimizer loads its selected calculation."""
    species.__dict__.pop('_anl_basis_point', None)
    species.__dict__.pop('anl_basis_reuse', None)
    seed = species.__dict__.pop('_anl_basis_seed', None)
    if getattr(species, 'wellorts', 1) != 0 or row.data.get('status') != 'normal':
        return
    symbols = list(map(str, species.atom))
    selected_hash = geometry_hash(symbols, species.geom)
    saved = row.data.get('anl_basis_reuse')
    if saved and saved.get('selected_geometry_sha256') == selected_hash:
        species.anl_basis_reuse = deepcopy(saved)
    route = getattr(qc, '_routes', {}).get(row.name, {})
    profile = getattr(qc, 'profiles', {}).get(route.get('level'))
    if (profile is None or route.get('source_profile_sha256', qc._fingerprint) != qc._fingerprint
            or profile.calculator != 'gaussian' or profile.method.upper() != 'B2PLYP'
            or canonical_basis(profile.basis) is None):
        return
    kwargs = {k.casefold(): str(v).upper() for k, v in profile.calculator_kwargs.items()}
    # Solvation, fields, custom basis/ECP and external guesses require their
    # own ANL contracts; numerical SCF/grid/symmetry choices do not change the
    # declared gas-phase method. An unsupported profile simply opts out.
    allowed = {'empiricaldispersion', 'scf', 'integral', 'symm', 'pop',
               'mem', 'nprocshared'}
    if kwargs.get('empiricaldispersion') != 'GD3BJ' or set(kwargs) - allowed:
        return
    if row.toatoms().get_chemical_symbols() != symbols:
        return
    values = [species.energy, species.zpe, *species.freq]
    if not np.isfinite(values).all():
        return
    point = seal(dict(schema=1, status='normal', stationary_point='minimum',
        symbols=symbols, positions=np.asarray(species.geom).tolist(),
        masses=row.toatoms().get_masses().tolist(), charge=int(species.charge),
        multiplicity=int(species.mult), geometry_sha256=selected_hash,
        theory=dict(calculator='gaussian', method='B2PLYP', basis=canonical_basis(profile.basis),
                    dispersion='GD3BJ'), source_job=row.name, source_row_id=row.id,
        profile_sha256=qc._fingerprint, database=str(Path(qc.db.filename).resolve()),
        energy_hartree=float(species.energy),
        harmonic_zpe_hartree=float(species.zpe), frequencies_cm_inverse=list(species.freq)))
    species._anl_basis_point = point
    if seed is not None and seed.get('refinement_job') == row.name:
        # Optimize.compare_structures has accepted the endpoint of this exact
        # seed's L2 optimization. No energy/RMSD guess across basis sets is used.
        species.anl_basis_reuse = dict(schema=1, selected_geometry_sha256=selected_hash,
            points=[seed, point], lineage='accepted L2 refinement of selected L1 minimum')


def seed(species, qc=None, job=None):
    """Bind reuse to the actual geometry handed to L2; HIR changes invalidate it."""
    species.__dict__.pop('_anl_basis_seed', None)
    point = getattr(species, '_anl_basis_point', None)
    if point is None and qc is not None and job and not getattr(species, 'wellorts', 1):
        # No conformer search: recover the accepted L1 seed directly, without
        # replacing the species' selected properties or projected frequencies.
        rows = list(qc.db.select(name=job, sort='-id', limit=1))
        if rows and np.array_equal(rows[0].positions, species.geom):
            row = rows[0]
            if all(row.data.get(k) is not None for k in ('energy', 'zpe', 'frequencies')):
                from kinbot import constants
                view = copy(species)
                view.energy = float(row.data['energy']) * constants.EVtoHARTREE
                view.zpe = float(row.data['zpe'])
                view.freq = list(row.data['frequencies'])
                remember(view, qc, row)
                point = getattr(view, '_anl_basis_point', None)
    if point and point['geometry_sha256'] == geometry_hash(list(map(str, species.atom)), species.geom):
        species._anl_basis_seed = deepcopy(point)


def bind_seed(species, qc, job):
    """Verify the rendered/reused L2 worker actually starts at this L1 point.

    Reading the worker also handles restart adoption: qc_opt may have found an
    existing job rather than rendering one from the current species geometry.
    Nothing is submitted, changed or invalidated if optional reuse cannot be
    established.
    """
    point = species.__dict__.pop('_anl_basis_seed', None)
    if point is None:
        return
    from kinbot.qc_reuse import definition
    backend = qc._backend_for(job) if hasattr(qc, '_backend_for') else qc
    worker = definition(backend, job, minimum_source=True)
    if (worker is None or worker['symbols'] != point['symbols']
            or worker['charge'] != point['charge']
            or worker['multiplicity'] != point['multiplicity']
            or not np.allclose(worker['positions'], point['positions'], atol=5.1e-9, rtol=0)):
        return
    point.update(refinement_job=job, refinement_worker_sha256=worker['worker_sha256'])
    species._anl_basis_seed = seal(point)


def exported(species):
    saved = deepcopy(getattr(species, 'anl_basis_reuse', None))
    if saved and saved['selected_geometry_sha256'] == geometry_hash(list(map(str, species.atom)), species.geom):
        return saved
    return None


def select(molecule, basis):
    source = molecule.get('source', {})
    saved = source.get('anl_basis_reuse')
    if not saved or saved.get('selected_geometry_sha256') != source.get('geometry_sha256'):
        return None
    for point in saved.get('points', []):
        validate(point, molecule)
        if point['theory'] == dict(calculator='gaussian', method='B2PLYP', basis=basis, dispersion='GD3BJ'):
            return deepcopy(point)
    return None


def remap(saved, permutation, symbols, positions):
    """Move all basis geometries with the same verified QC atom permutation."""
    result = deepcopy(saved)
    for point in result['points']:
        old = point['record_sha256']
        for key in ('symbols', 'positions', 'masses'):
            point[key] = [point[key][i] for i in permutation]
        point['geometry_sha256'] = geometry_hash(point['symbols'], point['positions'])
        point['mapping_source_sha256'] = old
        point['atom_permutation'] = list(permutation)
        seal(point)
    result['selected_geometry_sha256'] = geometry_hash(symbols, positions)
    return result
