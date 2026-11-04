"""Evidence for sharing new downstream QC of independently optimized minima.

This does not merge wells, optical populations, conformer ensembles or reactions.
Failure to establish identity merely leaves the ordinary calculation in place.
"""
import json
from pathlib import Path

import numpy as np

from kinbot import constants
from kinbot.execution import atomic_json, file_sha256
from kinbot.spatial import proper_alignment

# Same resolution as accepted stationary-point reuse, not a search cutoff.
LIMITS = dict(max_distance_angstrom=5.e-4, energy_hartree=1.e-6,
              zpe_hartree=1.e-5, frequency_cm_inverse=1.)


def capture(qc, species):
    """Freeze a selected, normally completed native minimum and its provenance.

    Native modes establish that the spectrum belongs to this geometry. A
    published alias is followed to its original script; mutable current input
    settings are never substituted for the calculation that actually ran.
    """
    if (getattr(species, 'wellorts', 0) or not getattr(species, 'source_job', None)
            or any(getattr(species, 'isotopes', []))):
        return None
    from kinbot.qc_reuse import definition, native_record
    job = species.source_job
    try:
        rows = list(qc.db.select(name=job))
        if not rows:
            return None
        row = rows[-1]
        if (row.id != species.source_row_id or row.data.get('status') != 'normal'
                or list(row.symbols) != list(map(str, species.atom))
                or not np.array_equal(row.positions, species.geom)):
            return None
        data = dict(row.data)
        freq = np.asarray(data['frequencies'], float)
        if (not np.isfinite([data['energy'], data['zpe']]).all()
                or not np.isfinite(freq).all() or np.any(freq <= 0)
                or (len(row.symbols) > 1 and len(freq) not in
                    (3*len(row.symbols)-6, 3*len(row.symbols)-5))):
            return None
        # A copied conventional name is not a separate calculation generation.
        original = row
        seen = set()
        while original.data.get('copied_from_job'):
            key = (original.data['copied_from_job'], original.data['copied_from_row_id'])
            if key in seen:
                return None
            seen.add(key)
            original = qc.db.get(id=key[1])
            if original.name != key[0] or not np.array_equal(original.positions, row.positions):
                return None
        shared = native_record(row)
        if shared:
            # A shared result's proof was captured before publication; native
            # files keep their original frame and ordered atoms.
            bundle, payload, permutation = shared
            script, log = bundle/'native.py', bundle/'native.log'
            symbols = payload['symbols']
            coordinates = payload['positions']
            source_job = payload['job']
        else:
            script, log = Path(original.name+'.py'), Path(original.name+'.log')
            symbols, coordinates = list(row.symbols), row.positions
            source_job = original.name
            if not Path(original.name+'_sella.log').read_text().rstrip().endswith('done'):
                return None
            ident = getattr(qc, 'job_ids', {}).get(original.name)
            if ident and Path(original.name+'.exitcode').read_text().split() != [str(ident), '0']:
                return None
        spec = definition(qc, source_job, path=script, minimum_source=True)
        if (spec is None or spec['kind'] != 'minimum'
                or spec['charge'] != species.charge or spec['multiplicity'] != species.mult):
            return None
        # Do not reuse a nondefault isotope calculation through a template
        # which has no isotope/mass contract.
        from ase import Atoms
        if not np.array_equal(row.toatoms().get_masses(), Atoms(row.symbols).get_masses()):
            return None
        before = {str(p): file_sha256(p) for p in (script, log)}
        from kinbot.gaussian_modes import read_normal_modes
        if len(symbols) > 1:
            native = read_normal_modes(log, symbols, coordinates)
            if not np.allclose(np.sort(native['frequencies_cm1']), np.sort(freq),
                               atol=.001, rtol=0):
                return None
        elif 'Normal termination of Gaussian' not in log.read_text()[-6000:]:
            return None
        if before != {str(p): file_sha256(p) for p in (script, log)}:
            return None
        return dict(schema=1, source_job=job, source_row_id=row.id,
            database=str(Path(qc.db.filename).resolve()), native_hashes=before,
            theory_program=spec['program_sha256'], symbols=list(row.symbols),
            masses=row.toatoms().get_masses().tolist(), positions=row.positions.tolist(),
            energy_hartree=float(data['energy'])*constants.EVtoHARTREE,
            zpe_hartree=float(data['zpe']), frequencies_cm1=freq.tolist())
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        # Missing evidence is not a failed minimum and never triggers a retry.
        return None


def attach(qc, species, job):
    """Bind optional minimum evidence to this exact newly rendered request."""
    from kinbot.qc_reuse import SUFFIXES, _marker
    if (getattr(qc, 'qc', None) != 'gauss' or not getattr(qc, 'use_sella', False)
            or job in getattr(qc, 'job_ids', {}) or qc.is_in_database(job)
            or _marker(job, 'consumer').exists() or _marker(job, 'producer').exists()
            or any(Path(job+s).exists() for s in (*SUFFIXES, '.pkl', '.sbatch', '.exitcode'))):
        return
    evidence = capture(qc, species)
    if evidence is not None and job.endswith('_well_high'):
        from kinbot.qc_reuse import definition
        request = definition(qc, job)
        if request is None or not np.array_equal(request['positions'], evidence['positions']):
            evidence = None
    if evidence is not None:
        atomic_json(Path(job+'.minimum_reuse.json'), dict(
            worker_sha256=file_sha256(job+'.py'), minimum=evidence))


def comparison(a, b, permutation):
    """Check two completed minima under ONE proposed request atom mapping.

    permutation indexes b in a's atom order. A successful comparison permits
    reuse within QC resolution; it is not permission to coalesce PES wells.
    """
    if not a or not b or a['theory_program'] != b['theory_program']:
        return None
    p = np.asarray(permutation, int)
    if (a['symbols'] != [b['symbols'][i] for i in p]
            or not np.array_equal(a['masses'], np.asarray(b['masses'])[p])):
        return None
    af, bf = np.sort(a['frequencies_cm1']), np.sort(b['frequencies_cm1'])
    if af.shape != bf.shape:
        return None
    measurements = dict(
        max_distance_angstrom=float(max(proper_alignment(a['positions'],
                                    np.asarray(b['positions'])[p]).errors)),
        energy_hartree=abs(a['energy_hartree']-b['energy_hartree']),
        zpe_hartree=abs(a['zpe_hartree']-b['zpe_hartree']),
        frequency_cm_inverse=float(np.max(np.abs(af-bf))) if len(af) else 0.)
    if any(not np.isfinite(v) or v > LIMITS[k] for k, v in measurements.items()):
        return None
    return dict(measurements=measurements, limits=LIMITS,
                a_source=[a['database'], a['source_job'], a['source_row_id']],
                b_source=[b['database'], b['source_job'], b['source_row_id']])


def canonical_rotors(species, qc):
    """Choose common terminal atoms for NEW scans of verified identical minima.

    Atom-order-dependent terminal choices describe different constrained
    optimizations even on the same axis. Agree on the definition before either
    claiming equivalence or submitting a scan. Keep local rotor numbering.
    Already staged/running/completed scans are never redefined by this cache.
    """
    from kinbot.qc_reuse import _scope, _locked, _read, matching, digest, NAME
    from kinbot.species_routing import routing_name
    import shutil
    import uuid
    root = _scope()
    if root is None or species.wellorts or not species.dihed:
        return
    name = routing_name(species)
    prefix = name+'_hir_'
    if (any(Path('hir').glob(prefix+'*'))
            or any(j.startswith('hir/'+prefix) for j in getattr(qc, 'job_ids', {}))
            or any(True for _ in qc.db.select(name='hir/'+prefix+'*'))):
        return
    evidence = capture(qc, species)
    if evidence is None:
        return
    wanted = dict(contract_sha256=evidence['theory_program'], symbols=evidence['symbols'],
                  positions=evidence['positions'], minimum=evidence, kind='minimum', constraint_atoms=[])
    key = digest([evidence['theory_program'], sorted(evidence['symbols'])])
    path = root/NAME/'minima'/(key+'.json')
    original = [list(map(int, rotor)) for rotor in species.dihed]
    with _locked(path.with_suffix('.lock')):
        records = _read(path) if path.exists() else []
        for record in records:
            match = matching(wanted, record['definition'])
            if match is None or match['minimum_comparison'] is None:
                continue
            p = match['permutation']
            inverse = np.argsort(p)
            canonical = [[int(inverse[i]) for i in rotor] for rotor in record['rotors']]
            by_axis = {tuple(sorted(r[1:3])): r for r in canonical}
            # A different selected axis set is a different rotor model.
            if (len(by_axis) != len(original)
                    or set(by_axis) != {tuple(sorted(r[1:3])) for r in original}):
                continue
            mapped = []
            for rotor in original:
                other = by_axis[tuple(sorted(rotor[1:3]))]
                mapped.append(other if other[1:3] == rotor[1:3] else other[::-1])
            species.dihed = mapped
            atomic_json(Path('hir')/(name+'.minimum_mapping.json'), dict(schema=1,
                source=record['definition']['minimum'], consumer=evidence, mapping=match,
                original_dihedrals=original, mapped_dihedrals=mapped,
                archived_native=record['archived_native']))
            return
        # Retain original native evidence for the canonical reference even if
        # its owning driver later explicitly revises that calculation.
        directory = root/NAME/'minima'/'evidence'/uuid.uuid4().hex
        directory.mkdir(parents=True)
        archived = {}
        for i, (filename, sha) in enumerate(evidence['native_hashes'].items()):
            target = directory/(str(i)+Path(filename).suffix)
            shutil.copyfile(filename, target)
            if file_sha256(target) != sha or file_sha256(filename) != sha:
                raise RuntimeError('Selected minimum changed while recording its rotor reference.')
            archived[str(target.relative_to(root))] = sha
        records.append(dict(definition=wanted, rotors=original, archived_native=archived))
        atomic_json(path, records)
