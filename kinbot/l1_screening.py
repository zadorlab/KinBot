"""Optional absolute L1 energy screen, with read-only counterfactual reporting.

An energetic screen is an explicit search restriction, not a rate prediction.
Missing or incomparable evidence retains a channel. Accepted channels are not
removed by this module; changing an existing network requires a selection.
"""
from collections import Counter
from dataclasses import asdict
import argparse
import hashlib
import json
import math
import re
from pathlib import Path
import sqlite3

from kinbot import constants
from kinbot.execution import atomic_json


def validate(limit):
    if limit is not None and (isinstance(limit, bool) or not isinstance(limit, (int, float))
                              or not math.isfinite(limit) or limit <= 0):
        raise ValueError('l1_pes_energy_ceiling must be null or positive kcal/mol.')


def theory(parameters):
    from kinbot.theory import resolve_profiles
    profiles = resolve_profiles(parameters)
    if 'l1' not in profiles:
        raise ValueError('Absolute L1 screening needs recorded per-level theory profiles.')
    return normalized_profile(asdict(profiles['l1']))


def normalized_profile(profile):
    profile = dict(profile)
    # Location/labels do not change the energy model; retain all other options.
    for key in ('command', 'label', 'device', 'name'):
        profile.pop(key, None)
    return profile


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def finite_number(value):
    return (not isinstance(value, bool) and isinstance(value, (int, float))
            and math.isfinite(value))


def comparable_family(parameters, profile, family):
    """Do legacy family overrides actually change the recorded L1 model?

    Profiled Gaussian inputs default scan/BLS theory to L1, whereas explicit
    overrides can still select a different model. Neural backends do not use
    these method/basis overrides. Include search-checkpoint family spellings.
    """
    if family == 'unknown':
        return False
    special = family in (*constants.mp2_list, 'radditionmultiplebond',
                          'biradrecombinationr', 'r12cycloaddition', 'r14iradscission')
    bls = family in ('barrierless_saddle', 'barrierlesssaddle')
    if not (special or bls) or profile.get('calculator') in ('fairchem', 'nn_pes'):
        return True
    prefix = 'scan' if special else 'barrierless_saddle'
    return all(str(parameters.get(prefix+'_'+key, '')).casefold()
               == str(profile.get(key, '')).casefold() for key in ('method', 'basis'))


def read_points(directory, names, profile=None, *, geometry=False, require_zpe=True):
    """Read latest rows without creating a database, journal, or ASE lock file.

    Geometry-only IRC endpoint identification must not require a Hessian. An
    absent ZPE stays absent, including when electronic energies are available.
    """
    from ase.db import connect
    directory = Path(directory)
    if profile is not None:
        routes_path = directory/'.kinbot_theory_jobs.json'
        if not routes_path.exists():
            raise ValueError('No recorded L1 theory routes; evidence remains unclassified.')
        routes = json.loads(routes_path.read_text())
        if normalized_profile(routes.get('profiles', {}).get('l1', {})) != profile:
            raise ValueError('Recorded L1 profile differs from the requested profile.')
        for name in names:
            entry = routes.get('jobs', {}).get(name, {})
            # IRC templates optimize the product in the same submitted L1
            # worker; the child row has no separate scheduler route.
            if not entry and re.search(r'_IRC_[FR]_prod$', name):
                entry = routes.get('jobs', {}).get(name[:-5], {})
            if (entry.get('level') != 'l1' or entry.get('source_profile_sha256', routes.get('profile_sha256'))
                    != routes.get('profile_sha256')):
                raise ValueError(f'{name}: L1 route/generation is unverified.')
    path = (directory/'kinbot.db').resolve(strict=True)
    db = connect(str(path), use_lock_file=False)
    db.connection = sqlite3.connect(path.as_uri()+'?mode=ro', uri=True)
    db.change_count = 0
    points = []
    try:
        for name in names:
            rows = list(db.select(name=name, sort='-id', limit=1))
            if not rows or rows[0].data.get('status') != 'normal':
                raise ValueError(f'{name}: no normally completed L1 row.')
            row = rows[0]
            energy, zpe = row.data.get('energy'), row.data.get('zpe')
            if not finite_number(energy):
                raise ValueError(f'{name}: missing/nonfinite electronic energy.')
            if require_zpe and not finite_number(zpe):
                raise ValueError(f'{name}: missing/nonfinite ZPE.')
            points.append(dict(name=name, row_id=row.id, symbols=row.toatoms().get_chemical_symbols(),
                energy_hartree=energy*constants.EVtoHARTREE,
                zpe_hartree=zpe if finite_number(zpe) else None))
            if not finite_number(zpe):
                points[-1]['zpe_status'] = 'missing' if zpe is None else 'invalid'
            if geometry:
                points[-1]['positions'] = row.positions.tolist()
                points[-1]['frequencies'] = list(map(float, row.data.get('frequencies', [])))
    finally:
        db.connection.close()
    return points


def initial(case):
    case = Path(case)
    phase = case/'pes_exploration.json'
    if phase.exists():
        name = json.loads(phase.read_text()).get('initial_worker')
        if name:
            return name
    names = (case/'chemids').read_text().splitlines()
    if not names:
        raise ValueError('Initial PES well is not recorded.')
    return names[0]


def reference(case):
    from kinbot.parameters import Parameters
    name = initial(case)
    worker = Path(case)/name
    par = Parameters(str(worker/(name+'.json'))).par
    point = read_points(worker, [name+'_well'], theory(par))[0]
    return dict(schema=1, worker=name, database=str(worker/'kinbot.db'), point=point,
                charge=par['charge'], theory=theory(par))


def evaluate(ref, profile, points, limit, *, comparable=True, charge=None):
    validate(limit)
    result = dict(decision='retain', limit_kcal_mol=limit, points=points,
                  reference_sha256=digest(ref))
    if not comparable or ref['theory'] != profile or charge != ref['charge']:
        return dict(result, decision='unknown', reason='Different or unverified L1 theory/charge.')
    symbols = Counter(s for p in points for s in p['symbols'])
    if not points or symbols != Counter(ref['point']['symbols']):
        return dict(result, decision='unknown', reason='Different stoichiometry; cannot compare absolute energies.')
    if any(not finite_number(p.get(key)) for p in [*points, ref['point']]
           for key in ('energy_hartree', 'zpe_hartree')):
        return dict(result, decision='unknown', reason='Missing/nonfinite energy or ZPE; no zero substitution.')
    value = sum(p['energy_hartree']+p['zpe_hartree'] for p in points)
    zero = ref['point']['energy_hartree']+ref['point']['zpe_hartree']
    relative = (value-zero)*constants.AUtoKCAL
    if not math.isfinite(relative):
        return dict(result, decision='unknown', reason='Nonfinite energy.')
    return dict(result, relative_kcal_mol=relative,
                decision='exclude' if limit is not None and relative > limit+1.e-7 else 'retain')


def screen(parameters, names, *, reaction, kind, comparable=True, directory=None):
    """Only call at a completed L1 stage, before any corresponding L2 submit."""
    limit = parameters.get('l1_pes_energy_ceiling')
    if limit is None:
        return False
    validate(limit)
    directory = Path(directory or Path.cwd()).resolve()
    # Recursive PES workers share a single immutable initial-well reference.
    case = directory.parent
    try:
        ref = reference(case)
        points = read_points(directory, names, theory(parameters))
        result = evaluate(ref, theory(parameters), points, limit,
                          comparable=comparable, charge=parameters['charge'])
    except (OSError, ValueError, sqlite3.Error) as error:
        result = dict(decision='unknown', reason=str(error), limit_kcal_mol=limit)
    path = directory/'l1_screening.json'
    record = json.loads(path.read_text()) if path.exists() else dict(schema=1, decisions={})
    key = reaction+':'+kind
    item = dict(result, reaction=reaction, kind=kind, jobs=names)
    previous = record['decisions'].get(key)
    if previous != item:
        if previous is not None:
            record.setdefault('history', []).append(previous)
        record['decisions'][key] = item
        atomic_json(path, record)
    return result['decision'] == 'exclude'


def report(case, limit):
    from kinbot.l1_audit import report as audit
    return audit(case, limit)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case')
    parser.add_argument('--limit', type=float, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(report(args.case, args.limit), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
