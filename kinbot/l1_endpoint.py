"""Owner-driven, fixed-geometry L1 endpoint ZPE for matched vdW screening."""
import json
from pathlib import Path
import sqlite3

import numpy as np

from kinbot import l1_screening
from kinbot.execution import atomic_json
from kinbot.frequencies import frequency_validation_error


def ensure(qc, species, parameters, *, family, directory=None):
    """Return complete/waiting/unresolved; never invent a missing endpoint ZPE.

    Each source database observation has one deterministic frequency job. The
    existing QC owner reconciles its allocation and any configured recovery.
    A terminal scientific failure remains unresolved rather than spawning a
    new attempt on every polling cycle or after an owner restart.
    """
    directory = Path(directory or Path.cwd()).resolve()
    name = species.name
    result = dict(status='unresolved', source_job=name)
    try:
        profile = l1_screening.theory(parameters)
        if not l1_screening.comparable_family(parameters, profile, family):
            raise ValueError('Endpoint L1 theory is not comparable to the configured fragment model.')
        source = l1_screening.read_points(directory, [name], profile, geometry=True,
                                          require_zpe=False)[0]
        result.update(source=source, source_sha256=l1_screening.digest(source), theory=profile)
        if (source['symbols'] != list(species.atom)
                or not np.allclose(source['positions'], species.geom, atol=1.e-6, rtol=0.)):
            raise ValueError('Selected IRC endpoint differs from its source observation.')
        threshold = parameters['imagfreq_threshold']
        valid = (l1_screening.finite_number(source['zpe_hartree']) and
                 (len(source['symbols']) == 1 and source['zpe_hartree'] == 0. or
                  not frequency_validation_error(source['frequencies'], 0, threshold)))
        if valid:
            result.update(status='complete', point=source)
        else:
            job = qc.qc_freq(species, name, high_level=0, endpoint=True)
            result['frequency_job'] = job
            status = qc.check_qc(job)
            if status == 'error':
                raise ValueError('Fixed-geometry endpoint frequency job failed: '+job)
            if status != 'normal':
                result.update(status='waiting', reason='Waiting for endpoint frequency completion.')
            else:
                point = l1_screening.read_points(directory, [job], profile, geometry=True)[0]
                if (point['symbols'] != source['symbols'] or not np.allclose(
                        point['positions'], source['positions'], atol=1.e-6, rtol=0.)):
                    raise ValueError('Endpoint frequency job changed the selected geometry.')
                if abs(point['energy_hartree']-source['energy_hartree']) > 2.e-6:
                    raise ValueError('Endpoint frequency energy differs from the source by more than 2 microhartree; inspect electronic state.')
                error = frequency_validation_error(point['frequencies'], 0, threshold)
                if error:
                    raise ValueError(error)
                result.update(status='complete', point=point)
    except (OSError, ValueError, sqlite3.Error) as error:
        result.update(status='unresolved', reason=str(error))
    path = directory/'l1_endpoint_checks.json'
    record = json.loads(path.read_text()) if path.exists() else dict(schema=1, checks={})
    old = record['checks'].get(name)
    if old != result:
        if old is not None:
            record.setdefault('history', []).append(old)
        record['checks'][name] = result
        atomic_json(path, record)
    return result
