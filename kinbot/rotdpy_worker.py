"""Execute a ROTD sample and atomically publish its completed pickle.

The generated input pickle stays RUNNING throughout execution and batch I/O
retries. The external controller can never ingest a half-written completion.
"""
import os
from pathlib import Path
import pickle
import random
import sys
import tempfile
import time

import numpy as np


def run(name):
    from rotd_py.system import FluxTag
    seed = os.environ.get('KINBOT_SLURM_IO_RETRY_SEED')
    if seed is not None:
        np.random.seed(int(seed))
        random.seed(int(seed))
    path = Path(name + '.pkl')
    for attempt in range(4):
        try:
            with path.open('rb') as stream:
                int_tag, flux, surf_id, face_id, length, sample_id, status = pickle.load(stream)
            break
        except (FileNotFoundError, EOFError, pickle.UnpicklingError):
            if attempt == 3:
                raise
            time.sleep(5)
    tag = FluxTag(int_tag)
    if tag == FluxTag.FLUX_TAG:
        flux.run(length, face_id, sample_id)
    elif tag == FluxTag.SURF_TAG:
        flux.run_surf(length)
    elif tag != FluxTag.STOP_TAG:
        raise ValueError(f'Unknown ROTD worker tag: {int_tag}.')
    descriptor, filename = tempfile.mkstemp(prefix=path.name + '.', suffix='.tmp', dir=path.parent)
    temporary = Path(filename)
    try:
        with os.fdopen(descriptor, 'wb') as stream:
            pickle.dump([tag.value, flux, surf_id, face_id, length, sample_id, 'COMPLETED'], stream)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == '__main__':
    run(sys.argv[1])
