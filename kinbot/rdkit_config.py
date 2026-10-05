"""The RDKit version and stereo-perception settings used by KinBot."""
import re

from rdkit import Chem, rdBase


MIN_RDKIT_VERSION = '2025.9.3'


def configure_rdkit():
    """Require the tested release floor and select the perception algorithms."""
    version = rdBase.rdkitVersion
    match = re.fullmatch(r'(\d+)\.(\d+)\.(\d+)', version)
    minimum = tuple(map(int, MIN_RDKIT_VERSION.split('.')))
    if match is None or tuple(map(int, match.groups())) < minimum:
        raise ImportError(f'KinBot requires RDKit >= {MIN_RDKIT_VERSION}; found {version}.')
    # These are RDKit algorithm settings, not optional KinBot fallbacks.
    # Retain the modes used to validate the supported molecular identities.
    Chem.SetUseLegacyStereoPerception(True)
    Chem.SetAllowNontetrahedralChirality(True)


def rdkit_runtime():
    configure_rdkit()
    return {'rdkit_version': rdBase.rdkitVersion,
            'use_legacy_stereo_perception': Chem.GetUseLegacyStereoPerception(),
            'allow_nontetrahedral_chirality': Chem.GetAllowNontetrahedralChirality()}


def log_rdkit(logger):
    runtime = rdkit_runtime()
    logger.info('RDKit %s; UseLegacyStereoPerception=%s; AllowNontetrahedralChirality=%s',
                runtime['rdkit_version'], runtime['use_legacy_stereo_perception'],
                runtime['allow_nontetrahedral_chirality'])
