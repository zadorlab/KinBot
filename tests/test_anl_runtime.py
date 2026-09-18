"""Native QC runtime checks are scoped to the launched program."""

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest

from kinbot.anl.runtime import qc_runtime_environment


def test_cfour_finds_only_missing_fortran_library_for_child():
    with TemporaryDirectory() as temporary:
        library = Path(temporary) / 'libgfortran.so.4'
        library.write_bytes(b'fixture')

        def missing(_, env):
            return [] if env.get('LD_PRELOAD') == str(library) else ['libgfortran.so.4']

        original = {'PATH': '/site/bin', 'LD_LIBRARY_PATH': '/site/current/lib'}
        with patch('kinbot.anl.runtime.shutil.which', return_value='/site/bin/xcfour'), \
                patch('kinbot.anl.runtime._missing_libraries', side_effect=missing), \
                patch('kinbot.anl.runtime._runtime_candidates', return_value=[library]):
            child, provenance = qc_runtime_environment('xcfour', 'cfour', original)
        assert original == {'PATH': '/site/bin', 'LD_LIBRARY_PATH': '/site/current/lib'}
        assert child['LD_LIBRARY_PATH'] == '/site/current/lib'
        assert child['LD_PRELOAD'] == str(library)
        assert provenance == {'runtime_library': str(library.resolve())}


def test_cfour_uses_module_runtime_when_available():
    with patch('kinbot.anl.runtime.shutil.which', return_value='/site/bin/xcfour'), \
            patch('kinbot.anl.runtime._missing_libraries', return_value=[]):
        child, provenance = qc_runtime_environment('xcfour', 'cfour', {'PATH': '/site/bin'})
    assert child == {'PATH': '/site/bin'}
    assert provenance == {}


def test_cfour_rejects_unresolved_runtime_before_launch():
    with TemporaryDirectory() as temporary:
        library = Path(temporary) / 'libgfortran.so.4'
        library.write_bytes(b'fixture')
        with patch('kinbot.anl.runtime.shutil.which', return_value='/site/bin/xcfour'), \
                patch('kinbot.anl.runtime._missing_libraries',
                      return_value=['libgfortran.so.4']), \
                patch('kinbot.anl.runtime._runtime_candidates', return_value=[library]):
            with pytest.raises(RuntimeError, match='CFOUR needs libgfortran.so.4'):
                qc_runtime_environment('xcfour', 'cfour', {'PATH': '/site/bin'})
