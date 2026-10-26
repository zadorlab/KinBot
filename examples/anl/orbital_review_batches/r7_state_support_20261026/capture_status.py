"""Use the shared read-only failure/request collector from the preceding panel."""
import runpy
from pathlib import Path

if __name__ == '__main__':
    runpy.run_path(str(Path(__file__).resolve().parent.parent/
                      'r7_root_spectrum_20261026/capture_status.py'), run_name='__main__')
