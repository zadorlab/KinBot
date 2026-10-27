# PESViewer offline installation artifact

This optional, platform-independent wheel contains unmodified PESViewer source
from upstream commit `ee0e18fcdb6bc8de6fddaa3e4726a0143c785a73`:

https://github.com/zadorlab/PESViewer/commit/ee0e18fcdb6bc8de6fddaa3e4726a0143c785a73

That revision replaces RDKit's removed drawing interfaces. Upstream still calls
the package 1.1.0, so the version alone cannot distinguish it from the older
conda package. `provenance.json` records the wheel checksum and hashes of every
Python source file, verified byte for byte against the pinned upstream source.
The BSD license is included both here and in the wheel.

The wheel was built with `pip` from the upstream source archive. It contains no
compiled code, dependency wheels or local source modifications. It is not
installed automatically with KinBot. Install it only in the separate viewer
environment, which must already contain PESViewer's runtime dependencies:

```bash
python -m pip --disable-pip-version-check install --no-index --no-deps \
  --force-reinstall external/pesviewer_wheel/pesviewer-1.1.0-py3-none-any.whl
```

See `docs/pesviewer_export.md` for headless plotting. The local wheel avoids a
source build and its setuptools download; it does not alter TLS trust settings.
