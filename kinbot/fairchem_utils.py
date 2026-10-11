"""Shared FairChem model loading."""

from contextlib import contextmanager
import os
from pathlib import Path
import re

from kinbot.network import without_invalid_generic_proxy


def _gated_model_error(model):
    return RuntimeError(
        f"FairChem model {model!r} is gated and the active Hugging Face "
        "account has not been granted access. While logged into that same "
        "account, request or accept access at "
        "https://huggingface.co/facebook/UMA, create a read token with "
        "permission for public gated repositories, and run "
        "`hf auth login --force` before retrying."
    )


def _enabled(name):
    return os.environ.get(name, '').strip().upper() in {
        '1', 'ON', 'YES', 'TRUE',
    }


def _huggingface_hub_roots():
    """Return standard Hugging Face Hub cache roots without network access."""
    configured = [os.environ.get('HF_HUB_CACHE'),
                  os.environ.get('HUGGINGFACE_HUB_CACHE')]
    hf_home = os.environ.get('HF_HOME')
    if hf_home:
        configured.append(str(Path(hf_home).expanduser() / 'hub'))
    else:
        cache_home = os.environ.get('XDG_CACHE_HOME')
        base = (Path(cache_home).expanduser() if cache_home
                else Path.home() / '.cache')
        configured.append(str(base / 'huggingface' / 'hub'))
    roots = []
    for value in configured:
        if value:
            root = Path(value).expanduser()
            if root not in roots:
                roots.append(root)
    return roots


def _cached_uma_checkpoint(model):
    """Resolve a downloaded UMA snapshot even when ``refs/main`` is absent."""
    if (not isinstance(model, str)
            or not re.fullmatch(r'uma-[A-Za-z0-9._-]+', model,
                                flags=re.IGNORECASE)):
        return None
    for hub in _huggingface_hub_roots():
        repository = hub / 'models--facebook--UMA'
        revision = repository / 'refs' / 'main'
        if revision.is_file():
            commit = revision.read_text(errors='replace').strip()
            candidate = (repository / 'snapshots' / commit / 'checkpoints'
                         / f'{model}.pt')
            if candidate.is_file():
                return candidate.resolve()
        candidates = [path for path in repository.glob(
            f'snapshots/*/checkpoints/{model}.pt') if path.is_file()]
        if candidates:
            return max(candidates, key=lambda path: path.stat().st_mtime).resolve()
    return None


@contextmanager
def _offline_hub_environment():
    """Hide invalid generic SOCKS proxies during an offline cache lookup."""
    if _enabled('HF_HUB_OFFLINE'):
        with without_invalid_generic_proxy():
            yield
    else:
        yield


def load_predictor(model, device):
    """Load a local checkpoint or a named pretrained model."""
    if os.path.isfile(model):
        from fairchem.core.units.mlip_unit import load_predict_unit
        return load_predict_unit(model, device=device)
    cached = _cached_uma_checkpoint(model)
    if cached is not None:
        from fairchem.core.units.mlip_unit import load_predict_unit
        return load_predict_unit(str(cached), device=device)
    from fairchem.core import pretrained_mlip
    try:
        # Some institutional environments export the nonstandard
        # ``socks://`` scheme. HTTPX constructs its proxy map before the Hub's
        # offline transport rejects network access, so even a cache hit fails.
        # HTTP(S)_PROXY is left untouched, and the original process environment
        # is restored as soon as model loading finishes.
        with _offline_hub_environment():
            return pretrained_mlip.get_predict_unit(model, device=device)
    except Exception as error:
        # huggingface_hub is optional and imported transitively by FairChem.
        # Match the public exception name without importing that optional
        # package in legacy KinBot installations.
        if type(error).__name__ == 'GatedRepoError':
            raise _gated_model_error(model) from error
        raise
