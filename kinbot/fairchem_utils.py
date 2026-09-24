"""Shared FairChem model loading."""

from contextlib import contextmanager
import os


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


@contextmanager
def _offline_hub_environment():
    """Hide invalid generic SOCKS proxies during an offline cache lookup."""
    removed = {}
    if _enabled('HF_HUB_OFFLINE'):
        for name in ('ALL_PROXY', 'all_proxy'):
            value = os.environ.get(name, '')
            if value.lower().startswith('socks://'):
                removed[name] = os.environ.pop(name)
    try:
        yield
    finally:
        os.environ.update(removed)


def load_predictor(model, device):
    """Load a local checkpoint or a named pretrained model."""
    if os.path.isfile(model):
        from fairchem.core.units.mlip_unit import load_predict_unit
        return load_predict_unit(model, device=device)
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
