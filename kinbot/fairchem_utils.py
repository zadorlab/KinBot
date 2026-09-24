"""Shared FairChem model loading."""

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


def load_predictor(model, device):
    """Load a local checkpoint or a named pretrained model."""
    if os.path.isfile(model):
        from fairchem.core.units.mlip_unit import load_predict_unit
        return load_predict_unit(model, device=device)
    from fairchem.core import pretrained_mlip
    try:
        return pretrained_mlip.get_predict_unit(model, device=device)
    except Exception as error:
        # huggingface_hub is optional and imported transitively by FairChem.
        # Match the public exception name without importing that optional
        # package in legacy KinBot installations.
        if type(error).__name__ == 'GatedRepoError':
            raise _gated_model_error(model) from error
        raise
