"""Shared FairChem model loading."""

import os


def load_predictor(model, device):
    """Load a local checkpoint or a named pretrained model."""
    if os.path.isfile(model):
        from fairchem.core.units.mlip_unit import load_predict_unit
        return load_predict_unit(model, device=device)
    from fairchem.core import pretrained_mlip
    return pretrained_mlip.get_predict_unit(model, device=device)
