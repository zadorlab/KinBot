"""Load FairChem models with KinBot's compilation setting."""

import os


def load_predictor(model, device, *, compile_model=False):
    """Load a local checkpoint or a named pretrained model."""
    from fairchem.core.units.mlip_unit.api.inference import inference_settings_default

    # FairChem 2.22 enabled compilation by default. Short KinBot jobs often
    # finish before faster force evaluations repay that startup cost.
    settings = inference_settings_default()
    settings.compile = compile_model
    if compile_model:
        # FairChem 2.14 can fail in compiled UMA force evaluations without
        # precomputed composition weights. A KinBot job keeps its composition,
        # charge and spin fixed, as required by this preparation.
        settings.merge_mole = True
    if os.path.isfile(model):
        from fairchem.core.units.mlip_unit import load_predict_unit
        return load_predict_unit(model, device=device, inference_settings=settings)
    from fairchem.core import pretrained_mlip
    return pretrained_mlip.get_predict_unit(model, device=device,
                                          inference_settings=settings)
