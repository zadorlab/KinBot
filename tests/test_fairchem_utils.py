"""Model-loading options without requiring a FairChem installation."""

import sys
from types import ModuleType, SimpleNamespace
from unittest.mock import Mock

import pytest

from kinbot.fairchem_utils import load_predictor


@pytest.mark.parametrize('local', [False, True])
@pytest.mark.parametrize('enabled', [None, False, True])
@pytest.mark.parametrize('default_merge', [False, True])
def test_compile_settings_for_old_and_new_fairchem_defaults(
        tmp_path, monkeypatch, local, enabled, default_merge):
    modules = {}
    for name in ('fairchem', 'fairchem.core', 'fairchem.core.units',
                 'fairchem.core.units.mlip_unit', 'fairchem.core.units.mlip_unit.api',
                 'fairchem.core.units.mlip_unit.api.inference',
                 'fairchem.core.pretrained_mlip'):
        modules[name] = ModuleType(name)
        monkeypatch.setitem(sys.modules, name, modules[name])
        parent, _, child = name.rpartition('.')
        if parent:
            setattr(modules[parent], child, modules[name])
    defaults = dict(compile=True, merge_mole=default_merge, activation_checkpointing=False)
    modules['fairchem.core.units.mlip_unit.api.inference'].inference_settings_default = (
        lambda: SimpleNamespace(**defaults))
    local_loader = Mock()
    named_loader = Mock()
    modules['fairchem.core.units.mlip_unit'].load_predict_unit = local_loader
    modules['fairchem.core.pretrained_mlip'].get_predict_unit = named_loader
    checkpoint = tmp_path / 'model.pt'
    checkpoint.touch()
    model = str(checkpoint) if local else 'uma-s-1p1'
    options = {} if enabled is None else {'compile_model': enabled}
    result = load_predictor(model, 'cuda:1', **options)
    chosen, unused = (local_loader, named_loader) if local else (named_loader, local_loader)
    chosen.assert_called_once()
    unused.assert_not_called()
    assert result is chosen.return_value
    assert chosen.call_args.args == (model,)
    assert chosen.call_args.kwargs['device'] == 'cuda:1'
    settings = chosen.call_args.kwargs['inference_settings']
    assert vars(settings) == dict(defaults, compile=bool(enabled),
                                 merge_mole=True if enabled else default_merge)
    assert defaults['compile'] is True
