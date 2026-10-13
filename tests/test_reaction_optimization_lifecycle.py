import json
from types import SimpleNamespace
from unittest.mock import Mock

from kinbot.reaction_generator import ReactionGenerator


def optimizer(*, shigh=0, shir=0, just_high=False, after=None, name='point'):
    point = SimpleNamespace(hir=None)
    result = SimpleNamespace(
        name=name,
        species=point,
        just_high=just_high,
        shigh=shigh,
        shir=shir,
        sconf=1,
        restart=0,
    )

    def poll():
        if after is not None:
            result.shigh, result.shir = after

    result.do_optimization = Mock(side_effect=poll)
    return result


def generator(reaction_type='test'):
    parent = SimpleNamespace(reac_type=[reaction_type])
    return ReactionGenerator(parent, {}, None, None)


def test_stage_four_uses_post_poll_completion_state():
    ts = optimizer(after=(1, 1), name='ts')
    product = optimizer(after=(1, 1), name='product')
    reaction = SimpleNamespace(
        ts_opt=ts,
        prod_opt=[product],
        do_vdW=False,
    )

    failures, pending, snapshots = generator()._poll_stage_four(0, reaction)

    assert failures == []
    assert pending == []
    assert [item['label'] for item in snapshots] == [
        'transition_state', 'product_0']
    ts.do_optimization.assert_called_once_with()
    product.do_optimization.assert_called_once_with()


def test_stage_four_deduplicates_shared_product_optimizer():
    product = optimizer(after=(1, 1))
    reaction = SimpleNamespace(prod_opt=[product, product], do_vdW=False)

    failures, pending, snapshots = generator('hom_sci')._poll_stage_four(
        0, reaction)

    assert failures == []
    assert pending == []
    assert len(snapshots) == 1
    product.do_optimization.assert_called_once_with()


def test_stage_four_drops_failed_optional_vdw_but_rejects_required_failure():
    ts = optimizer(shigh=-999, shir=-1, name='ts')
    vdw = optimizer(shigh=-999, shir=1, just_high=True, name='vdw')
    reaction = SimpleNamespace(
        ts_opt=ts,
        irc_prod_opt=vdw,
        irc_prod=SimpleNamespace(name='vdw'),
        prod_opt=[],
        do_vdW=True,
    )

    failures, pending, _ = generator()._poll_stage_four(0, reaction)

    assert failures == ['transition_state']
    assert pending == []
    assert reaction.do_vdW is False
    ts.do_optimization.assert_not_called()
    vdw.do_optimization.assert_not_called()


def test_stage_four_monitor_records_pending_hir_counts(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    pending = optimizer(shigh=1, shir=0)
    pending.species.hir = SimpleNamespace(hir_status=[[-1, 0, 1, 2]])

    snapshot = generator()._optimizer_snapshot('product_0', pending)
    generator()._write_optimization_monitor([{
        'reaction': 'test',
        'state': 4,
        'failures': [],
        'pending': ['product_0'],
        'optimizers': [snapshot],
    }])

    payload = json.loads((tmp_path / 'kinbot_optimization_monitor.json').read_text())
    assert payload['reactions'][0]['pending'] == ['product_0']
    assert payload['reactions'][0]['optimizers'][0]['hir_points'] == {
        'pending': 1,
        'complete': 1,
        'failed': 1,
        'skipped': 1,
    }
