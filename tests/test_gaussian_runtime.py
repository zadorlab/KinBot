"""Gaussian runtime selection must survive invalid site scratch defaults."""

import os
from pathlib import Path

from kinbot.ase_modules.calculators.gaussian import gaussian_scratch_environment


def test_invalid_gaussian_scratch_falls_back_to_site_scratch(tmp_path,
                                                              monkeypatch):
    blocked = tmp_path / 'blocked'
    blocked.write_text('not a directory')
    site = tmp_path / 'site-scratch'
    monkeypatch.setenv('GAUSS_SCRDIR', str(blocked))
    monkeypatch.delenv('SLURM_TMPDIR', raising=False)
    monkeypatch.setenv('SCRATCH', str(site))
    monkeypatch.delenv('TMPDIR', raising=False)

    with gaussian_scratch_environment() as (scratch, source):
        assert source == 'SCRATCH'
        assert scratch.parent == site
        assert scratch.is_dir()
        assert os.environ['GAUSS_SCRDIR'] == str(scratch)

    assert os.environ['GAUSS_SCRDIR'] == str(blocked)
    assert not scratch.exists()


def test_scheduler_local_scratch_precedes_generic_site_scratch(tmp_path,
                                                                monkeypatch):
    node = tmp_path / 'node-scratch'
    site = tmp_path / 'site-scratch'
    monkeypatch.delenv('GAUSS_SCRDIR', raising=False)
    monkeypatch.setenv('SLURM_TMPDIR', str(node))
    monkeypatch.setenv('SCRATCH', str(site))

    with gaussian_scratch_environment() as (scratch, source):
        assert source == 'SLURM_TMPDIR'
        assert scratch.parent == node

    assert 'GAUSS_SCRDIR' not in os.environ
    assert not scratch.exists()
