import pytest
from kinbot.qc import exclusive_slurm_script
from kinbot.parameters import Parameters


def test_exclusive_defaults_apply_to_ordinary_and_sampling_jobs(tmp_path):
    source = tmp_path/'input.json'
    source.write_text('{"barrier_threshold": 100}')
    par = Parameters(str(source)).par
    assert par['queue_exclusive'] == par['rotdpy_exclusive'] == 1


@pytest.mark.parametrize('header', [
    '#!/bin/bash\n#SBATCH -c 8\n#SBATCH --mem=10G\necho work\n',
    '#!/bin/bash\necho work\n',
    '#SBATCH -c 8\necho work\n'])
def test_exclusivity_precedes_executable_lines_and_keeps_resources(header):
    result = exclusive_slurm_script(header)
    assert result.index('#SBATCH --exclusive') < result.index('echo work')
    assert result.replace('#SBATCH --exclusive\n', '', 1) == header
    assert exclusive_slurm_script(header, False) == header


@pytest.mark.parametrize('option', ['--oversubscribe', '--share', '-s', '--exclusive=user', '--exclusive=mcs'])
def test_conflicting_custom_header_does_not_silently_weaken_exclusivity(option):
    with pytest.raises(ValueError, match='conflicts'):
        exclusive_slurm_script(f'#!/bin/bash\n#SBATCH -c 8 {option}\necho work\n')
