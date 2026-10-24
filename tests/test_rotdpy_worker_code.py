import os
from pathlib import Path
import subprocess
import sys

import pytest
from kinbot import rotdpy_worker_code as route
from kinbot.execution import atomic_json, file_sha256


def test_cached_template_gets_new_worker_without_changing_sampling_data(tmp_path):
    root = tmp_path/'run'; folder = root/'case/rotdPy'; folder.mkdir(parents=True)
    atomic_json(root/'prepared.json', dict(inputs=[str(root/'case/input.json')]))
    setup = folder/'site_setup.sh'
    base = '# Generated from the environment visible during prepare.\nexport SAMPLE_NATIVE=retained\n'
    setup.write_text(base)
    cached = 'source '+str(setup)+'\n'+sys.executable+' -B -m kinbot.rotdpy_worker\n'
    (folder/'qu.tpl').write_text(cached)
    pin = root/'.salvage_code.test'; (pin/'kinbot').mkdir(parents=True)
    (pin/'kinbot/__init__.py').write_text('')
    (pin/'kinbot/rotdpy_worker.py').write_text('import os; print("new worker",os.environ["SAMPLE_NATIVE"])')
    atomic_json(pin/'manifest.json', {str(p.relative_to(pin)):file_sha256(p) for p in (pin/'kinbot').iterdir()})
    for name in ('owner.pkl','samples.db','sample.py'):
        (folder/name).write_bytes(b'saved')
    for _ in range(2):
        route.install(root,pin)
        route.publish(folder, base)
        assert setup.read_text().count(route.MARKER) == 1
        assert route.render(folder,base) == setup.read_text()
        result = subprocess.run(['bash','-c',cached],cwd=folder,env=dict(os.environ,PYTHONPATH='/old/code'),
                                text=True,capture_output=True,check=True)
        assert result.stdout.strip() == 'new worker retained'
        assert all((folder/n).read_bytes() == b'saved' for n in ('owner.pkl','samples.db','sample.py'))
    # Interleaving replacement-PES setup publication with handoff publication
    # cannot drop the route after either writer has read a stale receipt.
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(route.install, root, pin) if i%2 else
                   executor.submit(route.publish, folder, base) for i in range(20)]
        for future in futures:
            future.result()
    assert setup.read_text() == route.render(folder, base)
    (pin/'kinbot/rotdpy_worker.py').write_text('changed')
    with pytest.raises(ValueError,match='changed'):
        route.install(root,pin)
