"""Update generated worker setup without replacing sampling owners or data."""
import json
from contextlib import contextmanager
import fcntl
from pathlib import Path
import shlex

from kinbot.execution import atomic_json, atomic_text, file_sha256

NAME = '.kinbot_worker_code.json'
MARKER = '# KinBot sampling worker code route\n'


@contextmanager
def _locked(folder):
    with (Path(folder)/'.kinbot_worker_code.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        yield


def validate(pin):
    pin = Path(pin).resolve(strict=True)
    hashes = json.loads((pin/'manifest.json').read_text())
    if 'kinbot/rotdpy_worker.py' not in hashes:
        raise ValueError('Sampling code snapshot lacks its worker.')
    for name, digest in hashes.items():
        p = pin/name
        if not p.resolve().is_relative_to(pin) or file_sha256(p) != digest:
            raise ValueError('Sampling code snapshot changed: '+str(p))
    return pin


def render(folder, base):
    receipt = Path(folder)/NAME
    if not receipt.exists():
        return base
    record = json.loads(receipt.read_text())
    pin = Path(record['code'])
    if file_sha256(pin/'manifest.json') != record['manifest_sha256']:
        raise ValueError('Sampling code manifest changed.')
    return base.split(MARKER)[0].rstrip()+'\n'+MARKER+(
        'export PYTHONPATH='+shlex.quote(str(pin))+':"${PYTHONPATH:-}"\n')


def publish(folder, base):
    # Replacement PES workers may regenerate the setup while a handoff installs
    # its route. Hold one lock across reading the receipt and replacing setup.
    with _locked(folder):
        return atomic_text(Path(folder)/'site_setup.sh', render(folder, base))


def preflight(root):
    root = Path(root).resolve(strict=True)
    plans = []
    prepared = json.loads((root/'prepared.json').read_text())
    for inp in prepared['inputs']:
        case = Path(inp).resolve().parent
        if not case.is_relative_to(root):
            raise ValueError('Sampling case is outside the prepared run.')
        folder = case/'rotdPy'
        setup = folder/'site_setup.sh'
        if not setup.exists():
            if list(folder.glob('*.campaign.json')):
                raise ValueError('Existing sampling campaign lacks generated setup: '+str(folder))
            continue
        text = setup.read_text()
        template = (folder/'qu.tpl').read_text()
        if (setup.is_symlink() or '# Generated from the environment visible during prepare.' not in text
                or shlex.quote(str(setup)) not in template or 'kinbot.rotdpy_worker' not in template):
            raise ValueError('Sampling workers do not use the generated setup: '+str(folder))
        plans.append((folder, setup, text))
    return plans


def install(root, pin):
    root = Path(root).resolve(strict=True)
    pin = validate(pin)
    if pin.parent != root or not pin.name.startswith('.salvage_code.'):
        raise ValueError('Expected this run\'s immutable source snapshot.')
    plans = preflight(root)
    # Validate every campaign before publishing any route. Old submitted scripts
    # source this same setup path, including templates cached by retained owners.
    for folder, setup, _ in plans:
        with _locked(folder):
            text = setup.read_text()
            archive = folder/'.worker_code_history'/file_sha256(setup)
            archive.mkdir(parents=True, exist_ok=True)
            if not (archive/'site_setup.sh').exists():
                atomic_text(archive/'site_setup.sh', text)
            atomic_json(folder/NAME, dict(schema=1, code=str(pin),
                manifest_sha256=file_sha256(pin/'manifest.json')))
            atomic_text(setup, render(folder, text))
    return [str(folder) for folder, _, _ in plans]


def main(argv=None):
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root'); p.add_argument('code')
    a = p.parse_args(argv)
    print(json.dumps(dict(worker_code=a.code, updated=install(a.root, a.code))))


if __name__ == '__main__':
    main()
