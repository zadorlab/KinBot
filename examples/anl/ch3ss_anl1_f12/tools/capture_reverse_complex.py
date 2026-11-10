from pathlib import Path
from datetime import datetime, timezone
from importlib.metadata import version, PackageNotFoundError
import hashlib
import json
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile

if len(sys.argv) != 4:
    raise SystemExit('Usage: capture_reverse_complex.py WORKER QC_JOB OUTPUT_PARENT')
worker = Path(sys.argv[1]).resolve(strict=True)
job = sys.argv[2]
if Path(job).name != job or not job.endswith(('_IRC_R_prod_high', '_IRC_F_prod_high')):
    raise SystemExit('QC_JOB must be an IRC endpoint L2 job name within WORKER.')
started = datetime.now(timezone.utc).isoformat()
parent = Path(sys.argv[3]).expanduser()
parent.mkdir(parents=True, exist_ok=True)
output = Path(tempfile.mkdtemp(prefix='reverse-complex-', dir=parent))
payload = output / 'capture'
payload.mkdir()
files, missing = [], []
print('Capturing this retry and its reaction context...', flush=True)

def take(source):
    source = Path(source)
    if not source.is_absolute():
        source = worker / source
    source = source.resolve()
    if not source.is_relative_to(worker):
        missing.append(str(source) + ' (outside worker)')
        return
    if not source.is_file():
        missing.append(str(source))
        return
    relative = source.relative_to(worker)
    target = payload / relative
    if target.exists():
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        for _ in range(3):
            before = source.stat()
            shutil.copyfile(source, target)
            after = source.stat()
            if (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns):
                break
    except OSError as error:
        if target.exists():
            target.unlink()
        missing.append(str(source) + ' (' + str(error) + ')')
        return
    files.append(dict(path=str(relative), source=str(source), bytes=target.stat().st_size,
        sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
        changed_during_copy=(before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns)))

for name in ('input.json', 'reverse_search_result.json', 'l1_endpoint_checks.json',
             'kinbot_status.json', 'kinbot.log', 'kinbot.err', 'job_ids.json',
             'kinbot_optimization_monitor.json'):
    take(name)
for source in worker.glob(job + '.*'):
    if source.suffix in ('.com', '.log', '.traj', '.py', '.json', '.sbatch',
                         '.exitcode', '.err', '.stdout', '.stderr'):
        take(source)
take(job + '_sella.log')
take(Path('perm') / (job + '.err'))
reaction = job.rsplit('_IRC_', 1)[0]
for stem in (reaction + '_high', reaction + '_IRC_F_prod_high',
             reaction + '_IRC_R_prod'):
    for suffix in ('.com', '.log', '.result.json'):
        take(stem + suffix)
for pattern in (reaction + '_IRC_R_prod_freq*.log', '*_well_high.log', '*_well_high.com'):
    for source in worker.glob(pattern):
        take(source)

token = hashlib.sha256(job.encode()).hexdigest()[:24]
receipt = worker / '.kinbot_timeout_retries' / token / 'current.json'
take(receipt)
if receipt.is_file():
    record = json.loads((payload / receipt.relative_to(worker)).read_text())
    for field in ('script', 'effective_worker'):
        if record.get(field):
            take(record[field])
    for attempt in record.get('attempts', []):
        for name in ('retry_worker.py', 'retry.sbatch', 'batch.original'):
            take(Path(attempt['archive']) / name)
    guard_token = hashlib.sha256(job.encode()).hexdigest()[:20]
    for job_id in [record['original_job_id'], *record.get('replacement_ids', [])]:
        take(Path('.kinbot_timeout_allocations') / (str(job_id) + '.json'))
        take(Path('.kinbot_slurm_io_retries') / guard_token / str(job_id) / 'execution.json')

database = worker / 'kinbot.db'
with sqlite3.connect(database.as_uri() + '?mode=ro', uri=True) as source:
    with sqlite3.connect(payload / 'kinbot.db') as target:
        source.backup(target)
db = payload / 'kinbot.db'
files.append(dict(path='kinbot.db', source=str(database), bytes=db.stat().st_size,
    sha256=hashlib.sha256(db.read_bytes()).hexdigest(), snapshot='SQLite read-only backup'))
versions = {}
for package in ('ase', 'sella', 'numpy', 'scipy'):
    try:
        versions[package] = version(package)
    except PackageNotFoundError:
        versions[package] = 'unavailable'
commit = subprocess.run(['git', 'rev-parse', 'HEAD'], text=True, capture_output=True).stdout.strip()
manifest = dict(worker=str(worker), job=job, commit=commit, versions=versions,
                files=files, unavailable_files=missing, started_utc=started,
                finished_utc=datetime.now(timezone.utc).isoformat())
(payload / 'manifest.json').write_text(json.dumps(manifest, indent=2))
archive = output / 'retry-evidence.tar.gz'
with tarfile.open(archive, 'w:gz', compresslevel=1) as tar:
    tar.add(payload, arcname='capture')
with tarfile.open(archive) as tar:
    for entry in files:
        data = tar.extractfile('capture/' + entry['path']).read()
        assert hashlib.sha256(data).hexdigest() == entry['sha256'], entry['path']
print('ARCHIVE', archive)
print('SHA256', hashlib.sha256(archive.read_bytes()).hexdigest())
print('FILES', len(files), 'BYTES', archive.stat().st_size)
print('LIVE_FILES_CHANGED', [f['path'] for f in files if f.get('changed_during_copy')])
print('MISSING_CURRENT_CORE', [job + suffix for suffix in ('.com', '.log', '.traj')
                             if not (payload / (job + suffix)).is_file()])
