"""Live resource overrides and explicit clean retries for ordinary QC workers."""
import ast
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import re

from kinbot import live_control as control
from kinbot.execution import atomic_json, file_sha256


def arguments(kwargs, backend, options):
    result = deepcopy(kwargs)
    if not options:
        return result
    if backend in ('gauss', 'gaussian'):
        if 'cores' in options:
            result['nprocshared'] = options['cores']
        if 'memory_mb' in options:
            result['mem'] = f"{max(1, int(options['memory_mb']*.7))}MB"
    elif backend == 'orca':
        blocks = result.get('orcablocks', '')
        if 'cores' in options:
            blocks = re.sub(r'%pal\s+nprocs\s+\d+\s+end', '', blocks, flags=re.I)
            blocks += f" %pal nprocs {options['cores']} end"
        if 'memory_mb' in options:
            match = re.search(r'%pal\s+nprocs\s+(\d+)\s+end', blocks, re.I)
            if not match:
                raise ValueError('ORCA memory edits require an explicit PAL rank count.')
            blocks = re.sub(r'%maxcore\s+\d+', '', blocks, flags=re.I)
            blocks += f" %maxcore {max(1,int(options['memory_mb']*.7/int(match[1])))}"
        result['orcablocks'] = blocks
    elif set(options)-{'walltime'}:
        raise ValueError(f'Native resource edits for ordinary {backend} workers are not implemented; use a dispatcher task.')
    return result


def worker_resources(text, backend, options):
    """Change only literal calculator keyword dictionaries, never chemistry."""
    if not set(options)-{'walltime'}:
        return text
    tree = ast.parse(text)
    source = text.encode("utf-8")
    lines = source.splitlines(keepends=True)
    offsets = [0]
    for line in lines:
        offsets.append(offsets[-1]+len(line))
    edits = []
    for node in tree.body:
        if (isinstance(node, ast.Assign) and len(node.targets)==1
                and isinstance(node.targets[0], ast.Name)
                and re.fullmatch(r'kwargs\d*', node.targets[0].id)):
            old = ast.literal_eval(node.value)
            if not isinstance(old, dict):
                raise ValueError('QC resource editing requires literal calculator kwargs.')
            edits.append((offsets[node.value.lineno-1]+node.value.col_offset,
                          offsets[node.value.end_lineno-1]+node.value.end_col_offset,
                          repr(arguments(old, backend, options)).encode("utf-8")))
    if not edits:
        raise ValueError('Generated QC worker has no supported calculator kwargs; no resource change made.')
    for start,end,value in sorted(edits, reverse=True):
        source = source[:start]+value+source[end:]
    text = source.decode("utf-8")
    compile(text, '<resource-updated worker>', 'exec')
    return text


def script_resources(text, options):
    """Keep custom site setup/exclusions, but make allocation and threads agree."""
    if not options:
        return text
    import shlex
    lines = text.splitlines(keepends=True)
    remove = set()
    if 'cores' in options:
        remove |= {'-c','--cpus-per-task','-n','--ntasks','--ntasks-per-node'}
    if 'memory_mb' in options:
        remove |= {'--mem','--mem-per-cpu','--mem-per-gpu'}
    if 'walltime' in options:
        remove |= {'-t','--time'}
    header, body = [], []
    in_header = True
    for line in lines:
        if in_header and line.strip() and not line.lstrip().startswith('#'):
            in_header = False
        if in_header and line.lstrip().startswith('#SBATCH'):
            words = shlex.split(line.lstrip()[len('#SBATCH'):])
            if len(words) > 2 or (len(words)==2 and words[1].startswith('-')):
                raise ValueError('Live resource edits require one SBATCH option per line.')
            if words and (words[0].split('=')[0] in remove or
                          any(words[0].startswith(x) and words[0][len(x):].isdigit()
                              for x in remove if len(x)==2)):
                continue
        (header if in_header else body).append(line)
    if 'cores' in options:
        header += [f"#SBATCH --ntasks=1\n#SBATCH --cpus-per-task={options['cores']}\n"]
        body = [re.sub(r'^(export\s+OMP_NUM_THREADS=).*$', rf'\g<1>{options["cores"]}', s.rstrip('\n'))+'\n'
                for s in body]
        header += [f"export OMP_NUM_THREADS={options['cores']}\n"]
    if 'memory_mb' in options:
        # All directives must precede the first export above.
        header.insert(next((i for i,s in enumerate(header) if s.startswith('export ')), len(header)),
                      f"#SBATCH --mem={options['memory_mb']}M\n")
    if 'walltime' in options:
        header.insert(next((i for i,s in enumerate(header) if s.startswith('export ')), len(header)),
                      f"#SBATCH --time={options['walltime']}\n")
    return ''.join(header+body)


def prepare_resources(qc, job, script, original_worker=None):
    options = control.overrides(Path.cwd()/job)
    if not options and original_worker is None:
        return script
    worker = Path(job+'.py')
    source = worker.read_text() if original_worker is None else original_worker
    updated = worker_resources(source, qc.qc, options) if options else source
    rewritten = script_resources(script, options)
    # Keep the original worker for audit; this is called only before sbatch.
    token = hashlib.sha256(job.encode()).hexdigest()[:24]
    archive = Path('.kinbot_controls/qc_resources')/token/file_sha256(worker)
    if not archive.exists():
        archive.parent.mkdir(parents=True, exist_ok=True)
        archive.write_bytes(worker.read_bytes())
    temporary = worker.with_name(worker.name+'.resources.tmp')
    temporary.write_text(updated)
    temporary.replace(worker)
    atomic_json(Path('.kinbot_controls/qc_resources')/(token+'.json'),
                dict(job=job, resources=options, source=str(archive), updated_at=control.now()))
    return rewritten


def process(qc, job, status, *, kinds=None, allow_completed=False):
    from kinbot.pool_drain import pulse
    pulse(Path.cwd())
    from kinbot import qc_timeout
    from kinbot.anl.live_control import stop_attempt
    target = (Path.cwd()/job).resolve()
    if job in qc.job_ids:
        control.observe_qc(target,qc.job_ids[job],status)
    for path, r in control.requests(Path.cwd(), kinds=kinds or {'qc'}):
        if Path(r['target']) != target:
            continue
        try:
            current = control.generation(target, r['kind'])
            if any(current.get(k)!=v for k,v in r['expected'].items() if k!='exitcode'):
                raise ValueError('QC input/attempt changed after this request.')
            if r['action'] == 'release':
                control.acknowledge(path, r, 'applied')
                continue
            if status not in ('error', 'running'):
                raise ValueError('Only a failed QC job or explicitly stopped allocation can be retried.')
            if job not in qc.job_ids:
                raise ValueError('No owned Slurm ID for this QC attempt; reconcile before retrying.')
            if r['expected'].get('job_id') != str(qc.job_ids[job]):
                raise ValueError('QC allocation generation is unknown or changed; inspect the current target and request again.')
            # Validate native resource rendering before cancellation.
            worker_resources(Path(job+'.py').read_text(), qc.qc, control.overrides(target))
            from kinbot import run_policy
            if run_policy.find() is None:
                raise ValueError('Explicit QC retries require shared admission receipts; enable a queue policy first.')
            token = hashlib.sha256(job.encode()).hexdigest()[:20]
            guard = Path('.kinbot_slurm_io_retries')/token/str(qc.job_ids[job])/'execution.json'
            if not r['pending_only']:
                if not guard.exists():
                    raise ValueError('Original allocation has no verified retry seeds; refusing to cancel it.')
                saved = json.loads(guard.read_text())
                if saved.get('identity') != job or str(saved.get('slurm_job_id')) != str(qc.job_ids[job]):
                    raise ValueError('Retry seed generation does not match.')
                for name,digest in saved['initial_artifacts'].items():
                    seed = guard.parent/'initial'/name
                    if not seed.resolve().is_relative_to(guard.parent.resolve()) or file_sha256(seed) != digest:
                        raise ValueError('Original allocation retry seed changed; refusing to cancel it.')
            if not stop_attempt(Path.cwd(), qc.job_ids[job], path, r,
                                allow_completed=allow_completed or status=='error'):
                return status
            control.acknowledge(path, r, 'applying')
            r['request_file'] = str(path)
            result = qc_timeout.recover(qc, job, manual=r)
            if result is None:
                raise ValueError('QC retry could not establish a replacement; existing work retained.')
            return result
        except (ValueError, RuntimeError, OSError) as error:
            control.acknowledge(path, r, 'blocked' if isinstance(error, ValueError) else r['status'], error=str(error))
    if status == 'error' and kinds is None and control.holding(target, kind='qc'):
        return 'running'
    return status
