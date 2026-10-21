"""Freeze resource rules in new MESS/ROTD allocations, not scientific inputs."""
import json
from pathlib import Path
import re
import shlex

from kinbot import live_control
from kinbot.execution import atomic_json
from kinbot.qc_live import script_resources


def callback(script, *, kind):
    script = Path(script).resolve()
    if live_control.find(script.parent) is None:
        return None
    original = script.read_text()
    target = script.with_suffix('')
    if kind == 'rotd':
        match = re.search(r'-m kinbot\.rotdpy_worker (surf\d+_face\d+_samp\d+)\s*(?:\n|$)',original)
        if not match:
            return None  # Custom worker: never guess its native settings.
        target = script.parent/match[1]
    def prepare():
        options = live_control.overrides(target)
        contents = original
        frozen = dict(options)
        if kind == 'rotd' and (set(options)-{'walltime'}):
            from kinbot.anl.dispatch import _molpro_stack_mw
            for key,pattern in [('cores',r'(?m)^#SBATCH --ntasks=(\d+)\s*$'),
                                ('memory_mb',r'(?m)^#SBATCH --mem=(\d+)M\s*$')]:
                match = re.search(pattern,original)
                if key not in frozen:
                    if not match: raise ValueError('ROTD resource edit needs explicit cores and memory_mb for this custom template.')
                    frozen[key] = int(match[1])
            frozen['stack_mw'] = _molpro_stack_mw(frozen)
            contents = script_resources(original,{k:v for k,v in frozen.items() if k!='stack_mw'})
            contents = contents.replace('#SBATCH --ntasks=1\n',f'#SBATCH --ntasks={frozen["cores"]}\n')
            contents = re.sub(r'(?m)^#SBATCH --cpus-per-task=\d+$','#SBATCH --cpus-per-task=1',contents)
            contents = re.sub(r'(?m)^export OMP_NUM_THREADS=\d+$','export OMP_NUM_THREADS=1',contents)
            rows=contents.splitlines(keepends=True)
            offset=next((i for i,s in enumerate(rows) if s.strip() and not s.lstrip().startswith('#')),len(rows))
            rows.insert(offset,'export KINBOT_ROTD_RESOURCES='+shlex.quote(json.dumps(frozen))+'\n')
            contents=''.join(rows)
        elif options:
            contents = script_resources(original,options)
        if contents != script.read_text():
            script.write_text(contents)
        atomic_json(script.with_suffix('.resources.json'),dict(target=str(target),kind=kind,
                    resources=frozen,updated_at=live_control.now()))
    return prepare
