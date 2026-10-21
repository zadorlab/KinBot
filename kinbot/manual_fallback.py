"""Explicit next-attempt fallbacks; never edit running workers or waive parsers.

The request embeds the JSON payload. Native setup commands are operator code,
executed only in the replacement allocation, never by the client/controller.
"""
import ast
from copy import deepcopy
import hashlib
import json
import re
import shlex
from pathlib import Path


def validate(value, kind):
    from kinbot.live_control import validate_resources
    common = {'schema', 'resources', 'environment', 'setup', 'policy'}
    allowed = common | ({'input_edits', 'command', 'calculator_kwargs', 'optimizer'}
                        if kind == 'task' else {'worker_parameters'})
    if not isinstance(value, dict) or value.get('schema') != 1 or set(value)-allowed:
        raise ValueError(f'Invalid {kind} manual fallback fields/schema.')
    validate_resources(value.get('resources', {}))
    env = value.get('environment', {})
    if not isinstance(env, dict):
        raise ValueError('Fallback environment must be an object.')
    for key, item in env.items():
        if (not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*', key)
                or key in {'HOME','USER','LOGNAME','UID','PWD','PYTHONPATH','PYTHONHOME'}
                or key.startswith(('SLURM_', 'SBATCH_'))
                or (key.startswith('KINBOT_') and not re.fullmatch(
                    r'KINBOT_(?:MOLPRO|GAUSSIAN|MRCC|CFOUR|FORMCHK|CUBEGEN)_(?:ROOT|EXECUTABLE|SETUP)', key))):
            raise ValueError(f'Fallback cannot replace controller/scheduler identity: {key}.')
        if item is not None and (not isinstance(item, str) or '\x00' in item):
            raise ValueError('Environment values must be strings or null (unset).')
    setup = value.get('setup', [])
    if not isinstance(setup, list) or any(not isinstance(s,str) or not s.strip()
            or any(c in s for c in '\n\r\x00') for s in setup):
        raise ValueError('Fallback setup must contain single shell command lines.')
    policy = value.get('policy', {})
    if (not isinstance(policy,dict) or set(policy)-{'optimization_step_retries'}
            or any(type(v) is not int or not 0 <= v <= 10 for v in policy.values())):
        raise ValueError('Fallback policy supports optimization_step_retries from 0 to 10.')
    edits = value.get('input_edits', [])
    if not isinstance(edits,list) or any(not isinstance(e,dict) or set(e)!={'old','new'}
            or not isinstance(e['old'],str) or not e['old']
            or not isinstance(e['new'],str) for e in edits):
        raise ValueError('Native input_edits require exact nonempty old and new strings.')
    for field in ('calculator_kwargs','optimizer','worker_parameters'):
        if field in value and not isinstance(value[field],dict):
            raise ValueError(f'{field} must be an object.')
    if set(value.get('optimizer',{}))-{'steps','fmax','sella_kwargs','minimum_recovery'}:
        raise ValueError('Fallback cannot change stationary-point order or constraints.')
    for options in (value.get('optimizer',{}).get('sella_kwargs',{}),
                    value.get('worker_parameters',{}).get('sella_kwargs',{})):
        if not isinstance(options,dict) or set(options)&{'order','constraints','atoms','trajectory','logfile','restart'}:
            raise ValueError('Sella fallback cannot replace geometry, order, constraints or result files.')
    for name, v in value.get('worker_parameters',{}).items():
        if not re.fullmatch(r'(?:kwargs\d*|prod_kwargs|sella_kwargs|steps|fmax)', name):
            raise ValueError(f'Unsupported generated-worker parameter: {name}.')
        if name in ('steps','fmax'):
            if isinstance(v,bool) or not isinstance(v,(int,float)) or v <= 0 or (name=='steps' and type(v)!=int):
                raise ValueError('steps/fmax must be positive; steps must be integer.')
        elif not isinstance(v,dict):
            raise ValueError(f'{name} must be a keyword object.')
    forbidden = {'method','basis','charge','mult','multiplicity','symbols','positions','atoms',
                 'geometry','restart','label','directory','chk','oldchk','nprocshared','mem','nproc','stack_mw'}
    for values in [value.get('calculator_kwargs',{}), *(v for k,v in value.get('worker_parameters',{}).items()
                   if k not in ('steps','fmax','sella_kwargs'))]:
        if {k.lower() for k in values} & forbidden:
            raise ValueError('Use resources for CPUs/memory; theory, geometry and checkpoint identity need a scientific migration.')
    if 'command' in value and (not isinstance(value['command'],list) or not value['command']
            or any(not isinstance(s,str) or not s or '\x00' in s for s in value['command'])):
        raise ValueError('Native command must be a nonempty argument list.')
    json.dumps(value,allow_nan=False)
    return deepcopy(value)


def shell_lines(value):
    lines = list(value.get('setup', []))
    for key, item in value.get('environment', {}).items():
        lines.append('unset '+key if item is None else 'export '+key+'='+shlex.quote(item))
    if 'optimization_step_retries' in value.get('policy', {}):
        lines.append('export KINBOT_OPTIMIZATION_STEP_RETRIES='+str(value['policy']['optimization_step_retries']))
    return lines


def digest(task):
    item=deepcopy(task)
    item.pop('resources',None)
    item.pop('manual_fallback',None)
    return hashlib.sha256(json.dumps(item,sort_keys=True,allow_nan=False).encode()).hexdigest()


def logical_task(task):
    marker=task.get('manual_fallback')
    if marker is None:
        return task
    if digest(task)!=marker['patched_sha256'] or digest(marker['original'])!=marker['original_sha256']:
        raise ValueError('Manual fallback task changed outside its approved receipt.')
    return deepcopy(marker['original'])


def task_options(task, value):
    value=validate(value,'task')
    original=logical_task(task)
    result=deepcopy(task)
    result.pop('manual_fallback',None)
    result['resources'].update(value.get('resources',{}))
    if value.get('resources'):
        result['resources']['live_override']=True
        if 'memory_mb' in value['resources']:
            result['resources']['use_all_node_memory']=False
    if result['kind']=='external':
        if value.get('calculator_kwargs') or value.get('optimizer'):
            raise ValueError('External native tasks accept input_edits, not ASE options.')
        for edit in value.get('input_edits',[]):
            if result['input_template'].count(edit['old']) != 1:
                raise ValueError('Each native input edit must match exactly once before cancellation.')
            result['input_template']=result['input_template'].replace(edit['old'],edit['new'])
        if 'command' in value:
            result['command']=value['command']
    else:
        if value.get('input_edits') or 'command' in value:
            raise ValueError('ASE tasks accept calculator_kwargs (including command), not native-deck edits.')
        kwargs=result['profile'].setdefault('calculator_kwargs',{})
        for k,v in value.get('calculator_kwargs',{}).items():
            if k == 'command':
                if not isinstance(v,str) or not v.strip():
                    raise ValueError('ASE profile command must be a nonempty string.')
                result['profile']['command']=v
                continue
            if v is None: kwargs.pop(k,None)
            else: kwargs[k]=v
        for k,v in value.get('optimizer',{}).items():
            if isinstance(v,dict): result['optimizer'].setdefault(k,{}).update(v)
            else: result['optimizer'][k]=v
    # Operator commands must never run in login-host environment preflight.
    result['runtime_setup']=[*result.get('runtime_setup',[]),*shell_lines(value)]
    if value.get('resources'):
        from kinbot.anl.live_control import validate_native_resources
        validate_native_resources(result,value['resources'])
    effective=deepcopy(value)
    effective['resources']={**task.get('manual_fallback',{}).get('options',{}).get('resources',{}),
                            **value.get('resources',{})}
    result['manual_fallback']=dict(original=original, original_sha256=digest(original),
                                   patched_sha256=digest(result), options=effective)
    return result


def worker_options(text, value):
    """Edit literal parameter assignments only, never execute a saved worker."""
    validate(value,'qc')
    requested=value.get('worker_parameters',{})
    source=text.encode('utf8'); offsets=[0]
    for line in source.splitlines(keepends=True): offsets.append(offsets[-1]+len(line))
    edits=[]; seen=set()
    for node in ast.parse(text).body:
        if not (isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name)):
            continue
        name=node.targets[0].id
        if name not in requested: continue
        if name in seen: raise ValueError(f'Ambiguous generated-worker assignment: {name}.')
        seen.add(name)
        old=ast.literal_eval(node.value); new=requested[name]
        if isinstance(new,dict):
            if not isinstance(old,dict): raise ValueError(f'{name} is not a literal keyword dictionary.')
            old.update({k:v for k,v in new.items() if v is not None})
            for k,v in new.items():
                if v is None: old.pop(k,None)
            new=old
        edits.append((offsets[node.value.lineno-1]+node.value.col_offset,
                      offsets[node.value.end_lineno-1]+node.value.end_col_offset,repr(new).encode('utf8')))
    if seen!=set(requested): raise ValueError('Worker has no literal assignments for: '+', '.join(sorted(set(requested)-seen)))
    for start,end,new in sorted(edits,reverse=True): source=source[:start]+new+source[end:]
    result=source.decode('utf8');compile(result,'<manual-fallback-worker>','exec')
    return result


def worker_script(text, worker_command, replacement, value):
    """Place runtime setup immediately before the uniquely identified worker."""
    pattern=re.compile(r'(?<![\w/.-])(?:"'+re.escape(worker_command)+r'"|\''+re.escape(worker_command)+r'\'|'+re.escape(worker_command)+r')(?=\s|$)')
    if len(pattern.findall(text))!=1:
        raise ValueError('Cannot uniquely identify the Python worker in its batch script.')
    lines=text.splitlines(keepends=True)
    for i,line in enumerate(lines):
        if pattern.search(line):
            if line.lstrip().startswith('#'):
                raise ValueError('Worker command occurs in a comment, not an executable line.')
            lines[i]=''.join(s+'\n' for s in shell_lines(value))+pattern.sub(lambda _:shlex.quote(str(replacement)),line)
            break
    return ''.join(lines)


def reaction_workers(target):
    target=Path(target).resolve()
    return [p for p in target.parent.glob(target.name+'*.py') if
            p.stem==target.name or p.stem.startswith(target.name+'_')]


def install_reaction(target, request):
    """Owner-only scope for newly generated attempts of one rejected search."""
    from kinbot.execution import atomic_json
    target=Path(target).resolve()
    value=validate(request['fallback'],'reaction')
    key=hashlib.sha256(str(target).encode()).hexdigest()
    atomic_json(target.parent/'.kinbot_reaction_fallbacks'/(key+'.json'),
                dict(target=str(target),request=request['id'],fallback=value))


def for_reaction_job(target):
    target=Path(target).resolve()
    matches=[]
    for p in (target.parent/'.kinbot_reaction_fallbacks').glob('*.json'):
        record=json.loads(p.read_text());reaction=Path(record['target'])
        if (reaction.parent==target.parent and
                (target.name==reaction.name or target.name.startswith(reaction.name+'_'))):
            matches.append((len(reaction.name),record))
    return max(matches,key=lambda x:x[0])[1]['fallback'] if matches else {}
