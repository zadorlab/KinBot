scratch = Path(os.environ.get('SCRATCH') or
               os.environ.get('SLURM_TMPDIR') or
               (Path.home() / '.kinbot_scratch'))
scratch.mkdir(parents=True, exist_ok=True)

calc = {{
'code': '{code}',
'scratch': str(scratch),
'method': '{method}',
'basis': '{basis}',
'mem': {mem},
'processors': {processors},
'queue': '{queue}',
'max_jobs': {max_jobs},
'max_retries': {max_retries},
'poll_interval': {poll_interval}
}}
