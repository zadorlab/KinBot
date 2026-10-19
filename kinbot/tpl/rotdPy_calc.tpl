# The adapter resolves verified node-local scratch inside each QC worker.
# This controller may execute on a login node; it must not select worker paths.

calc = {{
'code': '{code}',
'scratch': '',
'method': '{method}',
'basis': {basis},
'mem': {mem},
'processors': {processors},
'queue': '{queue}',
'max_jobs': {max_jobs},
'max_retries': {max_retries},
'poll_interval': {poll_interval}
}}
