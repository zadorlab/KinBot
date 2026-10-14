kinbot_exit_file={name}.exitcode
kinbot_exit_tmp="$kinbot_exit_file.tmp.$$"
rm -f "$kinbot_exit_file" "$kinbot_exit_tmp"
set +e
{python_executable} {python_file} {arguments}
kinbot_exit_code=$?
set -e
printf '%s %s\n' "${{SLURM_JOB_ID:-unknown}}" "$kinbot_exit_code" > "$kinbot_exit_tmp"
mv "$kinbot_exit_tmp" "$kinbot_exit_file"
exit "$kinbot_exit_code"
