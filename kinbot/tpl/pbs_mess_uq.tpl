cd "${{PBS_O_WORKDIR}}/me" || exit 1
if mess mess_{n}.inp; then
    mess_exit_code=0
else
    mess_exit_code=$?
fi
printf '%s\n' "$mess_exit_code" > mess_{n}.exitcode
exit "$mess_exit_code"
