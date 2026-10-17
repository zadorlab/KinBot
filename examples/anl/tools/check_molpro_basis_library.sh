#!/usr/bin/env bash
# Standalone Molpro 2024 library probe for the CH3S/CH3SS production workflows.
# No KinBot imports, SSH, scheduler submission, SCF, or correlation calculation.
# INT with GDIRECT forces basis loading and one-electron setup without storing
# the four-index electron-repulsion integrals. Auxiliary definitions are loaded
# as the probe's orbital basis: this checks their data, not F12 driver selection.
# References: Molpro manual, basis_input and integration; augmented F12 family:
# Sylvetsky et al., JCP 147, 134106 (2017), DOI 10.1063/1.4998332, Table II.
set -eo pipefail
module load molpro/molpro24
command -v molpro
command -v timeout >/dev/null
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MOLPRO_OPTIONS=''

KB_BASIS_REPORT=$(mktemp -d "$HOME/kinbot-basis-check.XXXXXX")
KB_BASIS_SCRATCH=$(mktemp -d "${TMPDIR:-/tmp}/kinbot-basis-check.XXXXXX")
printf 'Report: %s\nScratch: %s\n' "$KB_BASIS_REPORT" "$KB_BASIS_SCRATCH"
module list > "$KB_BASIS_REPORT/modules.txt" 2>&1
timeout -k 5s 15s molpro --version < /dev/null > "$KB_BASIS_REPORT/version.txt" 2>&1 || true
cat "$KB_BASIS_REPORT/version.txt"
printf 'ROLE\tBASIS\tSTATUS\tOUTPUT\n' > "$KB_BASIS_REPORT/summary.tsv"

KB_BASIS_INDEX=0
while read -r role basis; do
    [[ -z "$role" || "$role" == \#* ]] && continue
    KB_BASIS_INDEX=$((KB_BASIS_INDEX + 1))
    case_dir="$KB_BASIS_REPORT/$KB_BASIS_INDEX"
    scratch_dir="$KB_BASIS_SCRATCH/$KB_BASIS_INDEX"
    mkdir -p "$case_dir" "$scratch_dir"
    cat > "$case_dir/load.inp" <<EOF
***,Basis library availability only
memory,64,m
symmetry,nosym
orient,noorient
angstrom
geometry={
3
All elements in the two PES and their CBH fragments
H 0.0 0.0 0.0
C 10.0 0.0 0.0
S 0.0 10.0 0.0
}
basis=$basis
gprint,basis
gdirect
int
kbload_ok=1
show,kbload_ok
---
EOF
    rc=0
    (cd "$case_dir" && timeout -k 5s 30s molpro -n 1 -m 64 \
        -d "$scratch_dir" -I "$scratch_dir" -W "$scratch_dir" load.inp \
        < /dev/null > launcher.log 2>&1) || rc=$?
    if [[ "$rc" == 0 ]] &&
       grep -Eq '^[[:space:]]*KBLOAD_OK[[:space:]]*=[[:space:]]*1([.]0+)?([[:space:]]|$)' "$case_dir/load.out"; then
        status=LOADED
    elif [[ "$rc" == 124 || "$rc" == 137 ]]; then
        status=TIMEOUT
    else
        status="ERROR(rc=$rc)"
    fi
    printf '%-11s %-36s %s\n' "$role" "$basis" "$status"
    printf '%s\t%s\t%s\t%s\n' "$role" "$basis" "$status" "$case_dir/load.out" \
        >> "$KB_BASIS_REPORT/summary.tsv"
    if [[ "$status" != LOADED ]]; then
        if [[ -f "$case_dir/load.out" ]]; then
            grep -Ei -A 2 'BASIS.*(NOT FOUND|NOT AVAILABLE)|CANNOT FIND|UNKNOWN|ILLEGAL|\? *ERROR|problem occurs|GLOBAL ERROR' \
                "$case_dir/load.out" | head -n 14 || true
            tail -n 6 "$case_dir/load.out"
        else
            tail -n 8 "$case_dir/launcher.log"
        fi
    fi
done <<'BASIS_LIST'
# Exact orbital names requested by the current ANL1-F12 and VRC inputs.
ORBITAL cc-pVDZ
ORBITAL cc-pVTZ
ORBITAL cc-pVQZ
ORBITAL aug-cc-pV5Z
ORBITAL aug-cc-pV6Z
ORBITAL cc-pCVTZ
ORBITAL cc-pCVQZ
ORBITAL aug-cc-pCVTZ-DK
ORBITAL avdz-f12
ORBITAL avtz-f12
ORBITAL aug-cc-pVQZ-F12
ORBITAL aug-cc-pV5Z-F12
# F12 auxiliary library lookups. These do not run the F12 driver.
AUXILIARY vdz-f12/jkfit
AUXILIARY avdz-f12/mp2fit
AUXILIARY cc-pVDZ-F12/optri
AUXILIARY vtz-f12/jkfit
AUXILIARY avtz-f12/mp2fit
AUXILIARY cc-pVTZ-F12/optri
AUXILIARY vqz-f12/jkfit
AUXILIARY avqz-f12/mp2fit
AUXILIARY cc-pVQZ-F12/optri
# Probe 5Z mappings and the large MP2Fit CABS separately from orbital coverage.
AUX-5Z v5z-f12/jkfit
AUX-5Z av5z-f12/mp2fit
AUX-5Z aug-cc-pwCV5Z/mp2fit
# Diagnostics only: no fallback or production basis change is made.
DIAGNOSTIC av5z-f12
DIAGNOSTIC cc-pV5Z-F12
DIAGNOSTIC cc-pV5Z-F12/optri
BASIS_LIST

if [[ "$KB_BASIS_INDEX" != 27 ]]; then
    printf 'INCOMPLETE: attempted %s of 27 lookups.\n' "$KB_BASIS_INDEX" >&2
    exit 1
fi
printf '\nCompleted all %s lookups.\n' "$KB_BASIS_INDEX"
printf '\nSummary: %s/summary.tsv\n' "$KB_BASIS_REPORT"
printf 'Inputs, native outputs and launcher logs retained; no existing runs changed.\n'
printf 'LOADED means H/C/S basis data loaded, not electronic-structure validation.\n'
printf 'A missing diagnostic 5Z/OPTRI is not by itself a missing 5Z orbital basis.\n'
