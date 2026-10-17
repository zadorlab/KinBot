# CASPT2 production transition, 16 October 2026

The user selected CASPT2-F12 for both sampling and the multireference spin
gaps. MRCI remains an optional benchmark; the already-submitted build trial
is left alone. Its unresolved numerical problem is not a production gate.
Minimum spaces remain CH3S C-S (5,4), CH3SS C-S (7,5), CH3SS S-S (7,5).
Keep IPEA=0.25 and set the real CASPT2 shift to 0.1 Eh. Preserve all native
evidence; do not relabel old shift-0.3 results as shift-0.1 calculations.

## Matched shift check: ready to submit

The single geometry is the reviewed CH3SS S-S inner point (2.1 angstrom).
It has the largest shifted first-order total norm in the completed shared-
reference checks: 1.16660671 at aTZ and 1.15169295 at aDZ (SHIFT=0.3).
This is an evidence-based sensitivity proxy, not a proof of which uncomputed
geometry is globally most sensitive. Its fractional active occupations also
make it a useful interacting-state check. Run both bases at this one geometry
because the requested production correction extrapolates their spin gaps.

The plan keeps all six doublet and six quartet roots, the equal same-spin
reference density, IPEA, geometry, basis, core treatment and solver thresholds.
Each job starts from the same completed center checkpoint that fed its old
inner calculation. It makes no additional orbital permutation. The source
is the stopped recovery campaign, independently bound to the deposited audit.
No other campaign is stopped, edited or inspected for cancellation.

```bash
cd /home/hram/KinBot
git pull --ff-only origin composite
bash examples/anl/orbital_reviews/molpro_confirmation/check_caspt2_shift_on_blodgett.sh
```

```bash
watch -n 10 '/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py monitor /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_shift01 --all --roots'
```

After completion, export and transfer it to `examples/anl/orbital_reviews`:

```bash
/home/hram/KinBot/.venv/bin/python /home/hram/KinBot/examples/anl/orbital_reviews/molpro_confirmation/run_confirmation.py collect /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_shift01 --output /home/hram/kinbot-runs/minimum_active_space_molpro_confirmation_20261016_shift01.evidence.tar.gz
```

Acceptance examines native completion, root energies/reference character,
first-order norms, selected orbital subspaces and changes in the spin gaps
at each basis and in their CBS estimate. A one-geometry check does not validate
shift independence over a whole PES. Production builds its own shift-0.1
asymptotes and correction scans consistently; it must not mix them with old
shift-0.3 energies. Resources: two jobs maximum, eight ranks each, 1000 MW/rank,
104000 MB/node; existing eleven exclusions inherited.

## Production work still in progress

Install the accepted basis-specific state/orbital models; use CASPT2 alone in
the correction; implement and document the chosen aDZ/aTZ CBS estimate;
audit spin-orbit/degeneracy handling across thermochemistry and capture;
update general orbital reconnaissance and the full maintained checklist.
Prepare two independent ten-job production runs, retaining previous evidence
and the optional live MRCI comparison. The matched check and final production
input audit are not yet complete native production calculations.
