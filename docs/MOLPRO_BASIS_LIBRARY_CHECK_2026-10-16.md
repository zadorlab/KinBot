# Blodgett Molpro basis-library check — 16 October 2026

The user ran the standalone library probe with `/opt/molpro/24/bin/molpro`,
which reported version **2024.1**. All **27** lookups reported `LOADED` for
the H, C and S atom inventory. This is the user's pasted execution summary;
the individual native output files have not yet been transferred or inspected
locally. They remain under `/home/hram/kinbot-basis-check.AiPuxM`, with scratch
under `/tmp/kinbot-basis-check.l0uNzp`.

| Output indices | Role | Successfully loaded definitions, in order |
| --- | --- | --- |
| 1–8 | Conventional ANL orbital bases | cc-pVDZ, cc-pVTZ, cc-pVQZ, aug-cc-pV5Z, aug-cc-pV6Z, cc-pCVTZ, cc-pCVQZ, aug-cc-pCVTZ-DK |
| 9–12 | VRC/ANL F12 orbital bases | avdz-f12, avtz-f12, aug-cc-pVQZ-F12, aug-cc-pV5Z-F12 |
| 13–15 | DZ auxiliary lookups | vdz-f12/jkfit, avdz-f12/mp2fit, cc-pVDZ-F12/optri |
| 16–18 | TZ auxiliary lookups | vtz-f12/jkfit, avtz-f12/mp2fit, cc-pVTZ-F12/optri |
| 19–21 | QZ auxiliary lookups | vqz-f12/jkfit, avqz-f12/mp2fit, cc-pVQZ-F12/optri |
| 22–24 | 5Z auxiliary lookups | v5z-f12/jkfit, av5z-f12/mp2fit, aug-cc-pwCV5Z/mp2fit |
| 25–27 | Diagnostic aliases/base sets | av5z-f12, cc-pV5Z-F12, cc-pV5Z-F12/optri |

Each output is `<report directory>/<index>/load.out`; the corresponding
input and launcher log are retained beside it. The check uses `GDIRECT; INT`
and no SCF or correlation commands. Auxiliary definitions are loaded as the
probe's orbital basis to check library access. Success does not establish
the identity of aliases, the F12 driver's automatic auxiliary selection,
electronic convergence, or the basis installations of Gaussian/CFOUR/MRCC.

**Decision:** retain the requested augmented QZ-F12/5Z-F12 ANL1-F12 pair.
There is no observed basis-library reason to impose a QZ maximum or modify
the production inputs. The VRC aDZ/aTZ gap extrapolation remains separate,
with exponent 3.7.

The published augmented family includes aV5Z-F12: see Table II of
[Sylvetsky, Kesharwani and Martin, JCP 147, 134106 (2017)](https://doi.org/10.1063/1.4998332).
The library-loading procedure follows the Molpro
[basis-input](https://www.molpro.net/manual/doku.php?id=basis_input) and
[integration](https://www.molpro.net/manual/doku.php?id=integration) documentation.

The [retained script](../examples/anl/caspt2_production/check_molpro_basis_library.sh)
redirects launcher stdin from `/dev/null`. The original version allowed the
launcher to consume the loop's basis-list stdin and stopped after one lookup.
A local mock-launcher regression reproduced that failure and verified all
27 lookups after the repair, including continuation past a simulated missing
basis. The user's subsequent native run reported all 27 loaded.

This check changes no scientific settings and submits no scheduled calculations.
For the later Gaussian failure and the current clean-restart command, use the
[production handoff](../CASPT2_PRODUCTION_HANDOFF.md).
