# Paused at user request, 25 October 2026

Historical checkpoint. The user resumed on 26 October; current work and remaining
dependencies are recorded in [the resumed handoff](../r7_root_spectrum_20261026/README.md).

No further calculations, publication, cancellation or controller refresh were
performed. The user owns the Blodgett SSH session. Resume only when requested.

## Latest verified evidence

Local export: `examples/anl/tools/natural-union-20261026T012004Z/`.
Archive SHA-256:
`8621cbfed310ed09d149cba2d0b6c71cd2b8f6a0478dc5b8afb4d1208c7102fa`.
All 731 archived files verified against the manifest; zero collection issues.
592 text/XML files extracted safely to `/tmp/kinbot-natural-union-audit`.
This archive includes the preceding ten-job union panel and the four new jobs.

All four new CAS(9,6), two-doublet calculations completed. The explicit near
natural checkpoint was read correctly (`2140.2`, MCSCF/NATURAL, state averaged;
16 closed, 6 active, 9 active electrons). Nevertheless, the optimized near
space again excludes the intended SH-p direction. Minimum overlap singular
values against the intended old near natural 17–22 union are 0.121846 (DZ) and
0.122409 (TZ). New active 17 is principally old occupied 12–14; the SH-p
direction remains inactive 16. Near energies reproduce the preceding canonical
initialization result to about 2e-8 Eh. Another identical seed retry is exhausted.

The new near-to-far continuation also differs from the independently optimized
far branch. Far D1/D2 energies are -834.481898752232/-834.470943757414 Eh (DZ)
and -834.496164694357/-834.486039847537 Eh (TZ). Its leading D1 CI strings
`222aba`, `222baa`, `222aab` describe three singly occupied active orbitals;
D2 is dominated by `2222a0`. The SH pair is no longer represented by two
half-filled averaged directions. This suggests competing fragment-spin/state
character, but does not establish the isolated fragment ground-state ordering
or justify changing production roots/weights. Scientific acceptance is pending.

## Work remaining

1. Resolve this last C(SH2)+SH exit using the native root/CI evidence and both
   optimized branches. A bounded fixed-orbital CASCI root-spectrum diagnostic
   in the already identified six-orbital union is under consideration, including
   doublet/quartet comparison to distinguish singlet-fragment and triplet-fragment
   channels. No such recipe or code change has been implemented, approved as a
   production model, or launched. Do not automatically enlarge the CAS or repeat
   the completed initialization trials.
2. Finish the numerical audit and hash bindings for the preceding ten-job panel
   using the newly verified archive. Eleven new selections are prepared; the
   twelfth remains unresolved. The combined publication plan remains held.
3. Complete the requested single combined publication of all twelve new CH3SS
   channels and reconciliation of the five previously published models after
   scientific selection is complete. Do not silently release a partial set.
4. Finish the full handoff: existing and new active spaces, qualitative orbital
   characters and source versus packed indices, roots/weights versus reactive
   flux, evidence and limitations, export instructions, fixes and refresh commands.
5. Refresh/restart owning controllers to load the endpoint and publication fixes,
   preserving live Slurm jobs, limits, exclusions, operator decisions and RotdPy
   campaigns. Complete orchestration tests before that handoff. Verify remote
   worker adoption, actual scan dispatch and sampling adoption afterward; an
   applied receipt alone is insufficient. The audited endpoint bug indicates
   VRC-preparation retries, not earlier PES/IRC/L2/ANL/CBH reruns.
6. Reconcile the latest production monitor in the handoff. At 18:20 PDT: cap 23,
   19 allocations (CH3S 13, CH3SS 6), zero diagnostic allocations, policy revision
   14 and 16 applied live controls. CH3S remains in postprocessing; CH3SS still
   has five wells exploring and five in VRC. Final network pruning is premature.
   Two suprafacial migration searches logged rejected coordinate modifications
   at 18:09:40; inspect attempt evidence before deciding whether recovery is
   needed. Historical inactive failures must not be treated as fresh failures.

## Repository state at pause

Latest pushed commit is `a70d48fa` on `composite`; endpoint fix is `a5cb22ed`.
No installed-code edits were made during this review of the new export.
At the pause, this note was local and uncommitted. Preserve the user's unrelated edit
to `r7_progress_20261024/EXPORT.md` and all untracked native archives/local runs.
Previously reported tests are recorded in README.md; they are local simulated
validation, not evidence that a native retry or remote controller handoff worked.
