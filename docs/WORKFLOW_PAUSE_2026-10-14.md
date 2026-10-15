# Workflow pause — 14 October 2026

**Historical checkpoint, superseded:** the user immediately resumed work.
Reconnaissance has since been released; see
[the current runbook](Blodgett_workflow_runbook.md) and
[the ongoing audit](ANL_workflow_audit.md). The text below records the pause
state and is not an instruction to stop current work.

Paused at the user's explicit request. Resume only when the user asks. No
automatic half-hour restart or scheduled action has been created. There are
no cluster jobs or SSH sessions started by this work, and no local test or
download sessions remain running. The user performs all SSH.

## Objective and authoritative constraints

Finish the CBH/ANL → VRC-TST → ROTD_py → MESS audit and implementation, with
regression/compatibility checks against master, and prepare CH3S → CH3 + S,
CH3SS → CH3 + S2, and CH3SS → CH3S + S as separate channels. Then hand off
complete Blodgett cleanup, source update, and orbital-reconnaissance commands.
The user will supply the recon results before the full KinBot campaigns.

- Use **Molpro's default F12 bases**. The earlier rev2 instruction was
  superseded. No forced revision, vendored basis definitions, generated
  diffuse exponents, or auxiliary substitutions.
- General input-driven spins, term degeneracies, roots, active spaces,
  coordinates, and statistical weights. Molecule-specific physical data
  belong in examples, not implementation branches.
- At most **ten nodes total**, including both systems. Run the two recon
  controllers sequentially or enforce a genuinely shared cap.
- Recon must provide radial, state, natural-orbital and orientation evidence,
  adaptive refinement and asymptote checks. It must not certify a minimal
  active space from a few fixed distances or a population cutoff alone.
- Verify the multistate approach against the **actual NH2O supporting MESS
  files**, not only the paper's text/figure.
- Preserve the existing vdW-well detection/inclusion policy; electronic
  state accounting is a separate concern. A complex must not become
  mandatory merely because the NH2O SI retains one.

## vdW clarification checked at the pause

Read-only inspection of `origin/master` confirms:

- `vdW_detection` defaults to **0.5 kcal/mol** and recognizes deeper
  IRC-product complexes. The comparison is `vdW_depth > vdW_detection`.
- Master's reaction-generator check explicitly requires two products and
  **excludes `hom_sci` channels**. Therefore that existing check cannot be
  assumed to characterize the sulfur homolytic VRC entrance complexes.
- `me_skip_vdW` defaults to zero. Setting it to one omits optional complexes
  and their exits from the direct kinetic network while retaining the selected
  complex reference for inner-barrier tunneling. PES assembly has matching
  retained/omitted-complex handling.

Relevant current code: `kinbot/parameters.py` (`vdW_detection`, `me_skip_vdW`),
`kinbot/reaction_generator.py` (depth check), `kinbot/mess.py`
(`_include_complex`, `_product_complex`), and `kinbot/pes.py`
(`create_mess_input`). Do not replace those decisions with state counting.
For homolytic VRC channels, review how significant long-range minima are
detected/characterized and attach that evidence to the existing well machinery.

## Last fully verified state

Before the implementation work in the current turn, the full local suite
passed **1198 tests and 389 subtests**, with the existing `imp` deprecation and
Sella cache-permission warnings. All 156 then-existing non-template source
modules parsed; `git diff --check` passed.

The additional unmodified `origin/master` suite gave 602 passes, 41 failures
and 369 passing subtests. The failures were classified as older fixture
assumptions (missing real MESS inputs, missing scheduler-generation completion
markers, removed private executable mocks); see `ANL_workflow_audit.md`.
Do not represent that older suite as passing or universal compatibility as
proved.

The starting workspace already contained extensive uncommitted work. Preserve
it. Branch is `composite`; HEAD was b3f790a. `.git` is read-only. The untracked
`.kinbot_run.json` predated the audit and must not be deleted or deployed.
Use `.venv/bin/python` (3.11); system Python lacks the test environment.
Do not compile raw `kinbot/tpl/*.py` templates as ordinary Python modules.
There are no local licensed Molpro, Slurm or MESS executables. Native chemistry
execution has not been validated here.

## Implementation added this turn — unfinished and not release certified

`kinbot/vrc_multistate.py` is new. It currently provides:

- A reviewed external surface-model schema selected by one-based bond indices.
  States carry multiplicity, selected root, population weight, explicit
  energy offset/reference note, and optional resolved electronic levels.
- Independent native relaxed paths and correction records for scalar state
  leaves, with separate job names/cache/provenance; shallow copies retain the
  original parent electronic state. Highest-spin leaves use direct trusted
  energies rather than a lower-spin CC anchor.
- Recursive generation of scalar ROTD inputs, with reactive-pivot or COM
  dividing surfaces and stage-specific theory/orbital overrides.
- State-table alignment, bounded interpolation, physical flux summation,
  explicit-level or multiplicity counting once, and hash revalidation.
  Correction variants are never additional electronic states.
- Single-bottleneck output; optional reduced statistical/minimum two-TS MESS
  models; and a newly added explicit complex + inner/outer barrier renderer.

Integration edits are in `kinbot/vrc_tst_scan.py`, `kinbot/pes.py`,
`kinbot/mess.py`, `kinbot/vrc_correction.py`, and `kinbot/parameters.py`.
New parameters include `vrc_tst_surface_model_file`,
`vrc_tst_selected_roots`, and `vrc_tst_state_average_roots`.

**Do not treat the explicit-network path as complete.** It currently accepts
a reviewed complex RRHO template/depth. That path must be reconciled with
the inherited depth/skip policy and existing characterized complexes, rather
than forcing a new manually specified well into every multistate channel.
The final connection into network normalization, stereochemical population
metadata, CBH energy references, duplicate-complex handling and PES assembly
still needs review. `_rotd_model` currently defers PES workers, so full-PES
multistate result staging must be checked/implemented explicitly; avoid a
silent phase-space fallback when state-resolved output was requested.

`kinbot/molpro_reconnaissance.py` now also has optional rigid-fragment
orientation probes (`distance_scales`, `angles_degrees`, named bond/transverse
axes, optional asymptote). They preserve internal geometry/bond separation,
skip atomic/unchanged rotations, use independent AVAS starts and report
orientation profiles separately from radial refinement. The examples have
**not yet been updated to enable these probes**. Recompute counts after doing
so; the older 152/312 job counts apply only to the radial configuration.

The last code edit added direct trusted-level CBS pairs and preset F12 trusted
MR levels, to fix the high-spin leaf falling back to conventional `avtz` and
therefore demanding an absent reviewed WFU. New helpers are
`trusted_settings` and `assemble_trusted_reference`; direct correction jobs
use `trusted:0`/`trusted:1` roles. This change has **not yet been tested**.
Its extrapolation power is now a configurable parameter.

## Tests and exact latest results

- `tests/test_vrc_multistate.py` was added. Its first ten tests passed before
  later explicit-network/level-list/CBS changes.
- Earlier focused compatibility run: 55 passes and 10 subtests.
- Latest completed focused run:
  `pytest -q tests/test_vrc_multistate.py tests/test_molpro_reconnaissance.py
  tests/test_rotdpy_profiles.py tests/test_complex_mess.py`
  gave **50 passes, 1 failure, 29 subtests**. The failure was the generated
  high-spin scalar input asking for a conventional `avtz` orbital record;
  the new F12 trusted-pair change addresses its cause but has not been rerun.
- New orientation geometry/report tests were in that focused run.
- No full suite has been run after the multistate changes. `git diff --check`
  passed at the pause. A source parse passed earlier in this turn, before
  the latest edits; rerun it.

Next tests must cover the actual generated scalar templates/native root
selection, direct trusted CBS algebra, physical resolved-level summation,
explicit retained/omitted complex topology, and source invalidation. Exercise
the complete mocked calculation → correction → input → result → MESS path,
then run the full regression suite. Do not claim native MESS parsing without
actually executing MESS or clearly labelling a structural check.

## NH2O SI evidence downloaded and inspected

Paper: Klippenstein, Mulvihill, Glarborg (2023),
DOI https://doi.org/10.1021/acs.jpca.3c05181.

All **ten TXT SI files** (002–011) were downloaded through the publisher's
Figshare API and checked against supplied MD5 checksums. Locations:

- `/private/tmp/nh2o_si/jp3c05181_si_002.txt`: complete MESS network.
- `/private/tmp/nh2o_si/jp3c05181_si_003.txt` through `_011.txt`: numeric Ne
  tables. Do not guess the original named-table mapping from their numbering.
- `/private/tmp/nh2o_si/sources.json`: public URLs, DOIs, SHA-256 and licence.
- `/private/tmp/nh2o_si/24274749.json` etc.: publisher metadata.
- `/private/tmp/fetch_nh2o_si.py`: completed download/checksum utility.
- Author manuscript text/PDF: `/private/tmp/nh2o_klippenstein.txt` and `.pdf`.

Metadata endpoint:
`https://api.figshare.com/v2/articles?resource_doi=10.1021%2Facs.jpca.3c05181`.
Network SI DOI: https://doi.org/10.1021/acs.jpca.3c05181.s002.
Network SI SHA-256:
`f76e0cdb130618f55850ade3e61f8268ea38e95c1ce877fb209f3d2b9e648403`.
The SI licence is CC BY-NC 4.0; do not silently relicense full SI files as BSD
test fixtures. Store provenance/derived contract facts and use an external
SI directory for reproducible verification if needed.

Direct inspection (network lines roughly 338–455) shows:

- NH2 + O: outer B1v links P1 to W1v; inner B1 links W1v to W1. Both use
  unions of state-specific RRHO/ROTD contributions (three outer, two inner).
- NH + OH: outer B2v links P3 to W2v; inner B2 links W2v to W2. Two outer
  RRHO contributions, one inner contribution.
- Outer NH2 + O electronic levels across its three components sum to
  degeneracies 10 at 0, 6 at 158.5, 2 at 226.5 cm^-1: consistent with NH2
  doublet times the O J-level degeneracies 5/3/1. Outer NH + OH has two
  sixfold contributions with distinct zero energies.
- Inner tables are assigned unit electronic degeneracies in the SI. Their
  table normalization must be respected; do not add another spin factor
  without evidence. The new explicit-level override supports that convention.
- The SI uses retained complex wells and competing pathways, not merely a
  `Variational/2TSMethod statistical` block for those entrance channels.
  The reduced two-TS option must be described as a reduction, not identical
  reproduction of that network.
- A first outer block spells its core `Roti`; another complex frequency
  declaration appears inconsistent with its listed count. Do not blindly
  reproduce apparent source inconsistencies as generated syntax.

Still required: write a reproducible SI structural/level-accounting verifier,
run it on the downloaded files, save its provenance/result, and add meaningful
regressions against those derived facts. Native rate reproduction is a
separate validation and is not yet available locally.

## Remaining work after resume

1. Reconcile multistate accounting with inherited vdW detection/skip behavior.
   Preserve separate state accounting for retained and omitted complexes.
   Address homolytic VRC-path minima explicitly rather than assuming the
   IRC-only master check covered them. Do not force extra wells universally.
2. Finish SI-backed validation and the network/PES/provenance integration
   above. Check energy-zero alignment, long-/short-range SO treatment,
   electronic normalizations and common conserved vibrational assumptions.
3. Run focused tests, fix any failures, then run the complete suite and
   non-template source parse/diff checks. Review backwards compatibility,
   including scalar inputs and archived outputs.
4. Rebuild both sulfur examples: enable useful orientation probes, retain all
   eight requested native basis families, independent C–S/S–S contexts, and
   candidate active-space labels. Add a reviewed-surface-model gate/template
   so future full runs cannot silently use one root. Do not invent reactive
   subsets or molecular SO/Jahn–Teller levels before recon evidence.
5. Provide a portable, reproducible source update. The changes are uncommitted
   and unpushed, so `git pull` alone cannot deliver them. Prepare an archive or
   explicit rsync workflow including new modules. Preserve the pinned ROTD_py
   Git revision f3352c8eea9d51d26f631d3238e24e304b5be7de; its runtime verifier
   checks Git metadata, so copying the submodule's files without a usable
   repository is insufficient. A Git bundle for that dependency is one option.
   Do not overwrite the user's existing checkout or environment blindly.
6. Finish the Blodgett runbook and an executable recon driver: environment
   import/preflight, scoped acetaldehyde stop/archive, fresh recon directories,
   restart/progress/result collection and sequential ten-node controllers.
   Use the existing environment; no SSH from this agent. Return commands for
   this initial diagnostic campaign without circularly requiring its results
   to certify the production chemistry first.
7. Update `ANL_workflow_audit.md`, both READMEs, and implementation docs with
   the final verified state and remaining native/chemistry acceptance checks.
   They currently describe the pre-implementation gap and older radial job
   counts and must not be treated as the final release instructions.

## Acetaldehyde diagnosis/cleanup context

Run: `/home/hram/kinbot-runs/acetaldehyde_anl0_f12_day`.
The cluster restart explains the dead campaign and saved `running` branch
labels. HIR failures explicitly selected harmonic fallbacks; the failed vdW
optimization was optional. Neither proves a fatal exception.

The earlier first `relaxed_scan` attempt's cause remains unknown without its
native logs and execution record. Exact point:
`vrctst/molpro_scan/441041030570000000001_hom_sci_1_2_vts_pt00_5914412a7cd7e609/attempts/relaxed_scan/1`.
Preserve those logs, stdout/stderr, execution JSON and Slurm accounting before
archiving. `kinbot/run_cleanup.py` provides a run-scoped owned-process/Slurm
WorkDir stop and archive utility. It must not use user-wide kills or delete
unattributable shared scratch. The final cleanup/update/recon command set is
still to be delivered after resume.
