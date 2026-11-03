# Homolysis coverage and reverse abstraction discovery

`hom_sci` deletes a bond to propose fragments. That operation does not locate a
transition state and does not establish a barrierless dissociation. Ordinary
unimolecular reaction discovery also does not automatically launch bimolecular
abstraction searches from each separated-product pair.

## Forward checks

`homolytic_saddle_search` defaults to `true`. KinBot attaches a separate
`dissociation_saddle` search to every selected `hom_sci` proposal. It waits for
ordinary discovery and accepted-result reconstruction to finish. Existing
scientific rejections, including above-ceiling fragment asymptotes, are retained.
They are recorded as `parent_not_retained`, not retried as QC failures.

A normally validated saddle can cover a proposal only if its IRC breaks the
particular proposed bond without other connectivity edits and reaches the same
configured fragment states. A rearrangement reaching the same fragments does
not suffice. Otherwise the companion searches the actual L1 surface:

- Relax all coordinates except the cleavage distance, starting at the optimized
  parent bond length; default 40 points separated by 0.1 angstrom.
- Require a converged finite constrained minimum at each point. A failed guess
  is never a normal zero-energy observation.
- Use a sampled maximum as a TS guess, then run the normal independent saddle
  frequency, IRC, endpoint, energy-screen, and L2 refinement procedures.
- If no maximum is found, require the final fragment identities to match, the
  electronic energy to approach the separately computed fragment sum within
  0.5 kcal/mol, and the last three points to span at most 0.1 kcal/mol before
  labeling PST **provisional finite-path scouting evidence**. This is not a
  proof that no saddle or alternative electronic-state path exists.
- Failed, nonconverged, mismatched, or nonasymptotic checks remain unresolved;
  they cannot silently authorize a new combined MESS export.

The parameters are `homolytic_saddle_scan_points`, `homolytic_saddle_scan_step`,
`homolytic_saddle_asymptote_tolerance`, and `homolytic_saddle_flatness`.
`homolytic_saddle_search` also accepts `false` or a list of exact `hom_sci` names.
The durable evidence is in `*.dissociation_scan.json`, `homolysis_checks.json`,
native QC artifacts, and `.kinbot_search_state.json`. Reentry adopts live work,
retains a pending maximum across pool-full interruptions, and does not repeat
an exhausted search.

## Reverse searches

`python -m kinbot.reverse_search CASE` inventories the retained fragment pairs
without changing them. `--prepare` stages independent `bimol=1` workers below
`CASE/reverse_abstractions`. The actual bimolecular finder enumerates abstraction
motifs in both directions. Pairs are deduplicated by configured fragment identity
and total multiplicity; electron spin coupling and charge conservation are checked.
An explicit `l1_pes_energy_ceiling` and recorded L1 profiles are required. The
same initial-PES reference is used, rather than treating each pair as a new zero.

Staging copies the selected fragment observations and their native artifacts,
with source hashes and no old scheduler ownership. It does not create native
completion evidence absent from the source. Source changes during copying stop
staging. Interrupted preparations can be resumed without overwriting work.

The workers use inherited L1/L2 settings and shared ancestor resource controls.
ANL, ROTD, and their MESS exports are disabled. The ordinary refresh owner
launches/adopts them and includes them in subsequent controller replacements.
The Sella abstraction path now uses a saddle optimization, not a minimum search,
on its transfer guess; Gaussian's QST3 route remains available. Alignment works
with either fragment first in the atom array.

Accepted reverse reactions retain both IRC endpoints. A matched L1 Hessian and
fragment E+ZPE establish any retained entrance-complex depth; a retained complex
receives L2 characterization. `reverse_search_result.json` includes the native
TS, product, entrance, and path-identity evidence. Failed entrance checks are
unresolved, not silently promoted to a MESS well.

**Integration boundary:** this change stages, runs, and exports reverse
discovery. It does not yet implement automatic reverse-result integration into
the combined PES. Newly found saddles must be matched to actual entrance and
product wells using their returned native evidence. It is wrong to connect two
bimolecular product labels by an ordinary unimolecular MESS barrier or assume
that a paper's vdW complex is the calculated entrance. Until integration is
completed, the combined-MESS exporter stops explicitly; the previous MESS file
is stale for this discovery revision. This is not a claim that a native QC
calculation failed. `--status` lists controller progress and result locations.

## Controller restart

After staging both cases, use the existing owner:

```bash
python -B -m kinbot.run_refresh RUN --check-all-homolysis --apply --wait-idle 120
```

This loads the new code and covers every retained homolysis, preserving accepted
fragment/TS calculations and unrelated rejected searches. It does not mean
rerunning all fragment calculations. `--check-homolysis EXACT_NAME` provides
optional narrower scope. Do not combine this with a broad rejected-search retry
unless that separate change is intended.

## Validation and upstream review (2026-11-02)

- Fetched `origin/master` at `684e513c4`. It already contains the separate
  bidirectional abstraction finder. Its unrelated backend-removal refactor does
  not supply this forward/reverse orchestration; no wholesale merge was made.
- The captured scout has 9 retained CH3S and 51 retained CH3SS homolyses. The
  reverse planner finds 7 motifs in 5 active pairs and 26 motifs in 15 active
  pairs, respectively. These are proposed motifs, not promised accepted TSs.
- Tests exercise real finder/IRC/stereopath code, generated constrained and
  saddle workers, native-state checks, live-job and pending-submission adoption,
  crash recovery during staging, input-decision application, sampling retention,
  and repeated refresh/preparation.
- Captured databases were staged into a disposable local copy and all reverse
  inputs passed the real L1/L2 router startup. The capture lacks some native log
  completion stamps, so it cannot prove that remote fragment reuse or new native
  Gaussian searches will succeed. No chemistry was run in this local replay.

No finite test suite proves that all bugs or all chemical pathways are absent.
Native convergence and connectivity remain calculation outcomes to inspect.

## Discovery holds and scan-target repair (2026-11-03)

A saddle proposal that does not produce a validated saddle is an expected
discovery outcome. It is not automatically a software fault or a useful retry.
An unresolved homolysis now publishes an **export hold** and releases its
exploration slot. The other wells continue. The coordinator keeps the network
open for targeted, evidence-supported revisions and withholds new rate export.
Standalone and combined MESS guards remain in place. A failed search does not
prove either that a tight saddle exists or that the dissociation is barrierless.

The latest remote evidence identified a distinct implementation error: an
approximate coordinate modification was used as the frozen Sella coordinate,
even though the scan validator checked the requested distance. Quantitative
dissociation scans now pass the explicit coordinate target into Sella. Force,
frequency, IRC, electronic-state and scan-distance acceptance criteria remain
unchanged.

On owner reconstruction, a `scan_constraint_not_met` rejection is eligible for
one targeted repair only when the saved legacy worker requested a strict scan,
froze only the cleavage bond at its approximate input distance, and the normal
native result agrees with that frozen distance but not the requested target.
The current target and recorded search step must agree as well. The owner
archives the old worker/native artifacts with hashes under
`.kinbot_scan_repairs`, records the revision, and repeats only that point. All
previous completed scan points and database observations remain. Live or
unknown QC, failed QC, already-explicit targets, different constraints and
exhausted repairs are ineligible. Journal-before-checkpoint interruptions,
pool admission, and repeated refreshes preserve the same single repair intent.

Reverse discovery also preserves the two native IRC job identities before
single-fragment extraction can rename a product or L2 optimization can mutate
its coordinates. Export reads those L1 observations and characterizes a fresh
entrance object, so repeated exports cannot substitute product L2 geometry for
the original IRC endpoint. This repairs the observed species-name/L1-route
error; it does not waive provenance or approve a new reaction.

Validation includes generated-worker execution with the actual local Sella
optimizer and an analytic test potential, both IRC direction assignments,
single-fragment object mutation, repeated export, interrupted targeted repair,
live-job adoption, exploration-slot release, export blocking, and the complete
simulated-site refresh with sampling and decision preservation. These checks
do not claim a successful native Gaussian retry on Blodgett. The disappearance
of controllers whose saved status still says `running` is not explained by the
provided records; refresh can reconstruct their owned progress, but an external
termination cause has not been established.

## Routed mode-reading correction (2026-11-02)

The first remote restart exposed an interface error during accepted-TS
reconstruction: `ProfiledQuantumChemistry.read_imag_mode()` did not accept the
`geometry` and `symbols` arguments required by `Optimize.compare_structures()`.
The backend reader supported them, but the router had not been updated. This
was a controller exception, not a native QC failure or a rejected pathway.

The router now forwards the complete reader arguments to the recorded L1/L2
backend, including the selected Cartesian frame. No mode-fit tolerance,
frequency criterion, electronic-state definition, or scientific result changes.
The tests now join the actual router, native mode reader, L1/L2 acceptance,
database adoption and repeated restart; they also exercise a mixed ASE/Gaussian
profile, incomplete/mismatched evidence, and compatibility of all shared public
router method signatures. The simulated site handoff includes live reverse
allocations and sampling work. A captured spectrum reused under two fixture job
labels tests coordinate rotation; it is not two independently computed energies.

For a run already prepared with homolysis and reverse discovery, load this
correction with an ordinary controller refresh:

```bash
python -B -m kinbot.run_refresh RUN --apply --wait-idle 120
```

Preparation and homolysis re-selection are not needed for this correction.
Existing Slurm calculations and scientific search decisions are retained.
Local tests do not certify successful execution of the subsequent remote search.
