# Audit of the VRC/ANL handoff, 25 October 2026

The four new orbital-screen recipes handed off at `21ec57ca` are withdrawn.
The user confirmed that the launch block was **not run**. No associated jobs
need cancellation. This audit does not authorize a replacement QC campaign,
change an accepted production CAS, prune a network, or restart controllers.

The earlier chemical analysis was available in this repository and conversation.
The assistant failed to apply it when choosing the new diagnostics: it treated
broad orbital coverage as sufficient justification for a valence CI and checked
the operational machinery more thoroughly than the scientific choice. That was
an error in selecting the work, not evidence that full valence is required by VRC.

## Findings and corrections

| Finding | Evidence and consequence | Disposition |
| --- | --- | --- |
| The four new (19,12)/(19,15) spaces were not justified by channel-specific deficiencies | The recipes include all heavy-atom valence targets, or all valence targets including H. No compact-space comparison or orbital-character evidence for these requests established that these groups were needed. | Batch and all four recipes marked `withdrawn`; preparation rejects them. Original numbers remain solely as an audit record. |
| New low-level reconnaissance could silently infer a valence CI | `molpro_reconnaissance._load` called `full_valence_active_space` when `probe_active_space=true` and `probe_definition` was absent. A hydrogen anchor could also add light-atom targets. | New preparation and planning now require a counted definition or explicit `probe_space.mode=full_valence`. Legacy evidence remains readable and already prepared jobs remain resumable. No molecular special case or arbitrary CAS-size ceiling was introduced. |
| Bond-only identity is insufficient for two new channels | The same parent has `hom_sci_2_3` and `R_Addition_CSm_R_5_1_2`, both with anchors [2,3], but different products and geometries. `validate_surface_model`, `extend_surface_model`, `channel_for_reaction`, `_selection`, batch skipping and review acknowledgement currently key on the bond. | **Unresolved publication blocker.** Separate screen directories do not fix the production contract. Reaction-specific identity must be carried through selection, publication, acknowledgement, dispatch and reuse before publishing these channels. Do not install one bond model and let it acknowledge both requests. |
| Historical documents expose conflicting recommendations | Older handoff text retains an auxiliary quartet/CC spin-gap recommendation and single-TZ frequency-shift formula, superseded by direct CASPT2 and the later TZ/QZ extension. | Added prominent precedence notices. Current chemical selections: `r7_progress_20261024/ALL_ACTIVE_SPACES.md`. Current thermal/live-control behavior: `docs/anl_live_controls.md`. |
| Operational tests were being used too broadly as reassurance | Deck/state/scheduler tests can pass for an unreasonable scientific recipe; the original test explicitly asserted the large CAS sizes. | Geometry/admission regressions now use SCF-only temporary fixtures. Withdrawal and implicit-selection refusal have dedicated tests. Test success is not chemical approval. |

The inherited-model defect is a separate, real issue: accepted case models
survived the refresh, but the worker inputs did not inherit them across the
controller generation. The `21ec57ca` fix and `--reconcile-models` operation
are retained. Their local tests pass; the supplied evidence does **not** yet
verify their application on Blodgett. A case-level “applied” receipt alone is
insufficient—worker options, review acknowledgements and actual scan graphs
must agree.

## Scientific selection audit

The objective remains the minimum chemically adequate subspace for the
specified dissociation/association potential and electronic states. Include
reacting radical/bond orbitals and the complete fragment manifolds needed for
orientation sampling, degeneracy and relevant reorganization. Include an
additional donor, antibond, sigma or pi pair for a stated chemical reason and
check its effect where relevant. Closed bonds do not all have to enter the CAS
for a correlated MR treatment to account for dynamic correlation.

The [Molpro AVAS manual](https://www.molpro.net/manual/doku.php?id=automated_construction_of_atomic_valence_active_spaces)
allows specified atomic targets and separate orbital generation; the target
choice remains a chemical input. Its suggestion to examine CASCI occupations
does not establish a requirement to use a whole-valence CI. The
[AVAS paper](https://arxiv.org/abs/1701.07862) describes construction from selected
atomic valence targets, not a universal VRC active space. The detailed earlier
VRC rationale and literature links remain in `ORBITAL_REVIEW_WORKFLOW.md`.

The recorded compact selections were re-read; none has been silently expanded:

| Existing/previously published channel | CAS and source indices | Evidence boundary |
| --- | --- | --- |
| Initial CH3S → CH3 + S | (5,4), 11–14 | S 3p manifold plus methyl radical; retained native backend and sampled angular comparisons. |
| CH2SH → CH2 + SH | (5,4), 11–14 | CH2 frontier pair plus SH pi pair; extra S–H pair near 1.98/0.02 in the expansion. No complete two-basis angular validation. |
| CH + H2S | (5,4), 8,12,13,14 | CH frontier/sigma subspace plus sulfur donor, selected by native character. Donor index changes with geometry. Quantitative convergence of the omitted S–H antibond is not established. |
| C + SH3 | (3,4), 12–15 | Previously selected C 2p set and SH3 frontier orbital; reaction subsequently rejected after connectivity change. **Still excluded.** |
| Initial CH3SS C–S → CH3 + S2 | (7,5), packed 18–22 | S2 pi/pi-star span plus methyl radical; source swap was already performed. Do not repeat it. |
| Initial CH3SS S–S → CH3S + S | (7,5), distinct 18–22 seed | Bound-fragment frontier pair plus complete atomic S 3p set. Distinct channel-specific records. |
| CH2SHS C–S → CH2 + HSS | (5,4), 19–22 | Fragment frontier pairs; provisional geometry supplied the seed, with request-binding checks. No completed relaxed-potential validation claimed. |
| CH2SHS S–S → CH2SH + S | (5,4), distinct 19–22 seed | Radical frontier orbital plus atomic S 3p set; same provisional-source limitation. |
| CH2(SH)S homolysis → CH2S + SH | (7,6), 18–23 | Extra sigma pair changed the measured inner relative potential by about 1.1 kcal/mol at 2 Å in both bases. It was an explicit expansion motivated by that comparison, not by lower absolute energy alone. |
| HCS2 + H2 product exit | (7,6), 15,19,20,21,22,25 | Heavy-fragment frontier/pi spans plus H2 sigma/sigma-star. H2-pair inclusion was a chemical modeling choice; an energetic need over (5,4) was **not** demonstrated. |
| CH2S + SH product exit | (5,4), 19–22 | CH2S pi/pi-star and SH frontier pair. Separate approach/seed; not validated by the homolysis (7,6) comparison. |

These are checkpoint-local indices from the documented basis-specific seeds.
Neither their ordinals nor their CAS dimensions transfer automatically to the
new parent. The full source paths, weights, qualitative assignments and
limitations are in `ALL_ACTIVE_SPACES.md`. The last five rows were published to
the case input; worker inheritance was the later unresolved remote issue.

Additional state-averaged roots remain distinct from reactive flux components.
Direct CASPT2-F12 remains the installed correction model; this audit does not
restore the superseded CC spin anchor. No completed scan or successful solver
is by itself proof of orbital/state continuity or a converged rate.

For the four new requests, the correct questions are:

* CH + H2SS: which CH frontier components and sulfur donor/reorganization
  partners participate? The previous CH + H2S model informs that question but
  does not supply H2SS orbital indices or prove the same CAS adequate.
* CHSH2 radical + S: preserve the atomic S 3p manifold and identify the actual
  radical's frontier span, including required delocalization partners.
* C(SH2) + SH exit: preserve SH's frontier pair and identify the relevant
  singlet fragment frontier/reorganization orbitals at the product geometry.
* HCS + H2S exit: identify the HCS radical/delocalized frontier subspace and
  whether an H2S donor/acceptor partner is needed along this noncovalent approach.

The supplied geometries establish atom maps and physical events. They do not
provide basis-specific electronic orbital assignments. No new native orbital
output for these four screens was found in the inspected local file inventory.
Before another launch, reuse applicable existing orbital evidence or obtain a
small explicit orbital-character panel, then state the compact candidates,
excluded groups, electronic states and cost. Do not replace one unsupported
CAS with another guessed index list.

## Other directives checked

| Area | Code/test finding | Limit of the conclusion |
| --- | --- | --- |
| vdW exits | `product_complex_geometry` checks the accepted network, normal product-side L2 record, atom order, fragment maps and states. The new exits use supplied rows 843/1162 and preserve internal fragment distances during translation. | No new production scan has run for these requests. Their physical event must remain distinct from homolysis, including the same-anchor collision above. |
| Scans, correction and RotdPy | Product exits start from the product complex. Correction scans use mapped interfragment atom distances; COM sampling surfaces are separate. The RotdPy correction uses the recorded scan-anchor distance in each configuration. | Neither a chosen grid nor orbital selection proves interpolation coverage, a valid asymptote or flux convergence. Native scan review remains necessary. |
| MESS inner/complex/outer path | Integration tests keep the inner TS, bound complex and required outer flux. Complex soft modes remain in the well model; the outer barrier uses fragment internal modes with sampled relative motion. Missing required flux does not silently become phase-space theory. | Harmonic treatment of a floppy complex is still an approximation. This is not a completed MESS solver run or mode-by-mode scientific approval. |
| ANL geometry/frequencies | MESS takes retained accepted high-level geometry/harmonics. Individual CC frequencies are not CBS-extrapolated. TZ and QZ VPT2 shifts are separately projected onto each retained mode, with ambiguity/quality/provenance checks. Scalar harmonic and anharmonic ZPVE arithmetic remain separate. | The 3.7 power for the requested double-hybrid anharmonic extension is a documented modeling convention, not separately benchmarked for these species merely because ANL electronic extrapolation uses it. |
| Rotors | MESS projects accepted L2 rotor identities from the retained Hessian; rotor blocks keep their L2 geometry/potential. Demotion restores the vibration and recomputes shifts. CBH 0 K energy retains full harmonic ZPVE and scalar VPT2. | A quantum hindered-rotor ground-state replacement and torsion/vibration coupling were not implemented. The retained vibrational 0 K treatment is explicit, not an exact treatment of floppy torsions. |
| Complex binding and CBH | Complex formation uses fragment formation energies plus a common harmonic complex-minus-fragments binding difference. There is no covalent CBH fragmentation of the disconnected complex. Individual products retain their independent VPT2 thermochemistry. | Retaining a VPT2 job does not certify its floppy-mode result; complex VPT2 is not automatically inserted into binding/MESS. |
| Live task choices | Tests cover one-task omission, retained-TZ geometry with higher-level single points, tier choices, dependency invalidation, preserved unrelated/reference jobs, generation fencing and terminal cancellation. | Reference arithmetic must have a compatible common expression; references are not cancelled just because one target omits a contribution. Already-running old code must first have the relevant feature loaded. |
| Cross-PES sharing | Opt-in, directional, provenance-bound whole-graph aliasing before receiver execution; canonical owner and history remain intact. | This is **not arbitrary component-level reuse** between differing graphs or geometries. Exact scientific compatibility is required. No sharing route or simplification was activated by this audit. |
| Native retries | CC limit retries are bounded (Molpro 100/200; eligible CFOUR VCC 600/1000); unsupported MRCC compiled-limit cases stay visible. `hsdel` recovery varies solver P-space/workspace with two replacements, preserving CAS, roots, shifts and Hamiltonian. Sampling retries the same geometry and persists exhaustion. | These are failure-class-specific procedures, not “retry every QC error.” Native success remains to be verified from actual outcomes. |
| Resources and preservation | Shared review admission counts within 23, preserves exclusions/policy bytes and existing allocations, fences uncertain submissions and prevents duplicate owners. Cancellation tests wait through `COMPLETING`. | No live scheduler query was made in this audit. A historical monitor is not the current queue. |
| Generality | Reviewed installed modules contain no new parent/run identifiers or molecule-dependent recovery branches. Chemistry and request IDs remain in example data/tests. | This was a scoped code audit, not proof of correctness of every KinBot module. |

The latest supplied 08:47 PDT monitor had CH3S exploration complete but CH3SS
still exploring (five active exploration workers and one unstarted well).
Final kinetic pruning, orphan-product removal and shared CBH dependency review
therefore remain deferred as requested. No timescale prediction can be made
from that snapshot alone while new wells are still entering the network.

## Evidence and validation

The locally retained evidence export has an earlier verification receipt for
8,432 files, including 319 WFU files, with archive SHA-256
`9256c3a78c940582d9f0788d03a0ce99c955d50164c0920eb774bbeb6d96daf5`.
This audit read that receipt; it did not repeat the entire archive hash pass.
File integrity is not scientific validation. No native files, accepted
geometries or production results were modified here.

Targeted local regression runs: **801 passed, 1 skipped**, with 29 additional
subtests passing. Of those, 241 cover reconnaissance preparation, saved evidence,
product geometry and review transport; 444 cover model inheritance, retained
refresh, decisions, resource admission, cancellation, retry classes, sharing,
ANL arithmetic and MESS/vdW behavior. Another 116 cover campaign/confirmation
compatibility and the earlier follow-up/queue orchestration. The new transport fixtures are SCF-only;
the deliberate historical valence fixtures explicitly opt in for replay.
The final reconnaissance rerun needed local process inspection outside the
sandbox for two ownership tests; the complete 241-test rerun then passed.

These tests use native-output fixtures and simulated execution/scheduler paths.
They did not run Molpro, Gaussian, Slurm or production MESS. The repository-wide
suite was not rerun; historical whole-suite failures described in
`anl_live_controls.md` are not declared resolved by these targeted runs.

**Release state:** the prior new-screen launch is withdrawn. Existing selected
models remain unchanged. Replacement compact selections and reaction-specific
channel routing remain prerequisites to a new publication handoff. The accepted
model inheritance repair has local evidence but still needs remote worker-level
confirmation. No controller restart or new QC launch is part of this audit.
