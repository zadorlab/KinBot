# R9 requirements and early R8 validation review — 14 November 2026

This records the user's latest requirements for both CH3S and CH3SS. It
supersedes the older draft's L1/L2 profile and blanket direct VRC correction
choice. It is a staging requirement, not a submission input or a completed
native validation. Existing running calculations and campaign definitions have
not been changed.

## Required production profile

| Stage | Requested definition |
|---|---|
| L1 | B2PLYP-D3(BJ)/aug-cc-pVTZ |
| L2 | B2PLYP-D3(BJ)/aug-cc-pVQZ |
| L3 | ANL1-F12 with augmented components |
| VRC | Previously reviewed channel-specific active spaces, orbitals, roots and weights; CASPT2-F12 and the applicable reviewed spin-splitting correction |
| Exploration | Only selected searches after actual connections and pruning relevance are reconciled; include CH3 + S |

Here aTZ/aQZ mean aug-cc-pVTZ/aug-cc-pVQZ. D3(BJ) preserves the current
dispersion convention. CASPT2 spin gaps and the high-spin coupled-cluster
anchor are separate components. A fragment pair with only one coupled spin
surface cannot acquire a fictitious high-spin partner: retain a direct
CASPT2 correction unless a different reference state has explicit review.
Existing spin-orbit/electronic-level declarations remain separate from this
spin-free surface correction.

The exact R9 channel list remains provisional. The current paper comparison
and L2 IRC evidence have changed several former L1 endpoint assignments.
Previously reviewed orbitals remain useful evidence; they do not establish
that an unresolved H-loss path is barrierless. A verified tight path belongs
in the tight-TS model, not automatically in VRC.

## Early review of the running validation panel

The staged panel is conservative: 220 native tasks plus up to 176 reused-saddle
IRC/endpoint jobs. These are unequal-cost tasks, not 396 new saddle searches.
It need not finish in full before assessing which unsubmitted work remains
useful. The controller currently dispatches the staged dependencies; it has
no automatic scientific stop when enough evidence is obtained.

1. Prioritize actual L2 connections for likely retained routes, disputed
   TS2/TS3/TS5 paths, the separate TS6 candidates, and the intended H-loss
   connections. Use native paths and stopping criteria as well as endpoints.
2. Once a disputed B2PLYP path is stable to independent integration, decide
   whether the remaining question is integration error or a difference between
   electronic surfaces. Reserve costly coupled-cluster paths for unresolved
   differences that affect the retained network or the paper interpretation.
3. First candidates for deferral are checks of exclusions far above 150
   kcal/mol, benchmarks of routes outside the shortlist, and additional
   confirmation of already consistent paths. Near-cutoff exclusions remain
   higher priority. An observed method shift is not a rigorous error bound.
4. Require numerical relevance under the agreed thermal and photoexcited
   scenarios before production inclusion. A valid route can still be pruned;
   a method/basis change can alter a marginal branching conclusion.

No cancellation is authorized by this note. Review a current progress/native
snapshot before changing any pending work, and preserve completed and running
calculations. MESS remains paused until the connection reconciliation is ready.
Augmentation sensitivity can be assessed later on the retained weak complexes,
associated saddles/asymptotes and marginal branches rather than expanding this
entire scouting panel now.

## Reuse requirements and implementation work before submission

- Reuse accepted aTZ/aQZ geometries and harmonic results when method,
  dispersion, state, atom mapping, stationary-point identity and the target
  component's geometry convention match. L1/L2 labels alone neither enable nor
  prevent reuse.
- Harmonic frequencies do not supply a VPT2 anharmonic correction. Compute only
  missing VPT2 components at the accepted basis-specific minima, then
  extrapolate the anharmonic differences using the declared protocol. Do not
  add two harmonic ZPEs or count rotor contributions twice.
- Current `automatic_composite_spec` accepts a cc-pVTZ L2/VPT2 profile, and
  `anl.anharmonic.extend` builds an unaugmented cc-pVTZ/cc-pVQZ pair. The new
  L1=aTZ, L2=aQZ request requires a general basis-pair/provenance extension and
  tests of adoption from both levels before the R9 input is executable. The
  augmented ANL flag alone does not implement this change.
- B2PLYP electronic energies cannot replace CCSD(T)-F12 components. Reuse native
  coupled-cluster calculations only when their complete method/basis/state/
  geometry contracts match the requested ANL or VRC component.
- The running unaugmented TZ scouting calculations can provide structures,
  connections and starting guesses; they are not numerically identical to
  future augmented-basis results. Preserve prior search outcomes and audit
  their applicability under the changed theory instead of either indiscriminately
  re-exploring them or treating every old rejection as immutable.
- Independent-run import, final reviewed model installation, final selected
  search definitions and complete staged orchestration validation remain
  required before claiming R9 submission readiness.

The older `input.draft.json`, readiness records and scope tables are historical
drafts until these requirements and the new native connection findings are
integrated. No new calculation or MESS run was started to write this note.
