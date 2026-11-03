# Native MESS numerical investigation

These are **auditable experimental patches**, not an approved replacement for
every MESS installation or a claim that every finite-pressure rate is accurate.
KinBot does not install or select them automatically. The corrected chemical
models, all native inputs, unsuccessful controls, and signed outputs remain
under `local_runs/mess_review_20261102`.

The scouting reference uses `WellCutoff=-1` (no explicit truncation). Cutoff 10
is a comparison, not a required repair. A tiny signed coefficient need not
prevent a scouting decision when its effect on the relevant branching is
bounded. Conversely, positive coefficients can still be grossly wrong. Missing
`***` entries are never filled with zeros.

## Source accounting

Original source archive SHA256:
`911085be5322c7de2686ac789d0f3a2516f70b971c19ee88f55d542f6a45c908`.
`experimental_source_manifest.json` records patch hashes, individual source
hashes, executable hashes, and line counts. Counts below are added/deleted
lines against each patch's stated parent, not the size of the original file.

| Patch and parent | Changed files, +added / −deleted | Purpose and disposition |
|---|---|---|
| `compiler_compatibility.patch`, original source | `mpack.hh` +11/−4; `mess.cc` +10/−1 | QD/DD scalar conversions and compiler/library interoperability. Does not increase matrix-assembly precision. |
| `local_accelerate_lp64.patch`, original source | `lapack.hh` +1/−1 | Match the local LP64 library's integer ABI. **Never apply this to an ILP64 build.** |
| `finite_pressure_experimental.patch`, compatible LP64 baseline | `mess.cc` +55/−24; `mess_driver.cc` +24/−0; `cool.cc` +8/−2; `energy_grained_export.hh` +51/−0 | Initial-population normalization, consistent time integration, closed-system guard, hot-energy coverage, diagnostic export and output precision. |
| `direct_hybrid_experimental.patch`, finite-pressure variant | `mess.cc` +31/−5; `mess_driver.cc` +7/−0 | Exercise the commented low-eigenvalue substitution and report residual/orthogonality diagnostics. **Not accepted as a general fix.** |
| `direct_conservative_experimental.patch`, hybrid variant | `mess.cc` +16/−1; `mess.hh` +2/−0; `mess_driver.cc` +5/−0; `direct_conservative_eigen.hh` +72/−0 | Opt-in high-precision assembly that preserves edge balance before direct diagonalization. Native CH3S results improve; broader qualification remains separate. |
| `flux_cutoff_experimental.patch`, conservative variant | `mess.cc` +11/−0; `mess.hh` +1/−0; `mess_driver.cc` +10/−0; `well_cutoff_flux.hh` +49/−0 | Opt-in bound on thermal barrier flux removed by explicit truncation. A diagnostic prototype, unnecessary for the chosen cutoff −1 scouting reference. |
| `flux_labels_diagnostic.patch`, finite-pressure variant | `energy_grained_export.hh` +4/−4 | Add the native barrier identifiers to exported edges and sinks. All exported numerical records were checked against the earlier full-model export and are unchanged. |

The compatible baseline was compared against all 233 matching upstream source
files. Only `lapack.hh`, `mpack.hh`, and `mess.cc` differed. Build recipes and
link commands are saved under `finite_pressure/build_provenance`; local absolute
library paths are evidence of the tested build, not portable installation
instructions. Native `mplapack_dd.cc` is an upstream file, not a new solver
written in this investigation.

## What the changes mean mathematically

**Library ABI.** A Fortran LAPACK `INTEGER` argument must have the width expected
by the linked library. The original local build passed 64-bit storage to an
LP64 interface. The independent small-matrix probe reproduced an invalid
status value, `4294967296`, for a valid positive-definite problem. Matching
32-bit declarations to that library fixes the interface, not the chemistry.
Both Cholesky and pivoted-solve probes pass with the corrected interface.

**Initial normalization.** The native reduced propagator computed `nfac` but
did not apply it to its retained initial state and eliminated-mode escape
contribution. Both must describe the same unit initial population:

```cpp
state /= nfac;
if (bimolecular_size()) escape /= nfac;
```

This is an initial-condition correction. It does not renormalize the final
answer to conceal mass loss.

**Product integration.** For population equation `p'=-Kp` and product flux
`y'=Fp`, an explicit midpoint step must use the same midpoint for both:

```cpp
midpoint = state - (tstep / 2.) * (kin_mat * state);
escape += tstep * (escape_operator * midpoint);
state  -= tstep * (kin_mat * midpoint);
```

The old escape update used the starting population while the well update used
the midpoint. If the operators conserve probability, the consistent updates
cancel their losses/gains algebraically. The empty-escape guard permits a
closed multiwell system, which previously attempted to scale an uninitialized
escape vector. Caching the columns of the linear `escape_flux` operator changes
evaluation cost, not its definition. Native time-step checks and a closed
three-well regression are retained.

**Hot-energy support.** Include each requested initial energy in the grid,
subject to the existing model-energy limit, and reject an out-of-range energy
explicitly. The previous silent omission of a requested energy is not a valid
zero branching fraction. The control retains headroom above that initial
energy. It changes the numerical domain, not the stationary-point energies.

**Direct matrix assembly.** Write each reversible edge as its nonnegative
conductance contribution in the symmetric representation, using high precision
for the paired diagonal and off-diagonal terms. With `b_i` the equilibrium
square-root weight and symmetric off-diagonal magnitude `c`, add

```text
A_ii += c b_j/b_i;  A_jj += c b_i/b_j;
A_ij -= c;         A_ji -= c.
```

Then this edge annihilates `b` exactly in real arithmetic. Add physical product
losses separately. Converting an already rounded double matrix to QD cannot
recover this relation when slow chemical eigenvalues are far below collision
frequencies. The opt-in implementation reconstructs balance from native edge
rates, reports the difference from the original matrix, and neither clips
eigenvalues nor edits output signs. Double-double is not automatically adequate;
CH3S at 300 K required the stronger QD control in this investigation. Returned
eigenvectors and subsequent species reduction still use ordinary doubles.
Radiative models are rejected by this prototype rather than silently omitted.

The scale separation is concrete: the 300 K, 760 torr CH3S control has a
collision frequency about `1.21e10 s^-1`, while its slow interwell rate is
about `4.1e-13 s^-1`. Multiplying that collision frequency by ordinary double
machine epsilon gives about `2.7e-6 s^-1`. This illustrates why rounding at the
fast scale can overwhelm the slow rate; it is not a rigorous eigensolver error
bound. A sign error near that absolute scale can be irrelevant to a picosecond
hot dissociation and still invalidate a claimed thermal lifetime. Convergence
has to be attached to the observable and condition, not to one universal
absolute negative-rate tolerance.

**Explicit cutoff protection.** The optional prototype sums native
`N(E) exp(-E/kT)` on the original energy grid and extends the well's retained
range until the omitted fraction is bounded. Independent sums over exported
reaction and product edges passed for both systems at 300 and 1000 K, for
requested bounds `1e-4` and `1e-6`. This bounds removed thermal flux only; it does
not prove eigenvalue accuracy, Markovian species reduction, or all physical
models. Global cutoff, additional escape flux and well extension are unsupported;
the prototype has not been qualified with `MicroRateMax`. It remains experimental.

## Interpretation and outstanding scope

- The restored hybrid substitution still gives enormous spurious CH3SS entries
  in a tested control. Its source is preserved to document a rejected option;
  it is not silently selected as a cure.
- A separate native alternative-driver calculation gave positive coefficients
  but disagreed substantially with direct and unreduced-generator benchmarks.
  Positivity alone is insufficient to approve it.
- Native high-precision direct CH3S controls at 300/1000 K and 1 atm remove the
  gross rate/loss inconsistency. Remaining tiny negative couplings are retained
  as signed numerical evidence. This is a bounded numerical result, not a
  guarantee of exact chemical rate constants.
- Both QD well-reduction controls completed the requested 10 temperatures
  (300–2000 K) and six pressures (0.01, 0.1, 1, 10, 100, 760 torr). The additional
  CH3S conservative-QD direct control completed all ten temperatures at 760 torr,
  combining the earlier 300/1000 K controls with the other eight temperatures.
  Direct loss rates at 400/500/600/800/1200 K are respectively
  `1.981250959e-7`, `0.001171249295`, `0.4107767022`, `540.439687`, and
  `386056.3649 s^-1`. The well-reduction loss at 400 K is `8.66e-7 s^-1`
  despite an outgoing sum near `1.98e-7 s^-1`: the precision control has not
  repaired that loss summary. Execution, population reduction and significant-rate
  convergence remain separate assessments.
- The CH3SS conservative-DD direct controls with cutoffs −1, 10, 20 and 40 each
  reached the 1800-second test limit. These are unsuccessful bounded controls,
  not evidence of accurate direct rates or a reason to discard channels.
- At 1800/2000 K the CH3S direct control selects a chemical subspace of dimension
  zero; the original-well rate entries are consequently `***`. The CH3SS
  well-reduction control likewise selects dimension zero at 1500–2000 K for the
  tested pressures. These entries must not be populated with invented individual
  well rates. This describes the chosen reduction, not proof that no alternative
  population definition is possible.
- The 193 nm, 300 K, 1 atm initial-population calculation uses the stated helium
  collision model. It is a statistical ground-state-PES scouting calculation,
  not a guarantee of an experimental photochemical quantum yield.

The diagnostic exporter exposes native unreduced transitions before any
eigensolver or species reduction. Independent propagation checks its edge
signs, detailed balance, collision diagonal, probability conservation, energy
grain convergence, and time-step convergence against native results. It does
not replace the scientific model or make a missing phenomenological rate
defined. Existing native initial-projection clipping was not removed by these
patches; the independent unreduced propagation does not clip probabilities.

The ordinary-double direct hot-energy controls at 300 K produced grossly
nonconserving product fractions, despite successful program exits. They are
preserved as unsuccessful controls, not rounded away. The native QD
well-reduction hot-energy results agree with the independently propagated
unreduced generator, for the full models and both pruned candidates. This
supports that particular scouting observable; it does not validate every
thermal rate produced by the same executable.

`audit_export.py` is a standalone read-only check of the experimental export.
It refuses negative/nonfinite rates, invalid equilibrium weights, missing
records and invalid state references; it reports collision-diagonal and
detailed-balance errors without repairing them. Detailed balance is compared
in log space to avoid overflowing `rate * weight**2`. For example:

```bash
python validation/native_mess/audit_export.py /path/to/generator.txt
```

Ten regression cases cover a balanced model with physical product loss,
malformed and signed evidence, conservation violations and very large
equilibrium weights. Passing this file-format/conservation audit alone does
not establish time-step convergence or chemical accuracy.

For the independent propagation, each exported transition `i -> j` enters
`Q[j,i]` and each product escape enters `F[product,i]`. Its diagonal is assembled
as minus the outgoing sum. Thus `1^T Q + 1^T F = 0` by construction in exact
arithmetic; the separately exported native collision diagonal is also checked
against this sum. This is an independent conservative assembly of the native
physical rates, not a claim of bit-identical agreement with the native rounded
matrix and not a repair of a printed phenomenological rate table.

With nonnegative off-diagonals, this finite-state generator preserves positive
populations. The implicit step `p_next = (I-hQ)^(-1) p` is paired with
`y_next = y + h F p_next`, so the population loss and product gain cancel
algebraically. The implementation still checks signs and mass after solving;
floating-point arithmetic is not assumed to preserve this proof exactly.
Time-step and energy-grain controls then test the actual observable. No
comparison can remove errors in the supplied chemical model.

Species reduction is an additional approximation. Projecting fast internal
relaxation onto a small set of well populations need not produce a useful
positive, memoryless species-level generator at every condition. Accordingly,
good unreduced dynamics do not make every reduced coefficient defined, and
the absence of separate slow populations must not be patched by filling
`***` entries. This is distinct from an incorrect numerical-library interface
or loss of floating-point conservation during matrix assembly.

`verify_patches.py` reconstructs all seven recorded source layers in a temporary
directory from the hash-identified upstream archive. It checks patch hashes
and the full source hashes for each build variant, without changing an
installed source tree:

```bash
python validation/native_mess/verify_patches.py /path/to/mess_source.tar.gz --report /tmp/mess-patch-check.json
```

This verifies reproducible source reconstruction, not portable library linking
or numerical correctness. See the [R8 scouting report](../../docs/R8_FINITE_PRESSURE_SCOUTING.md)
for the full-model and combined-pruning results.
