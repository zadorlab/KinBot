# Reverse-search complex optimization recovery

A reverse bimolecular search can retain a tight transition state while the
associated weak entrance or product complex fails to optimize. Such an optimizer
failure is distinct from a proposed saddle search that finds no reaction.

Once the reverse TS has passed L2 validation, its owning driver enables one
additional recovery for an associated complex: an exhausted, unconstrained
internal-coordinate Sella minimum can restart in Cartesian coordinates. Both
entrance and product complexes use the same recovery procedure. If the product
complex finishes unsuccessfully before its TS, it waits for the TS verdict
before deciding whether this recovery applies.

The recovery requires a matching generated minimum worker, finite trajectory
geometry/energy/forces, agreement with the last evaluated optimizer log point,
and retained allocation ownership. It starts at that verified evaluated point.
The new attempt has the usual bounded optimizer budget (by default three chunks
of 250 steps). It does not extend or repeat an already exhausted Cartesian
attempt. Existing SCF and ODE recovery procedures remain available for their
respective failure classes.

The original attempt is archived by the ordinary QC retry transaction. Method,
basis, charge, multiplicity, force threshold, frequency validation, and final
geometry/connectivity checks are unchanged. The accepted TS, fragments, sibling
calculations and live allocations are preserved. The owning optimizer adopts
the replacement result; no synthetic successful database row is installed.

Missing saddles, Hessian rejection, unverified trajectories and arbitrary native
errors do not qualify for this optimizer fallback. Exhausting the fallback is
still an unsuccessful optimization, not evidence that the proposed complex or
reaction is physically absent.

An explicit QC retry can also use the existing live control with
`worker_parameters.sella_kwargs.internal=false`. After Slurm has purged a
finished attempt, this requires a successful scheduler query confirming absence,
a generation-tagged nonzero exit marker, the matching failed execution guard,
and the retained allocation identity. An unavailable scheduler does not satisfy
these checks.

Validation covers failed/unfinished TSs, unrelated failure classes, repeat
invocation, delayed admission, retained evidence, and result adoption. Local
replay of an exported reverse worker reached only its eligible complex retry
through the native controller, with submission intercepted. This verifies the
orchestration, not the convergence of a calculation on the remote cluster.
