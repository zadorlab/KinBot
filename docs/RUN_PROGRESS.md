# Read-only phase and timing evidence

Run `python -m kinbot.run_progress RUN --output /tmp/fresh-report.json`.
The output file must be new and outside RUN. The compact terminal summary
groups ANL tasks by method, reports the saved VRC workload and current optimizer
progress, and shows whether correction graphs and PES-level ROTD outputs exist.
The JSON retains individual graph dependencies, native timing text, optimizer
tails, the saved CBH network request, and the PES worker summary.

This command uses the standard library and one read-only `squeue` query. It
does not invoke dispatcher reconciliation/preflight, load pickles or QC
databases, hash wavefunctions, launch QC, or alter the running controllers.
Native/log reads are bounded to 256 KiB per selected file; only one current
native tail is read per task. Reported completions are recorded status, not a
new scientific validation. Read errors and scheduler errors are explicit.

The reporter reads the saved input and printed effective scan grid. Missing
scan settings stay unknown; current library defaults must not be substituted
silently for an older run. Scan graph counts are records, not necessarily
unique accepted points: alternate input fingerprints can coexist. Correction
task counts use `current_dispatch.json`, rather than summing obsolete graphs.
The VRC restart journal identifies retired scan graphs. Those remain in the
JSON as comparison evidence but are excluded from current scan counts and
native-log reads. MR fragment-asymptote optimizations and fixed-geometry
orbital tasks are counted separately from finite-distance scan points.
Each configured bond, electronic root and scan basis is printed explicitly.
An optimizer step number is an iteration within one point, not a scan-point
index. Task counts remain observations, not an exact remaining-work total.

## Interpreting timing

Task fractions are not wall-time fractions. Geometry optimization, numerical
gradients, basis sizes and correlation ranks can have very different costs.
File modification time measures output freshness, not computation duration.
Native REAL TIME lines can describe individual modules, not whole allocations;
the report preserves their source text rather than assuming a total runtime.

In the current full-PES workflow, each multireference VRC leaf performs its
relaxed scan sequentially, carrying geometry/orbitals from one point to the
next. It then dispatches the independent correction single points and assembles
that leaf's correction. State leaves/channels are processed sequentially by
the worker. With `vrc_tst_fragment_geometry=asymptotic_scan`, each leaf first
optimizes the separated supermolecule at the configured asymptotic distance
using the MR scan level, then extracts fragment geometries and prepares
fixed-geometry auxiliary orbitals when HOMO pivots are requested. It reuses
that optimized structure for the scan's asymptote. The legacy fragment mode
instead constructs an explicit separated structure at that endpoint.
ROTD dividing-surface distances and relaxed scan distances are separate grids.

CBH/ANL proceeds concurrently with VRC. For `pes` runs, workers prepare ROTD
inputs with sampling deferred. Actual sampling is launched by PES
postprocessing after the workers finish; discovered waiting wells can extend
this gate. Within a PES, scalar ROTD inputs are executed in sequence, with
their sampling allocations using the shared capacity. The separate PES runs
need not finish together. Final MESS rate calculations need the accepted
thermochemistry and required ROTD results.

An ETA needs representative accepted-point durations, optimizer progress,
remaining correlation tiers, and (once sampling starts) measured sample
throughput and statistical convergence. A scenario such as remaining points
times an assumed mean point duration is useful for planning, but is not a
predicted completion date. Future PES discoveries and adaptive Monte Carlo
work are not fully enumerated by today's materialized graphs.
