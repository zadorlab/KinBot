# Explicit reference validation

This directory is outside the installed `kinbot` package. Reference energies,
geometries, rate fits, paper labels and example-specific reproduction tasks are
opt-in validation data. Runtime discovery, selection, recovery and acceptance
never load these catalogues.

From a checkout, compare any PES case with a schema-1 reference network:

```bash
python -B -m validation.pes_compare /absolute/path/to/case \
  --reference /absolute/path/to/reference.json
```

Use `--evidence` for readable comparison plus native geometries, frequencies and
energies, or `--json` for the complete snapshot, including native records and unresolved
identities. The command prints to stdout, opens ASE databases with SQLite
`mode=ro`, disables ASE lock files, and never submits, cancels or alters a run.
Run with `-B` to prevent Python bytecode writes as well.

`data/d4cp04873c_fig6.json` transcribes the supplied Figure 6, with paper labels
and energies kept solely here. Its connectivity comparison distinguishes
complexes from separated products, and reports pending/rejected candidates
separately from accepted discoveries and outer-PES selection. A matching pair
of endpoints is not proof of identical saddle geometry, spin or conformer.

The generic installed ANL comparator now requires an explicit reference:

```bash
python -m kinbot.anl.literature --reference validation/data/anl_2017.json \
  compare-run BENCHMARK /absolute/path/to/composite
python -m validation.rate_compare /absolute/path/to/mess-run \
  --reference validation/data/methyl_recombination.json \
  --reactant REACTANT_LABEL --product PRODUCT_LABEL
```

`anl_source_reproduction.source_reference_spec` is an opt-in historical UHF
reproduction helper. The chosen catalogue task must explicitly declare UHF;
ordinary production ANL retains its RHF/ROHF reference conventions.
