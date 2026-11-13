# Targeted alternative reaction guesses

`run_refresh --reaction-guesses FILE --apply` gives each explicitly selected
settled search one additional attempt. The JSON maps each selected source
name to `{"bond_scale": 1.6, "angle": "parent"}` or a finite list of distinct
guesses in that format.
This option currently applies to three-atom `R_Addition_MultipleBond` motifs.

The constrained search retains the parent motif angle and stretches the
breaking bond toward the selected multiple of its parent length. The final
saddle optimization releases these search constraints. Native frequency,
IRC, L1 screening and retained-channel L2 refinement use the existing checks.
The scale is a starting-guess choice and gives no assurance of a new saddle.

The stopped-controller handoff records the options in the inputs and archives
the former input. The reaction finder gives the alternative an independent
name derived from the source and the options' hash. Existing accepted and
rejected searches, shared product calculations and active allocations retain
their histories. The same completed handoff is idempotent. Its deterministic
alternative names adopt saved terminal outcomes across subsequent restarts.

The example uses scales of 1.3, 1.6 and 1.9 for each of four targets, giving
twelve independent searches. The CH2SH and CH3SS H-loss earlier ordinary
searches converged to a loose fragment saddle and a same-well saddle,
respectively. The alternative changes both the imposed 90-degree approach
and the fixed elemental distance guess. A successful target pathway still
requires a saddle whose IRC connects the intended parent and products.

The CH3SS dissociation companion converged to a separate CH2(SH)S-to-CH2SS+H
pathway. That pathway already exists in the full scout as `rxn_26`, together
with its mirrored counterpart. The targeted retry preserves those routes.

The example also includes HSCH2S to CH2S + SH and HSCH2S to HC(=S)SH + H,
corresponding to the connectivity of paper TS3 and TS5. The latter source
search already supplied a different accepted isomerization. Giving the
alternative its own identity preserves that valid result. A finite unsuccessful
search leaves the requested connection unestablished. It cannot certify a
barrierless pathway or rule out a saddle on a higher-level surface.
