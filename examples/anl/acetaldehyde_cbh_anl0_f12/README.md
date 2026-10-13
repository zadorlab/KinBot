# Acetaldehyde automatic CBH-ANL0-F12 validation

This is one ordinary KinBot input. The user requests `ANL0-F12`; KinBot
generates and balances the CBH network, chooses the highest structurally valid
CBH rung, applies size based electronic structure fallbacks, runs the accepted
tight and barrierless channels, and hands the resulting thermochemistry and
ROTD_py flux to MESS.

The three atom motifs deliberately limit this validation to C-C homolysis,
the adjacent-H molecular channel to methane plus CO, and H2 elimination to
ketene. They scope the PES search and do not choose a CBH rung.

The input contains no executable path, module name, FairChem cache path, CBH
rung, CBH reference list, or per-program resource table. Site software remains
discoverable through the loaded environment and `site_setup.sh`; one
`queue_job_limit` caps all automatic calculation pools unless an expert
explicitly overrides a subsystem.

If an optional ANL component exhausts its execution policy, KinBot preserves
the completed reference and corrections, omits only whole additive families,
and records the reduced expression explicitly. It tries that common reduced
recipe across the current CBH reaction before lowering the electronic tier or
the CBH rung. Restarting the same input reuses completed dispatcher records;
it does not resubmit successful calculations.
