# Acetaldehyde automatic CBH-ANL0-F12 validation

This is one ordinary KinBot input. The user requests `ANL0-F12`; KinBot
generates and balances the CBH network, chooses the highest structurally valid
CBH rung, applies size based electronic structure fallbacks, runs the accepted
tight and barrierless channels, and hands the resulting thermochemistry and
ROTD_py flux to MESS.

Reaction discovery uses the ordinary KinBot defaults: reaction search is on,
all supported unimolecular families are considered, and no family or atom
motif is selected by the user. Optimization, IRC, connectivity, and barrier
checks determine which discovered channels enter the network. The automatic
VRC branch selects accepted homolytic channels from that completed network.
The default `pes: 0` keeps this as a single-starting-well validation; recursive
searches from newly formed wells are intentionally outside this small test.

The input contains no executable path, module name, FairChem cache path, CBH
rung, CBH reference list, or per-program resource table. KinBot checks explicit
configuration, `PATH`, conventional user and system installation prefixes, and
both loaded and available environment modules. It records the resolved paths
and any required module loads in `site_setup.sh`; the user does not need to
load program modules before starting the driver.

After the accepted reaction network is frozen, VRC-TST and CBH/ANL execute as
independent branches. With this eight-slot input, KinBot reserves four slots
for VRC scanning, correction potentials, and ROTD_py, and four for reference
L2 calculations and composite graphs. Reference drivers use a shared
submission lock and count only jobs inside their reserved branch. MESSPF and
MESS are an explicit join: neither starts until both branches have completed
and their hash-bound artifacts are available.

This validation selects `day-long-cpu`, which gives the B2PLYP L2 saddle and
frequency calculations their intended 24-hour execution window. Bundled Slurm
templates, including ROTD_py sample jobs and VRC correction dispatchers,
discover and request the selected partition's maximum wall time. A timeout is
recorded by the normal retry and recipe fallback machinery; work is not
silently routed to a different partition. Existing custom scheduler templates
retain their own directives.

If an optional ANL component exhausts its execution policy, KinBot preserves
the completed reference and corrections, omits only whole additive families,
and records the reduced expression explicitly. It tries that common reduced
recipe across the current CBH reaction before lowering the electronic tier or
the CBH rung. Restarting the same input reuses completed dispatcher records;
it does not resubmit successful calculations.
