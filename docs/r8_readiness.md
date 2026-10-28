# R8: inexpensive L1/L2 PES scouting

R8 now uses the existing KinBot PES driver, B3LYP-D3(BJ)/cc-pVDZ L1,
B2PLYP-D3(BJ)/cc-pVTZ L2, and automatic MESS input assembly with default PST.
No ANL, Molpro, VRC scans, reviewed CAS models, or rotdPy execution is requested.

Inputs and staging/launch commands are in
[the scouting handoff](../examples/anl/ch3s_anl1_f12/scouting_r8/README.md).
MESS execution and branching analysis follow exploration. The rate analysis
will guide manual channel selection for R9; R8 does not promise that its
full discovered PES is already kinetically trimmed.

The previous production work ledger is preserved as
[R9 readiness](r9_readiness.md), with unresolved items still unresolved.
The original ANL example inputs and orbital evidence remain for that later
production selection. R7 stays at its existing path, untouched.
