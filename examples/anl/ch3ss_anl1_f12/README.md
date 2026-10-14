# CH3SS ANL1-F12 and VRC-TST example

This input starts from neutral doublet methyl perthiyl and otherwise uses the
same hands-off automatic network, ANL ladder, ROTD_py, MESSPF, and MESS path as
the CH3S example.

The sulfur multireference model uses the targeted 9-electron, 7-orbital CH3SS
space. Doublet and quartet surfaces are allowed for a doublet plus triplet
fragment pair. Sampling uses CASPT2-F12b with a TZ/QZ CBS extrapolation; the
trusted correction uses QZ/5Z RHF/ROHF-UCCSD(T)-F12b, CASPT2-F12b, and
Davidson-corrected MRCI-F12b+Q branches with the common 3.7-power expression.
The ANL1-F12 reference sets `anl_f12_augmented` so its F12 branch is the
fully augmented QZ-F12/5Z-F12 pair; the conventional 5Z/6Z branch is also
fully augmented, including hydrogen.

Natural occupations and orbital continuity must still be reviewed for a new
dissociation topology. KinBot validates electron parity, spin coupling,
method identity, F12 variant B, Davidson selection, hashes, and complete CBS
pairs; it does not claim that a chemically chosen active space is universal.
