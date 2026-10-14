# CH3S ANL1-F12 and VRC-TST example

This ordinary KinBot input starts from neutral doublet methylthiyl. Reaction
families, accepted products, the CBH rung, and feasible ANL tiers are selected
by KinBot. The large barrier cutoff prevents energy pruning during this
research validation.

The VRC settings are intentionally expensive. The relaxed path and frozen
sampling surface use a 13-electron, 8-orbital CASPT2-F12b model. Sampling is
extrapolated from aug-cc-pVTZ-F12/aug-cc-pVQZ-F12. The trusted doublet curve is
constructed from a quartet RHF/ROHF-UCCSD(T)-F12b reference plus doublet minus
quartet CASPT2-F12b and Davidson-corrected MRCI-F12b+Q gaps. Each branch uses
the 3.7-power QZ/5Z extrapolation and the distance-dependent CASPT2/MRCI
weighting stored in the final correction record.

The ANL1-F12 electronic reference also opts into fully augmented F12 bases:
aug-cc-pVQZ-F12/aug-cc-pV5Z-F12 is averaged with the conventional fully
augmented aug-cc-pV5Z/aug-cc-pV6Z limit. Both branches use the 3.7 power.

The 13e,8o space is a deliberately broad high-accuracy CH3S space. Before a
production calculation, inspect the natural occupations over the full scan;
an active space is a scientific input and cannot be inferred safely from atom
counts alone.
