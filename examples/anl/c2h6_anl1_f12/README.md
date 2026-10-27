# Ethane / methyl recombination

`c2h6.json` stages an ethane PES exploration including C–C homolysis to
CH3 + CH3. It uses ordinary reaction discovery and CBH selection defaults;
there is no molecule-specific reaction list, geometry, orbital index or root
assignment in the input. This example is separate from the sulfur campaigns.
KinBot requires a barrier threshold; this demo sets the L1 threshold to
150 kcal/mol.

The requested calculation uses augmented ANL1-F12 with the existing
B3LYP/B2PLYP search/refinement preset. Its VRC correction reconstructs the
singlet potential from the RHF-UCCSD(T)-F12b triplet reference at CBS(T,Q) plus
the CASPT2-F12 singlet–triplet gap at CBS(T,Q). These are the
`avtz-f12` / `avqz-f12` basis pair, with the existing 3.7 extrapolation power.
The VRC accuracy preset supplies CASPT2-F12/AVTZ-F12 relaxed scans and
CASPT2-F12/AVDZ-F12 sampling. The explicit gap list is needed because the
current preset defaults to CBS(D,T), and lists are replaced as a whole.

The ten-job live pool covers the PES workers and downstream calculations
together. Ordinary QC, composite and VRC tasks request up to eight processors;
ROTD_py samples request eight. Whole-node exclusivity is already the default.
No node names or executable locations are embedded. This input leaves
pressure-dependent MESS kinetics disabled: bath-gas collision parameters have
not been specified for this demo.

## Native orbital verification is still required

**This is a staged input, not a native-validated VRC example.** The existing
`radical_pair` option requests CAS(2,2) and checks the two-doublet spin contract.
It does **not** implement the requested automatic orbital-probe, verification
and promotion workflow. In particular, it must not be treated as proof that
frontier orbitals retain the intended dissociation character.

The required `orbital_review/surfaces.json` is deliberately absent. The VRC
driver requires that reviewed file before submitting its fragment or scan work;
it cannot silently fall back to the unverified shortcut. Do not remove that
requirement to launch production. No artificial `reviewed: true` model or
dummy checkpoint is supplied. PES discovery and ANL setup are independently
configured, but this example is not ready for a complete VRC run yet.

For CH3 + CH3 the candidate subspace is the two fragment radical orbitals
(the C–C bonding/antibonding pair on association). The eventual native review
must establish their identities and continuity, the matched singlet/triplet
asymptote and CC reference, and radial, angular, return-path and relevant
expanded-space comparisons under this T/Q protocol. Root counts and weights
belong to that evidence. Other channels discovered by the default PES search
need their own applicability checks and models; this request does not approve
every two-doublet fragmentation.

The general automatic CAS(2,2) verification feature remains separate unfinished
work. Once it is implemented and validated, it can replace this file
prerequisite without adding manual orbital recipes to the molecular input.

To check input parsing without starting any calculations, run from the
repository root:

```bash
.venv/bin/python -B - <<'PY'
from kinbot.parameters import Parameters
from kinbot.vrc_correction import spin_splitting_settings
p = Parameters('examples/anl/c2h6_anl1_f12/c2h6.json').par
s = spin_splitting_settings(p)
print('Composite:', p['composite_method'], 'augmented:', p['anl_f12_augmented'])
print('CC reference:', s['reference'])
print('MR gaps:', s['gap_levels'])
print('Shared job cap:', p['queue_live_policy']['max_jobs'])
print('VRC prerequisite:', p['vrc_tst_surface_model_file'])
PY
```
