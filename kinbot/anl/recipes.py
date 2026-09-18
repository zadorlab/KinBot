"""Declarative ANL expressions over independently validated components.

CBS, core-valence, and scalar-relativistic components are named inputs here.
Verified Molpro task pairs can provide CBS references; the remaining native
QC providers still require implementation.
"""

from __future__ import annotations

from kinbot.anl.extrapolation import ORIGINAL_ANL_POWER
from kinbot.anl.model import ComponentRequirement, CompositeRecipe, ExpressionTerm


def _required(key, quantity, method, basis, role='l3',
              backends=('composite',), **settings):
    return ComponentRequirement(key, quantity, method, basis, role,
                                tuple(backends), settings)


def recipe(name: str, *, vpt2_method: str | None = None) -> CompositeRecipe:
    """Return the original equation or an explicitly labeled VPT2 variant.

    The published ANL1 correction uses B3LYP. The selected higher tier uses
    B2PLYP-D3(BJ), so that choice gets a distinct profiled recipe label.
    """
    if name not in ('ANL0', 'ANL0-F12', 'ANL1'):
        raise ValueError(f'Unsupported ANL recipe {name!r}.')
    if vpt2_method is None:
        vpt2_method = 'B3LYP'
    if vpt2_method not in ('B3LYP', 'B2PLYP-D3BJ'):
        raise ValueError(f'Unsupported VPT2 surface {vpt2_method!r}.')
    variants = []
    if name == 'ANL0-F12':
        # The 2017 article names F12b but does not specify SCALE_TRIP.
        # Our explicit SCALE_TRIP=1 must not be labeled an exact reproduction.
        variants.append('scaled-triples')
    if vpt2_method == 'B2PLYP-D3BJ':
        variants.append('B2PLYP-D3BJ')
    label = (f'profiled:{name}:' + ':'.join(variants)) if variants else name

    if name == 'ANL0-F12':
        reference = _required('reference_cbs', 'electronic', 'CCSD(T)-F12b',
                              'CBS(cc-pVTZ-F12,cc-pVQZ-F12)', scale_trip=1,
                              extrapolation_power=ORIGINAL_ANL_POWER,
                              upper_cardinal=4)
    elif name == 'ANL1':
        reference = _required('reference_cbs', 'electronic', 'CCSD(T)',
                              "CBS(a'5Z,a'6Z)//cc-pVQZ",
                              extrapolation_power=ORIGINAL_ANL_POWER,
                              upper_cardinal=6)
    else:
        reference = _required('reference_cbs', 'electronic', 'CCSD(T)',
                              "CBS(a'QZ,a'5Z)//cc-pVTZ",
                              extrapolation_power=ORIGINAL_ANL_POWER,
                              upper_cardinal=5)

    harmonic_basis = ('CBS(cc-pVTZ,cc-pVQZ)' if name == 'ANL1'
                      else 'cc-pVTZ')
    harmonic_settings = ({'extrapolation_power': ORIGINAL_ANL_POWER,
                          'upper_cardinal': 4}
                         if name == 'ANL1' else {})
    vpt2_settings = ({'dispersion': 'GD3BJ'} if vpt2_method == 'B2PLYP-D3BJ'
                     else {'dispersion': ''})
    common = [
        reference,
        _required('harmonic_zpe', 'zpe', 'CCSD(T)', harmonic_basis,
                  backends=(('composite',) if name == 'ANL1' else ('molpro',)),
                  **harmonic_settings),
        _required('vpt2_correction', 'correction', vpt2_method, 'cc-pVTZ',
                  role='l2', backends=('gaussian',), **vpt2_settings),
        _required('core_valence_cbs', 'correction',
                  'CCSD(T,full)-CCSD(T,frozen-core)',
                  'CBS(cc-pcVTZ,cc-pcVQZ)',
                  extrapolation_power=ORIGINAL_ANL_POWER,
                  upper_cardinal=4),
        _required('scalar_relativistic', 'correction', 'CCSD(T)-DKH delta',
                  'aug-cc-pcVTZ-DK'),
        _required('dboc', 'correction', 'HF', 'cc-pVTZ',
                  backends=('cfour',)),
        _required('spin_orbit', 'correction', 'SO', 'state-specific',
                  role='state', backends=('known_zero', 'table',
                                          'calculated', 'manual')),
    ]
    if name == 'ANL1':
        higher = [
            _required('hoe_tz_high', 'electronic', 'CCSDT(Q)', 'cc-pVTZ',
                      backends=('cfour', 'mrcc')),
            _required('hoe_tz_low', 'electronic', 'CCSD(T)', 'cc-pVTZ',
                      backends=('molpro',)),
            _required('hoe_dz_high', 'electronic', 'CCSDTQ(P)', 'cc-pVDZ',
                      backends=('mrcc',)),
            _required('hoe_dz_low', 'electronic', 'CCSDT(Q)', 'cc-pVDZ',
                      backends=('cfour', 'mrcc')),
        ]
        hoe_terms = (ExpressionTerm('hoe_tz_high'),
                     ExpressionTerm('hoe_tz_low', -1),
                     ExpressionTerm('hoe_dz_high'),
                     ExpressionTerm('hoe_dz_low', -1))
    else:
        higher = [
            _required('hoe_high', 'electronic', 'CCSDT(Q)', 'cc-pVDZ',
                      backends=('cfour', 'mrcc')),
            _required('hoe_low', 'electronic', 'CCSD(T)', 'cc-pVDZ',
                      backends=('molpro',)),
        ]
        hoe_terms = (ExpressionTerm('hoe_high'), ExpressionTerm('hoe_low', -1))
    electronic = (ExpressionTerm('reference_cbs'), *hoe_terms,
                  ExpressionTerm('core_valence_cbs'),
                  ExpressionTerm('scalar_relativistic'),
                  ExpressionTerm('dboc'), ExpressionTerm('spin_orbit'))
    zero_point = (ExpressionTerm('harmonic_zpe'),
                  ExpressionTerm('vpt2_correction'))
    return CompositeRecipe(label, 1, tuple(common + higher),
                           electronic, zero_point)
