"""Declarative ANL expressions over independently validated components.

CBS, core-valence, and scalar-relativistic components are named inputs here.
Every requirement is assembled only from hash-verified native results or an
explicit exact-zero/state datum.  The recipe layer contains no program path,
scheduler, or molecule-specific value.
"""

from __future__ import annotations

from kinbot.anl.extrapolation import ORIGINAL_ANL_POWER
from kinbot.anl.model import ComponentRequirement, CompositeRecipe, ExpressionTerm


# KinBot's production profiles use the fully augmented Dunning bases on every
# atom, including hydrogen.  Keep the executable Molpro spellings beside the
# declarative recipe so generation and assembly cannot drift apart.
ANL0_NATIVE_REFERENCE_BASES = (
    'aug-cc-pVQZ',
    'aug-cc-pV5Z',
)

ANL1_NATIVE_REFERENCE_BASES = (
    'aug-cc-pV5Z',
    'aug-cc-pV6Z',
)

ANL1_F12_NATIVE_REFERENCE_BASES = (
    'cc-pVQZ-F12',
    'cc-pV5Z-F12',
)


# These are independent additive corrections to the mandatory reference CBS
# energy.  Automatic production may omit a whole family after a real task
# failure, but it must never keep only one side of a difference.  The strict
# published recipes remain the default.
ANL_REDUCIBLE_FAMILIES = (
    'higher_order',
    'dboc',
    'scalar_relativistic',
    'core_valence',
)


def _required(key, quantity, method, basis, role='l3',
              backends=('composite',), **settings):
    return ComponentRequirement(key, quantity, method, basis, role,
                                tuple(backends), settings)


def recipe(name: str, *, vpt2_method: str | None = None,
           vpt2_cbs: bool = False, multiplicity: int = 1,
           atomic: bool = False,
           rank_exact_higher_order: bool = False,
           harmonic_only: bool = False,
           omitted_families: tuple[str, ...] = (),
           fallback_geometry_method: str | None = None,
           fallback_harmonic_method: str | None = None,
           fallback_basis: str = 'cc-pVTZ') -> CompositeRecipe:
    """Return the original equation or an explicitly labeled VPT2 variant.

    The published ANL1 anharmonic correction uses B3LYP/cc-pVTZ. The selected
    higher tier uses B2PLYP-D3(BJ); optional TZ/QZ VPT2 extrapolation is also
    a profiled variant. Both choices receive distinct recipe labels.
    """
    if name not in ('ANL0', 'ANL0-F12', 'ANL1', 'ANL1-F12'):
        raise ValueError(f'Unsupported ANL recipe {name!r}.')
    if vpt2_method is None:
        vpt2_method = 'B3LYP'
    if vpt2_method not in ('B3LYP', 'B2PLYP-D3BJ'):
        raise ValueError(f'Unsupported VPT2 surface {vpt2_method!r}.')
    if not isinstance(vpt2_cbs, bool):
        raise ValueError('vpt2_cbs must be boolean.')
    if (isinstance(multiplicity, bool) or not isinstance(multiplicity, int)
            or multiplicity < 1):
        raise ValueError('multiplicity must be a positive integer.')
    state_reference = 'RHF' if multiplicity == 1 else 'ROHF'
    if not isinstance(atomic, bool):
        raise ValueError('atomic must be boolean.')
    if not isinstance(rank_exact_higher_order, bool):
        raise ValueError('rank_exact_higher_order must be boolean.')
    if not isinstance(harmonic_only, bool):
        raise ValueError('harmonic_only must be boolean.')
    if (not isinstance(omitted_families, tuple)
            or len(omitted_families) != len(set(omitted_families))
            or any(item not in ANL_REDUCIBLE_FAMILIES
                   for item in omitted_families)):
        raise ValueError('Invalid omitted ANL component families.')
    geometry_fallback = fallback_geometry_method is not None
    harmonic_fallback = fallback_harmonic_method is not None
    if ((geometry_fallback and not fallback_geometry_method)
            or (harmonic_fallback and not fallback_harmonic_method)
            or ((geometry_fallback or harmonic_fallback)
                and not fallback_basis)):
        raise ValueError('Invalid L2 geometry or harmonic fallback identity.')
    variants = []
    if name in ('ANL0-F12', 'ANL1-F12'):
        # The 2017 article names F12b but does not specify SCALE_TRIP.
        # Our explicit SCALE_TRIP=1 must not be labeled an exact reproduction.
        variants.append('scaled-triples')
    if name == 'ANL0-F12':
        variants.append('dual-CBS-high-reference')
    if name == 'ANL1-F12':
        variants.append('Ram-2025-high-reference')
    if vpt2_method == 'B2PLYP-D3BJ' and not harmonic_only:
        variants.append('B2PLYP-D3BJ')
    if vpt2_cbs and not harmonic_only:
        variants.append('VPT2-CBS')
    if harmonic_only:
        variants.append('harmonic-ZPE')
    if geometry_fallback:
        variants.append(f'{fallback_geometry_method}/{fallback_basis}-geometry')
    if harmonic_fallback:
        variants.append(
            f'{fallback_harmonic_method}/{fallback_basis}-ZPE')
    if omitted_families:
        ordered_omissions = tuple(
            item for item in ANL_REDUCIBLE_FAMILIES
            if item in omitted_families)
        variants.append('reduced[without=' + '+'.join(ordered_omissions) + ']')
    label = (f'profiled:{name}:' + ':'.join(variants)) if variants else name

    if name == 'ANL0-F12':
        # User-selected profiled extension: use the same dual-CBS reference
        # construction as ANL1-F12, but retain the ANL0 geometry, harmonic,
        # and higher-order components.
        reference = _required(
                              'reference_cbs', 'electronic',
                              'mean[CCSD(T)-F12b,CCSD(T)]',
                              ('mean[CBS(QZ-F12,5Z-F12),'
                               'CBS(aug-5Z,aug-6Z)]'),
                              scale_trip=1,
                              correlation='unrestricted',
                              reference=state_reference,
                              reference_convention='RHF/ROHF-UCC',
                              f12_native_basis_pair=
                                  ANL1_F12_NATIVE_REFERENCE_BASES,
                              f12_upper_cardinal=5,
                              f12_extrapolation_power=ORIGINAL_ANL_POWER,
                              conventional_native_basis_pair=
                                  ANL1_NATIVE_REFERENCE_BASES,
                              conventional_upper_cardinal=6,
                              conventional_extrapolation_power=
                                  ORIGINAL_ANL_POWER,
                              component_weights=(0.5, 0.5),
                              geometry_method=(fallback_geometry_method
                                               or 'CCSD(T)'),
                              geometry_basis=(fallback_basis if geometry_fallback
                                              else 'cc-pVTZ'))
    elif name == 'ANL1-F12':
        # This is a KinBot extension, not a fourth named method in the 2017
        # paper.  Its reference follows the highest-level robust practice in
        # Ram et al.: average the QZ-F12/5Z-F12 and conventional 5Z/6Z
        # CBS estimates.  KinBot deliberately retains the ANL inverse-l^3.7
        # total-energy extrapolation for both estimates.
        reference = _required(
                              'reference_cbs', 'electronic',
                              'mean[CCSD(T)-F12b,CCSD(T)]',
                              ('mean[CBS(QZ-F12,5Z-F12),'
                               'CBS(aug-5Z,aug-6Z)]'),
                              scale_trip=1,
                              correlation='unrestricted',
                              reference=state_reference,
                              reference_convention='RHF/ROHF-UCC',
                              f12_native_basis_pair=
                                  ANL1_F12_NATIVE_REFERENCE_BASES,
                              f12_upper_cardinal=5,
                              f12_extrapolation_power=ORIGINAL_ANL_POWER,
                              conventional_native_basis_pair=
                                  ANL1_NATIVE_REFERENCE_BASES,
                              conventional_upper_cardinal=6,
                              conventional_extrapolation_power=
                                  ORIGINAL_ANL_POWER,
                              component_weights=(0.5, 0.5),
                              geometry_method=(fallback_geometry_method
                                               or 'CCSD(T)'),
                              geometry_basis=(fallback_basis if geometry_fallback
                                              else 'cc-pVQZ'))
    elif name == 'ANL1':
        reference = _required('reference_cbs', 'electronic', 'CCSD(T)',
                              'CBS(aug-5Z,aug-6Z)//cc-pVQZ',
                              correlation='unrestricted',
                              reference=state_reference,
                              reference_convention='RHF/ROHF-UCC',
                              extrapolation_power=ORIGINAL_ANL_POWER,
                              upper_cardinal=6,
                              native_basis_pair=ANL1_NATIVE_REFERENCE_BASES,
                              geometry_method=(fallback_geometry_method
                                               or 'CCSD(T)'),
                              geometry_basis=(fallback_basis if geometry_fallback
                                              else 'cc-pVQZ'))
    else:
        reference = _required('reference_cbs', 'electronic', 'CCSD(T)',
                              'CBS(aug-QZ,aug-5Z)//cc-pVTZ',
                              correlation='unrestricted',
                              reference=state_reference,
                              reference_convention='RHF/ROHF-UCC',
                              extrapolation_power=ORIGINAL_ANL_POWER,
                              upper_cardinal=5,
                              native_basis_pair=ANL0_NATIVE_REFERENCE_BASES,
                              geometry_method=(fallback_geometry_method
                                               or 'CCSD(T)'),
                              geometry_basis=(fallback_basis if geometry_fallback
                                              else 'cc-pVTZ'))

    harmonic_basis = (fallback_basis if harmonic_fallback else
                      ('CBS(cc-pVTZ,cc-pVQZ)' if name in ('ANL1', 'ANL1-F12')
                       else 'cc-pVTZ'))
    harmonic_method = fallback_harmonic_method or 'CCSD(T)'
    harmonic_backends = (('kinbot',) if harmonic_fallback else
                         (('known_zero',) if atomic else
                          ('composite',)
                          if name in ('ANL1', 'ANL1-F12') else ('molpro',)))
    harmonic_settings = ({} if harmonic_fallback else
                         ({'extrapolation_power': ORIGINAL_ANL_POWER,
                           'upper_cardinal': 4,
                           'geometry_mode': 'basis_optimized',
                           'correlation': 'unrestricted',
                           'reference': state_reference,
                           'reference_convention': 'RHF/ROHF-UCC'}
                          if name in ('ANL1', 'ANL1-F12')
                          else {'correlation': 'unrestricted',
                                'reference': state_reference,
                                'reference_convention': 'RHF/ROHF-UCC'}))
    vpt2_settings = ({'dispersion': 'GD3BJ'} if vpt2_method == 'B2PLYP-D3BJ'
                     else {'dispersion': ''})
    if vpt2_cbs:
        vpt2_settings.update(extrapolation_power=ORIGINAL_ANL_POWER,
                             upper_cardinal=4,
                             geometry_mode='basis_optimized')
    vpt2_basis = ('CBS(cc-pVTZ,cc-pVQZ)' if vpt2_cbs else 'cc-pVTZ')
    vpt2_backends = ('composite',) if vpt2_cbs else ('gaussian',)
    common = [
        reference,
        _required('harmonic_zpe', 'zpe', harmonic_method, harmonic_basis,
                  role=('l2' if harmonic_fallback else 'l3'),
                  backends=harmonic_backends,
                  **harmonic_settings),
        _required('core_valence_cbs', 'correction',
                  'CCSD(T,full)-CCSD(T,frozen-core)',
                  'CBS(cc-pCVTZ,cc-pCVQZ)',
                  backends=('composite', 'known_zero'),
                  correlation='unrestricted',
                  reference=state_reference,
                  reference_convention='RHF/ROHF-UCC',
                  extrapolation_power=ORIGINAL_ANL_POWER,
                  upper_cardinal=4),
        _required('scalar_relativistic', 'correction', 'CCSD(T)-DKH delta',
                  'aug-cc-pCVTZ-DK', correlation='unrestricted',
                  reference=state_reference,
                  reference_convention='RHF/ROHF-UCC'),
        _required('dboc', 'correction', 'HF', 'cc-pVTZ',
                  backends=('cfour',)),
        _required('spin_orbit', 'correction', 'SO', 'state-specific',
                  role='state', backends=('known_zero', 'table',
                                          'calculated', 'manual')),
    ]
    if not harmonic_only:
        common.insert(2, _required(
            'vpt2_correction', 'correction', vpt2_method, vpt2_basis,
            role='l2', backends=(('known_zero',) if atomic else
                                  vpt2_backends), **vpt2_settings))
    high_backend = 'mrcc'
    high_settings = {
        'correlation': 'unrestricted', 'core': 'frozen',
        'reference': state_reference,
        'reference_convention': 'RHF/ROHF-UCC',
    }
    qp_settings = {
        'correlation': 'unrestricted', 'core': 'frozen',
        'reference': state_reference,
        'reference_convention': 'RHF/ROHF-UCC',
    }
    if rank_exact_higher_order:
        method = (('CCSDT(Q)-CCSD(T)@cc-pVTZ + '
                   'CCSDTQ(P)-CCSDT(Q)@cc-pVDZ')
                  if name in ('ANL1', 'ANL1-F12') else
                  'CCSDT(Q)-CCSD(T)@cc-pVDZ')
        basis = ('cc-pVTZ/cc-pVDZ'
                 if name in ('ANL1', 'ANL1-F12') else 'cc-pVDZ')
        higher = [
            _required(
                'higher_order_correction', 'correction', method, basis,
                backends=('known_zero',),
                rank_exact_electron_maximum=2),
        ]
        hoe_terms = (ExpressionTerm('higher_order_correction'),)
    elif name in ('ANL1', 'ANL1-F12'):
        higher = [
            _required('hoe_tz_high', 'electronic', 'CCSDT(Q)', 'cc-pVTZ',
                      backends=(high_backend, 'known_zero'), **high_settings),
            _required('hoe_tz_low', 'electronic', 'CCSD(T)', 'cc-pVTZ',
                      backends=('molpro',), correlation='unrestricted',
                      core='frozen', reference=state_reference,
                      reference_convention='RHF/ROHF-UCC'),
            _required('hoe_dz_high', 'electronic', 'CCSDTQ(P)', 'cc-pVDZ',
                      backends=('mrcc', 'known_zero'), **qp_settings),
            _required('hoe_dz_low', 'electronic', 'CCSDT(Q)', 'cc-pVDZ',
                      backends=(high_backend, 'known_zero'), **high_settings),
        ]
        hoe_terms = (ExpressionTerm('hoe_tz_high'),
                     ExpressionTerm('hoe_tz_low', -1),
                     ExpressionTerm('hoe_dz_high'),
                     ExpressionTerm('hoe_dz_low', -1))
    else:
        higher = [
            _required('hoe_high', 'electronic', 'CCSDT(Q)', 'cc-pVDZ',
                      backends=(high_backend, 'known_zero'), **high_settings),
            _required('hoe_low', 'electronic', 'CCSD(T)', 'cc-pVDZ',
                      backends=('molpro',), correlation='unrestricted',
                      core='frozen', reference=state_reference,
                      reference_convention='RHF/ROHF-UCC'),
        ]
        hoe_terms = (ExpressionTerm('hoe_high'), ExpressionTerm('hoe_low', -1))
    electronic = (ExpressionTerm('reference_cbs'), *hoe_terms,
                  ExpressionTerm('core_valence_cbs'),
                  ExpressionTerm('scalar_relativistic'),
                  ExpressionTerm('dboc'), ExpressionTerm('spin_orbit'))
    zero_point = ((ExpressionTerm('harmonic_zpe'),)
                  if harmonic_only else
                  (ExpressionTerm('harmonic_zpe'),
                   ExpressionTerm('vpt2_correction')))
    family_keys = {
        'higher_order': {item.key for item in higher},
        'core_valence': {'core_valence_cbs'},
        'scalar_relativistic': {'scalar_relativistic'},
        'dboc': {'dboc'},
    }
    omitted_keys = set().union(*(
        family_keys[family] for family in omitted_families)) \
        if omitted_families else set()
    requirements = tuple(
        item for item in common + higher if item.key not in omitted_keys)
    electronic = tuple(
        item for item in electronic if item.component not in omitted_keys)
    return CompositeRecipe(label, 1, requirements, electronic, zero_point)
