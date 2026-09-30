"""Configured species names layered on top of the existing connectivity ID.

Ordinary names are unchanged. A stereo suffix contains the first 16 hexadecimal
characters of the full identity. Saved references retain the complete identity.
"""
import json
import hashlib
import copy
import logging
from pathlib import Path
import re
from kinbot.stereo_identity import require_supported_identity, UnsupportedStereochemistry


_KEY = re.compile(r'^(\d+)(?:-s([0-9a-f]{16}))?$')
_JOB = re.compile(r'(^|/)(\d+-s[0-9a-f]{16})(?=_|$)')
_MOTIF = re.compile(r'_m\d+(?:-\d+)*(?=_|$)')


def routing_key(species):
    """Return the connectivity key or a configured key; never change ``chemid``."""
    if getattr(species, 'wellorts', 0):
        return species.name
    identity = getattr(species, 'optical_reference', None)
    if identity is None:
        identity = require_supported_identity(species)
    if identity.get('status') != 'assigned':
        raise UnsupportedStereochemistry(
            f'{species.name}: unsupported stereochemical reference: '
            f'{identity.get("reason", "no assigned identity")}')
    if any(
            any(tag in graph for tag in ('@', '/', '\\'))
            for graph in identity.get('canonical_graphs', ())):
        if (identity['id'] != identity['mirror_id']
                and identity['id'][:16] == identity['mirror_id'][:16]):
            raise ValueError(f'{species.name}: the stereoisomer and its mirror have '
                             'different full identities but the same 16-character name suffix.')
        return f"{species.chemid}-s{identity['id'][:16]}"
    return species.chemid


def routing_name(species):
    return str(routing_key(species))


def mess_filename(key, iteration):
    """Keep complete network keys in content without overlong file components."""
    stem = str(key)
    if len(stem.encode()) > 220:
        stem = 'network-' + hashlib.sha256(stem.encode()).hexdigest()
    return f'{stem}_{int(iteration):04d}.mess'


def connectivity_name(name):
    match = _KEY.fullmatch(str(name))
    return match.group(1) if match else str(name)


def is_species_name(name):
    return _KEY.fullmatch(str(name)) is not None


def matches_name(name, choices):
    """Honor connectivity/configured selectors and unsplit reaction names.

    Removing a motif suffix is for selection only, never for QC result reuse.
    """
    connectivity = _JOB.sub(lambda match: match.group(1) + match.group(2).split('-s')[0], str(name))
    return any(candidate in choices for value in (str(name), connectivity)
               for candidate in (value, _MOTIF.sub('', value)))


def configured_selection(mapping, species, default=None):
    """Prefer the configured input selector, retaining connectivity defaults."""
    key = routing_name(species)
    return mapping[key] if key in mapping else mapping.get(str(species.chemid), default)


def same_species(first, second, context='species reuse'):
    """Distinct supported configurations are independent objects, not errors."""
    if first.chemid != second.chemid:
        return False
    left, right = require_supported_identity(first), require_supported_identity(second)
    if routing_key(first) != routing_key(second):
        return False
    return (left['id'] == right['id'] or
            getattr(first, 'optical_population', 'specified') ==
            getattr(second, 'optical_population', 'specified') == 'racemic'
            and left['mirror_id'] == right['id'])


def reusable_cached_product(qc, species):
    """Adopt a verified product's saved atom ordering before optimization.

    A product reached by another H-transfer path may have the same configured
    identity but different H indices. Rebuild an unoptimized product using the
    trusted input order; subsequent normal reads then load coordinates and all
    indexed tensors in that order. Discovery endpoint snapshots stay untouched.
    """
    if not hasattr(qc, 'db') or getattr(species, 'wellorts', 0):
        return species
    from kinbot.stereo_routing import _row_species
    job = routing_name(species) + '_well'
    references = list(qc.db.select(name='stereochemistry/' + job))
    rows = list(qc.db.select(name=job))
    reference = references[-1] if references else None
    row = reference or (rows[-1] if rows and rows[-1].data.get('status') == 'normal' else None)
    if row is None:
        return species
    saved = _row_species(species, row, job, reference)
    if not same_species(species, saved, 'cached product'):
        return species
    import numpy as np
    if (list(species.atom) == list(saved.atom)
            and np.array_equal(species.bond01, saved.bond01)):
        return species
    # Rebuild rather than reorder an existing Hessian, rotor, or conformer
    # object. This helper is used only for initial isolated product objects.
    saved.characterize(bond_mx=saved.bond)
    saved.name = species.name
    for attribute in ('optical_reference', 'optical_population', 'optical_counting_scope'):
        if hasattr(species, attribute):
            setattr(saved, attribute, copy.deepcopy(getattr(species, attribute)))
    logging.getLogger('KinBot').warning(
        '%s: reusing the verified product in its saved atom order; the original '
        'IRC endpoint is retained separately.', job)
    return saved


def require_cache_atom_order(requested, observed, context):
    """Identity is permutation invariant; indexed cache tensors are not remapped."""
    if routing_name(requested) == str(requested.chemid):
        return
    import numpy as np
    elements_match = list(requested.atom) == list(observed.atom)
    graph_match = (requested.chemid != observed.chemid or
                   np.array_equal(requested.bond01, observed.bond01))
    if not elements_match or not graph_match:
        from kinbot.stereo_routing import refuse_routing
        refuse_routing(f'{context}: cached atom indexing differs; restart with the '
                       'original atom ordering or use a fresh calculation directory',
                       [requested, observed])


def input_species(filename):
    from kinbot.stationary_pt import StationaryPoint
    data = json.loads(Path(filename).read_text())
    species = StationaryPoint('saved input', data.get('charge', 0), data.get('mult', 1),
        structure=data.get('structure'), smiles=data.get('smiles') or None)
    species.characterize()
    apply_input_reference(species, data)
    return species


def apply_input_reference(species, parameters):
    """Restore the declared population even if its selected member is a mirror."""
    population = parameters.get('optical_population', 'specified')
    reference = parameters.get('stereo_reference')
    if reference is not None:
        observed = require_supported_identity(species)
        allowed = {observed.get('id')}
        if population == 'racemic':
            allowed.add(observed.get('mirror_id'))
        if reference.get('status') != 'assigned' or reference.get('id') not in allowed:
            from kinbot.stereo_routing import refuse_routing
            refuse_routing('Serialized stereo reference disagrees with the input geometry/population', [species])
        # Recompute all reference fields instead of trusting serialized hashes
        # or a supplied mirror_id independently of the actual input geometry.
        if reference['id'] != observed['id']:
            import numpy as np
            observed = require_supported_identity(species, np.asarray(species.geom) * [-1., 1., 1.])
        species.optical_reference = observed
    species.optical_population = population


def prepare_pes_directory(root, species):
    """Use only the current stereoisomer directory and calculation format."""
    from kinbot.run_format import ensure_current_run
    target = Path(root) / routing_name(species)
    ensure_current_run(target, create=True)
    return target


def expand_pes_names(root, choices):
    """Resolve connectivity selectors to the configured directories on disk."""
    root = Path(root)
    names = []
    for choice in choices:
        candidates = sorted(path for path in root.iterdir()
                            if path.is_dir() and is_species_name(path.name)
                            and matches_name(path.name, [str(choice)]))
        if not candidates:
            names.append(str(choice))  # Preserve the existing missing-well error.
        for path in candidates:
            from kinbot.run_format import ensure_current_run
            ensure_current_run(path)
            if path.name not in names:
                names.append(path.name)
    return names


def configured_result_identity(db, name, row, population=None):
    """Return the full declared identity for a verified configured PES result.

    Ordinary connectivity names are returned unchanged. None rejects a result
    whose geometry or saved full identity does not agree with its name.
    """
    match = _KEY.fullmatch(str(name))
    if match is None or match.group(2) is None:
        return str(name)
    from kinbot.stationary_pt import StationaryPoint
    from kinbot.stereo_routing import _row_species
    source = f'{name}_well'
    references = list(db.select(name=f'stereochemistry/{source}'))
    reference = references[-1] if references else None
    if reference is None:
        # The shortened filename cannot supply the missing full identity.
        return None
    scope = reference.data.get('chemical_context', {})
    requested = reference.data.get('identity', {})
    declared = scope.get('optical_reference') or requested
    if (declared.get('status') != 'assigned'
            or declared.get('id', '')[:16] != match.group(2)):
        return None
    allowed = {declared['id']}
    if (population or scope.get('optical_population')) == 'racemic':
        allowed.add(declared['mirror_id'])
    if requested.get('status') != 'assigned' or requested.get('id') not in allowed:
        return None
    template = StationaryPoint('energy input',
        reference.data['input_charge'],
        int(match.group(1)[-1]), atom=row.symbols, geom=row.positions)
    observed = _row_species(template, row, str(name), reference)
    identity = require_supported_identity(observed)
    return declared['id'] if identity['id'] in allowed else None


def configured_result_matches(db, name, row, population=None):
    """Compare a PES result with its complete saved identity, not its short name."""
    return configured_result_identity(db, name, row, population) is not None
