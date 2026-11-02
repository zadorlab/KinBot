"""Mirror action on a native first-order dividing surface, separate from sites."""
from pathlib import Path
import hashlib
import json
import base64
import re
import numpy as np
from kinbot.calculation import array_fingerprint
from kinbot.molecular_symmetry import spatial_mappings


SURFACE_MARKER = re.compile(r'^! kinbot_transition_symmetry (.+)$', re.M)


def _documented_reactive_mode(species):
    """Read only the selected, atom-bound Hessian in documented native units."""
    from kinbot import constants, frequencies, geometry
    g = np.asarray(species.geom, float)
    h = np.asarray(getattr(species, 'hess', []), float)
    ref = getattr(species, 'optical_hessian_reference', None) or {}
    weighted = ref.get('hessian_massweighted')
    unit = ('hartree / (bohr^2 * amu)' if weighted is True else
            'hartree / bohr^2' if weighted is False else None)
    if not (h.shape == (g.size, g.size) and np.isfinite(h).all()
            and ref.get('geometry_sha256') == array_fingerprint(g)
            and ref.get('hessian_sha256') == array_fingerprint(h)
            and ref.get('atoms') == list(map(str, species.atom))
            and unit is not None and ref.get('hessian_unit') == unit
            and all(getattr(species, k, None) is not None
                    and ref.get(k) == getattr(species, k)
                    for k in ('source_job', 'source_row_id'))):
        return None
    mass = np.repeat([constants.exact_mass[str(a)] for a in species.atom], 3)
    hm = h if weighted else h / np.sqrt(np.outer(mass, mass))
    # Project rigid motions before counting negative vibrational modes. Tiny
    # translational/rotational eigenvalues do not define another saddle order.
    centered = g - geometry.get_center_of_mass(g, species.atom)
    translations, rotations = frequencies.rigid_body_vectors(centered, species.atom)
    external = np.asarray([v/np.linalg.norm(v) for v in (*translations, *rotations)
                           if np.linalg.norm(v) > 1.e-10])
    _, singular, vectors = np.linalg.svd(external, full_matrices=True)
    rank = int(np.sum(singular > 1.e-8*singular[0]))
    vibration = vectors[rank:]
    values, modes = np.linalg.eigh(vibration @ ((hm+hm.T)/2) @ vibration.T)
    negative = [i for i, v in enumerate(values)
                if frequencies.convert_to_wavenumbers(v) < -1.]
    if len(negative) != 1:
        return None
    mode = ((modes[:, negative[0]] @ vibration)/np.sqrt(mass)).reshape(g.shape)
    return mode, dict(kind='selected_hessian', **ref)


def mirror_mode_parity(species, mode, tolerance=.01):
    """Fit one whole geometry and transform its reactive displacement with it.

    An even unstable coordinate is preserved by the mirror; an odd one
    reverses reaction direction. Neither rotational sigma nor an endpoint's
    optical factor determines this property. The 0.99 parity separation
    tolerates the two-decimal native Gaussian displacement table.
    """
    g = np.asarray(species.geom, float)
    q = np.asarray(mode, float)
    if q.shape != g.shape or not np.isfinite(q).all() or not np.linalg.norm(q):
        raise ValueError('A finite nonzero reactive displacement is required.')
    mirror = g * [-1., 1., 1.]
    from kinbot import constants
    from kinbot.spatial import proper_alignment, mode_overlap
    masses = [constants.exact_mass[str(a)] for a in species.atom]
    results = []
    for order in spatial_mappings(species, g, mirror, tolerance):
        fit = proper_alignment(mirror[order], g)
        error = float(np.max(fit.errors))
        if error > tolerance:
            continue
        mapped = (q * [-1., 1., 1.])[order] @ fit.rotation
        parity = mode_overlap(q, mapped, masses=masses, signed=True)
        results.append(dict(atom_mapping=order, maximum_atom_error_angstrom=error,
                            reactive_mode_overlap=parity))
    if not results:
        return dict(status='distinct', reason='The reflected rigid surface is distinct.')
    parities = {'even' if r['reactive_mode_overlap'] >= .99 else
                'odd' if r['reactive_mode_overlap'] <= -.99 else 'mixed' for r in results}
    if 'mixed' in parities or len(parities) != 1:
        return dict(status='unresolved', mappings=results,
                    reason='Reactive mode has no unique mirror parity at native precision.')
    return dict(status='self', parity=parities.pop(), mappings=results,
                geometry_sha256=array_fingerprint(g),
                mode_sha256=array_fingerprint(q), method='native-reactive-mode-mirror-v1')


def selected_surface_symmetry(species):
    """Use an existing documented Hessian or matching native frequency output."""
    g = np.asarray(species.geom, float)
    source = None
    if any(getattr(species, 'isotopes', ())):
        return dict(status='unavailable', reason='Isotope-specific native masses are not bound to this mode.')
    documented = _documented_reactive_mode(species)
    if documented is not None:
        mode, source = documented
    if source is None:
        from kinbot.gaussian_modes import read_normal_modes
        for job in dict.fromkeys((getattr(species, 'hessian_source_job', None),
                                  getattr(species, 'source_job', None))):
            if not job or not Path(str(job)+'.log').is_file():
                continue
            try:
                native = read_normal_modes(str(job)+'.log', species.atom, g)
            except ValueError:
                continue
            negative = [i for i, v in enumerate(native['frequencies_cm1']) if v < -1.]
            if len(negative) == 1:
                mode = native['modes'][negative[0]]
                source = {k: v for k, v in native.items() if k != 'modes'}
                break
    if source is None:
        return dict(status='unavailable', reason='No matching native reactive mode was retained.')
    result = mirror_mode_parity(species, mode)
    result['source'] = source
    return result


def annotate_surface_symmetry(block, species):
    if SURFACE_MARKER.search(block):
        return block
    result = selected_surface_symmetry(species)
    if result.get('status') == 'self':
        # Reflected copies carry this same source ID; it identifies one native
        # dividing surface even when its endpoint populations are different.
        result['surface_id'] = hashlib.sha256(json.dumps(result, sort_keys=True).encode()).hexdigest()
        payload = base64.urlsafe_b64encode(json.dumps(result, sort_keys=True).encode()).decode()
        starts = list(re.finditer(r'^\s*RRHO[^\S\n]*\n', block, re.M))
        if len(starts) != 1:
            raise ValueError('Native dividing-surface evidence must be attached to its own RRHO member.')
        end = starts[0].end()
        block = block[:end] + '! kinbot_transition_symmetry ' + payload + '\n' + block[end:]
    return block


def annotate_cached_surface(block, worker, reaction, par):
    """Upgrade a legacy single-TS template in memory, without altering its files.

    Bind the native log to the exact serialized geometry and selected job
    level. Do not guess an MC member's calculation from its parent name.
    Newly written MC templates already carry member-specific evidence.
    """
    from kinbot.mess_symmetry import upgrade_cached_rigid_rotor
    block = upgrade_cached_rigid_rotor(block, par['charge'], par['mult'], saddle=True)
    from kinbot.mess_racemates import _RRHO
    if SURFACE_MARKER.search(block) or len(_RRHO.findall(block)) != 1:
        return block
    if re.search(r'Core\s+(?:PhaseSpaceTheory|Rotd)\b|^\s*Variational\b', block, re.M):
        return block
    job = Path(worker)/(reaction+('_high' if par['high_level'] else ''))
    if not Path(str(job)+'.log').is_file():
        return block
    geometry = re.search(r'Geometry\[angstrom\]\s+(\d+)\s*\n', block)
    if geometry is None:
        return block
    rows = block[geometry.end():].splitlines()[:int(geometry[1])]
    atoms, positions = [], []
    try:
        for row in rows:
            atom, *xyz = row.split()
            if len(xyz) != 3:
                return block
            atoms.append(atom); positions.append(list(map(float, xyz)))
    except ValueError:
        return block
    from kinbot.stationary_pt import StationaryPoint
    point = StationaryPoint(reaction, par['charge'], par['mult'],
                            atom=atoms, geom=np.asarray(positions), wellorts=1)
    point.characterize()
    point.source_job = str(job)
    return annotate_surface_symmetry(block, point)
