from kinbot.species_routing import routing_name
import time
import copy
import logging
import numpy as np
import matplotlib.pyplot as plt
from kinbot import constants
from kinbot import geometry
from kinbot import zmatrix
from kinbot.frequencies import skip_rotor
from kinbot.stationary_pt import StationaryPoint
from kinbot.calculation import geometry_reference
from kinbot.reaction_path import path_geometry_allowed

logger = logging.getLogger('KinBot')


def use_harmonic_model(species, par, reason):
    """Omit HIR consistently while retaining its calculated observations."""
    from kinbot.frequencies import thermochemical_frequencies
    logger.warning('%s: %s; retaining full harmonic frequencies and omitting '
                   'hindered rotors.', species.name, reason)
    hir = getattr(species, 'hir', None)
    if hir is not None:
        hir.projection_failure = reason
        hir._harmonic_fallback_reference = (
            geometry_reference(species), copy.deepcopy(species.dihed))
    modes = thermochemical_frequencies(species.freq, species.wellorts,
                                      par.get('imagfreq_threshold', 50.))
    species.kinbot_freqs = list(modes)
    species.reduced_freqs = list(modes)
    reference = copy.deepcopy((getattr(species, 'rotor_projection', None) or {}).get('reference'))
    species.rotor_projection = dict(method='harmonic_fallback', internal_rank=0,
        reference=reference,
        reason=reason, rotors=[dict(rotor_index=i, projected=False, reason=reason)
                             for i in range(len(species.dihed))])


def recover_hir_model(species, qc, par):
    """Attempt existing scan recovery, then verify the complete thermal model."""
    hir = getattr(species, 'hir', None)
    if hir is None or par.get('multi_conf_tst') or not par.get('rotor_scan'):
        return False
    reference = (geometry_reference(species), species.dihed)
    if getattr(hir, '_harmonic_fallback_reference', None) == reference:
        return True  # Do not resubmit a failed recovery for the same input.
    try:
        changed = _recover_hir_model(species, qc, par)
        if changed and species.dihed and not any(hir.is_valid_rotor(i)
                                                for i in range(len(species.dihed))):
            use_harmonic_model(species, par, 'No usable HIR scans remain after recovery')
            return True
        from kinbot.counting_contract import optical_counting
        from kinbot.thermochemistry import hir_evidence
        counting = optical_counting(species, hir_evidence(species))
        if (counting['status'] == 'unresolved'
                and counting.get('remaining_multiplier') is None):
            use_harmonic_model(species, par, 'HIR recovery: ' + counting['reason'])
            return True
        return changed
    except (ValueError, OSError, RuntimeError, AttributeError, IndexError) as error:
        use_harmonic_model(species, par, f'HIR recovery failed: {error}')
        return True


def _recover_hir_model(species, qc, par):
    """Repair stale scan bookkeeping or replace scans, then restore projection.

    Existing calculation files are never discarded. A replacement uses the
    requested input's stable job name, so a failed replacement is a terminal
    failed scan rather than another reason to submit the same work again.
    """
    from kinbot import frequencies
    from kinbot.calculation import array_fingerprint

    hir = getattr(species, 'hir', None)
    if hir is None or par.get('multi_conf_tst') or not par.get('rotor_scan'):
        return False
    if not getattr(species, 'source_job', None) or getattr(species, 'source_row_id', None) is None:
        source = getattr(species, 'source_job', None) or getattr(species, 'hessian_source_job', None)
        rows = list(qc.db.select(name=source)) if source and getattr(qc, 'db', None) is not None else []
        if rows:
            row = rows[-1]
            energy = row.data.get('energy')
            if (isinstance(energy, (int, float, np.number)) and np.isfinite(energy)
                    and np.allclose(row.positions, species.geom, atol=1.e-6, rtol=0.)
                    and np.isclose(energy * constants.EVtoHARTREE,
                                   species.energy, atol=1.e-10, rtol=0.)
                    and np.array_equal(row.data.get('frequencies'), species.freq)):
                species.source_job, species.source_row_id = row.name, row.id
    reference = geometry_reference(species)
    reference.update(backend=getattr(qc, 'qc', None),
                     dihedrals=[list(map(int, rotor)) for rotor in species.dihed])
    previous = getattr(hir, 'scan_reference', None) or {}
    keys = ('geometry_sha256', 'atoms', 'dihedrals', 'source_job', 'source_row_id')
    scan_changed = any(previous.get(key) != reference.get(key) for key in keys)
    from kinbot.thermochemistry import hir_evidence
    observations_changed = any(
        point['status'] == 'successful' and (point['geometry_angstrom'] is None
                                            or point['electronic_energy_hartree'] is None)
        for rotor in hir_evidence(species)['rotors'] if rotor['usable']
        for point in rotor['points'])
    scan_changed |= observations_changed
    projection = getattr(species, 'rotor_projection', None) or {}
    projected_reference = projection.get('reference', {})
    projection_changed = any(projected_reference.get(key) != reference.get(key) for key in keys)
    entries = projection.get('rotors', [])
    flags = [entry.get('projected') for entry in entries]
    indices = [entry.get('rotor_index') for entry in entries]
    valid_indices = (all(isinstance(i, (int, np.integer)) and not isinstance(i, (bool, np.bool_))
                         for i in indices)
                     and sorted(indices) == list(range(len(species.dihed))))
    hessian = np.asarray(getattr(species, 'hess', []))
    projection_changed |= (not valid_indices
        or any(not isinstance(flag, (bool, np.bool_)) for flag in flags)
        or projection.get('internal_rank') != sum(bool(flag) for flag in flags)
        or len(species.freq) - len(species.reduced_freqs) != projection.get('internal_rank')
        or (valid_indices and any(entry.get('projected') and not hir.is_valid_rotor(entry['rotor_index'])
                                  for entry in entries))
        or bool(hessian.size and projected_reference.get('hessian_sha256') != array_fingerprint(hessian))
        or not reference['source_job'] or reference['source_row_id'] is None)
    if not scan_changed and not projection_changed:
        return False
    hir.species, hir.qc = species, qc
    if scan_changed:
        # A saved complete reference with identical Cartesian input proves a
        # mere source-name/row change. Missing or different geometry evidence
        # must instead be checked against the literal initial scan inputs.
        trusted = not observations_changed and all(previous.get(key) == reference.get(key)
                      for key in ('geometry_sha256', 'atoms', 'dihedrals'))
        if not trusted and not observations_changed:
            try:
                for rotor, definition in enumerate(species.dihed):
                    qc._check_hir_definition(hir.point_job(rotor, 0),
                        [int(atom) + 1 for atom in definition], True,
                        geometry=species.geom, atoms=species.atom)
                trusted = True
            except (ValueError, OSError, AttributeError):
                trusted = False
        if trusted:
            logger.warning('%s: recovered HIR reference from matching saved input.', species.name)
            hir.scan_reference = reference
        else:
            attempt = (reference['geometry_sha256'], repr(reference['dihedrals']))
            if getattr(hir, '_recovery_attempt', None) == attempt:
                logger.warning('%s: replacement HIR input remains inconsistent; retaining harmonic modes.', species.name)
                hir.hir_status = [[1] * hir.nrotation for _ in species.dihed]
                hir.scan_reference = reference
            else:
                hir._recovery_attempt = attempt
                history = getattr(hir, 'recovery_observations', [])
                history.append({key: copy.deepcopy(getattr(hir, key, None)) for key in
                    ('scan_reference', 'scan_jobs', 'hir_status', 'hir_raw_energies',
                     'hir_energies', 'hir_geoms', 'point_observations')})
                hir.recovery_observations = history
                logger.warning('%s: HIR inputs changed or lack supporting records; '
                               'recomputing affected scans and preserving previous results.', species.name)
                try:
                    hir.generate_hir_geoms(np.asarray(species.geom).copy(), par.get('rigid_hir', False))
                    hir.check_hir(wait=1)
                except (ValueError, OSError, RuntimeError) as error:
                    logger.warning('%s: replacement HIR failed (%s); retaining harmonic modes.',
                                   species.name, error)
                    hir.hir_status = [[1] * hir.nrotation for _ in species.dihed]
                    hir.scan_reference = reference
    # Always start with the selected Cartesian Hessian, never the previously
    # reduced frequencies. Failed replacement rotors therefore keep their mode.
    hir.projection_status = []
    try:
        hessian = np.asarray(getattr(species, 'hess', []), dtype=float)
        shape = (3 * species.natom, 3 * species.natom)
        source = getattr(species, 'source_job', None)
        if not source or getattr(species, 'source_row_id', None) is None:
            raise ValueError('selected calculation source is not established')
        if (hessian.shape != shape or not np.all(np.isfinite(hessian))
                or getattr(species, 'hessian_source_job', None) != source):
            if not source:
                raise ValueError('selected Hessian source is unknown')
            hessian = np.asarray(qc.read_qc_hess(source, species.natom), dtype=float)
        if hessian.shape != shape or not np.all(np.isfinite(hessian)):
            raise ValueError('selected Hessian is unavailable')
        species.hess = hessian.tolist()
        species.hessian_source_job = source
        hir.__dict__.pop('projection_failure', None)
        species.kinbot_freqs, species.reduced_freqs = frequencies.get_frequencies(
            species, hessian, np.asarray(species.geom), massweighted=qc.hessian_is_massweighted())
    except (ValueError, OSError, RuntimeError, AttributeError) as error:
        reason = f'HIR recovery has no usable selected Hessian: {error}'
        use_harmonic_model(species, par, reason)
    for attr in ('kinbot_freqs', 'reduced_freqs'):
        setattr(species, attr, frequencies.thermochemical_frequencies(
            getattr(species, attr), species.wellorts, par.get('imagfreq_threshold', 50.)))
    return True


class HIR:
    """
    Class that does all the steps for the HIR calculations of one species
    """
    def __init__(self, species, qc, par):
        """
        species: instance of StationaryPoint
        qc: instance of QuantumChemistry
        par: instance of Parameters
        """
        self.species = species
        self.qc = qc

        # number of points along one scan
        self.nrotation = par['nrotation']
        # boolean tells if profiles should be plotted
        self.plot_hir_profiles = par['plot_hir_profiles']
        # if 1 (default), remove bad rotors
        self.rotor_0_test = par['rotor_0_test']

        # -1 (not finished), 0 (successful), 1 (failed), or 2 (skipped)
        # for each HIR scan point
        self.hir_status = []
        # energies of all the HIR scan points
        self.hir_energies = []
        # Results before failed points are filled from a Fourier fit.
        self.hir_raw_energies = []
        # Fourier fit of each scan
        self.hir_fourier = []
        self.hir_fit_diagnostics = []
        self.projection_status = []
        self.scan_jobs = []
        self.rigid_scan = None
        # number of terms for Fourier
        self.n_terms = 6
        # all the geometries of the HIR scan points
        self.hir_geoms = []
        self.scan_reference = None
        self.point_observations = []
        self.demotion_reason = None
        self.statuses_before_demotion = None

    def generate_hir_geoms(self, cart, rigid):
        """
        Generate the initial geometries of the points along the scans
        """
        # re-initialize the lists in case of a restart of the HIR scans
        self.hir_status = []
        self.hir_energies = []
        self.hir_raw_energies = []
        self.hir_fourier = []
        self.hir_fit_diagnostics = []
        self.projection_status = []
        self.scan_jobs = []
        self.rigid_scan = bool(rigid)
        self.hir_geoms = []
        self.scan_reference = geometry_reference(self.species, cart)
        self.scan_reference['backend'] = getattr(self.qc, 'qc', None)
        self.scan_reference['dihedrals'] = [list(map(int, rotor)) for rotor in self.species.dihed]
        self.point_observations = []
        self.demotion_reason = None
        self.statuses_before_demotion = None

        while len(self.hir_status) < len(self.species.dihed):
            self.hir_status.append([-1 for i in range(self.nrotation)])
            self.scan_jobs.append([None for i in range(self.nrotation)])
            self.hir_energies.append([-1 for i in range(self.nrotation)])
            self.hir_geoms.append([[] for i in range(self.nrotation)])
            self.point_observations.append([None for i in range(self.nrotation)])

        for rotor in range(len(self.species.dihed)):
            if skip_rotor(self.species.name, self.species.dihed[rotor]) == 1:
                self.hir_status[rotor] = [2 for i in range(self.nrotation)]
                logger.info('\tFor {} rotor {} was skipped in HIR.'.format(self.species.name, rotor))
                continue

            cart = np.asarray(cart)
            zmat_atom, zmat_ref, zmat, zmatorder = zmatrix.make_zmat_from_cart(self.species, rotor, cart, 0)

            # first element has same geometry
            cart_new = zmatrix.make_cart_from_zmat(zmat,
                                                   zmat_atom,
                                                   zmat_ref,
                                                   self.species.natom,
                                                   self.species.atom,
                                                   zmatorder)
            fi = [(zi + 1) for zi in zmatorder[:4]]
            self._submit_point(cart_new, rotor, 0, fi, rigid)
            for ai in range(1, self.nrotation):
                ang = 360. / float(self.nrotation)
                zmat[3][2] += ang
                for i in range(4, self.species.natom):
                    if zmat_ref[i][2] == 4:
                        zmat[i][2] += ang
                    if zmat_ref[i][2] == 1:
                        zmat[i][2] += ang
                cart_new = zmatrix.make_cart_from_zmat(zmat,
                                                       zmat_atom,
                                                       zmat_ref,
                                                       self.species.natom,
                                                       self.species.atom,
                                                       zmatorder)
                self._submit_point(cart_new, rotor, ai, fi, rigid)
        return 0

    def _submit_point(self, geom, rotor, point, dihedral, rigid):
        try:
            self.scan_jobs[rotor][point] = self.qc.qc_hir(
                self.species, geom, rotor, point, [dihedral], rigid)
        except (ValueError, OSError, RuntimeError) as error:
            logger.warning('%s rotor %s point %s could not be submitted: %s',
                           self.species.name, rotor, point, error)
            self.hir_status[rotor][point] = 1
            self.point_observations[rotor][point] = {
                'accepted': False, 'rejection_reason': str(error),
                'qc_geometry_status': 'not_requested', 'qc_energy_status': 'not_requested',
                'geometry_angstrom': None, 'electronic_energy_hartree': None,
            }

    def point_job(self, rotor, point):
        """Use the calculation actually submitted, including replacement scans."""
        jobs = getattr(self, 'scan_jobs', [])
        if rotor < len(jobs) and point < len(jobs[rotor]) and isinstance(jobs[rotor][point], str):
            return jobs[rotor][point]
        name = self.species.name if self.species.wellorts else routing_name(self.species)
        return f'hir/{name}_hir_{rotor}_{point:02d}'

    def test_hir(self):
        for rotor in range(len(self.species.dihed)):
            for ai in range(self.nrotation):
                if ai and self.hir_status[rotor][0] == -1:
                    # A completed point cannot be screened against an unknown
                    # reference energy. Revisit it after point zero finishes.
                    continue
                success = None
                if self.hir_status[rotor][ai] == -1:
                    energy, energy_status, reason = None, None, None
                    job = self.point_job(rotor, ai)
                    err, geom = self.qc.get_qc_geom(job, self.species.natom)
                    geometry_status = err
                    if err == 1:  # still running
                        continue
                    elif err == -1:  # failed
                        success = -1
                        reason = 'QC geometry failed'
                    else:
                        # check if all the bond lenghts are within
                        # 15% of the original bond lengths
                        temp = StationaryPoint('temp',
                                               self.species.charge,
                                               self.species.mult,
                                               atom=self.species.atom,
                                               geom=geom)
                        temp.characterize()
                        if not path_geometry_allowed(self.species, geom):
                            success = -1
                            reason = 'scan crosses into a different stereochemical pathway'
                        elif geometry.equal_geom(self.species,
                                               temp,
                                               0.15):
                            err, energy = self.qc.get_qc_energy(job)
                            energy_status = err
                            if err or not np.isfinite(energy):
                                success = -1
                                reason = 'QC energy unavailable or nonfinite'
                            elif ai == 0:
                                success = 1
                            # cut off barriers above 20 kcal/mol to prevent the Fourier fit to oscillate
                            elif (energy - self.hir_energies[rotor][0]) < 20. / constants.AUtoKCAL:
                                success = 1
                            else:
                                success = -1
                                reason = 'scan point exceeds the 20 kcal/mol acceptance limit'
                        else:
                            success = -1
                            reason = 'geometry failed the existing bond-length check'
                    # Acceptance arrays use failure sentinels. Keep the actual
                    # available observation before those arrays are overwritten.
                    while len(self.point_observations) <= rotor:
                        self.point_observations.append([None] * self.nrotation)
                    coordinates = np.asarray(geom)
                    self.point_observations[rotor][ai] = {
                        'source_job': job,
                        'qc_geometry_status': 'failed' if geometry_status == -1 else 'returned',
                        'qc_energy_status': {None: 'not_requested', -1: 'failed', 1: 'running',
                                             0: 'returned'}.get(energy_status, 'unknown'),
                        'electronic_energy_hartree': (float(energy) if energy_status == 0 and energy is not None
                                                     and np.isfinite(energy) else None),
                        'geometry_angstrom': (coordinates.tolist()
                                              if geometry_status == 0 and coordinates.shape == (self.species.natom, 3)
                                              and np.all(np.isfinite(coordinates)) else None),
                        'accepted': success == 1, 'rejection_reason': reason,
                    }
                if success == 1:
                    self.hir_status[rotor][ai] = 0
                    self.hir_energies[rotor][ai] = energy
                    self.hir_geoms[rotor][ai] = geom
                elif success == -1:
                    logger.warning("Hindered rotor optimization not successful for {}".format(job))
                    self.hir_status[rotor][ai] = 1
                    self.hir_energies[rotor][ai] = -1
                    self.hir_geoms[rotor][ai] = geom

        return 0

    def invalid_rotor_reason(self, rotor):
        """Why a rotor is excluded from the hindered-rotor treatment.

        Returns None for a usable scan. Failed non-reference points retain
        master's Fourier-fill policy; stereochemical pathway guards still apply.
        """
        if getattr(self, 'projection_failure', None):
            return self.projection_failure
        observations = getattr(self, 'point_observations', [])
        if rotor < len(observations) and any(
                point and point.get('rejection_reason') ==
                'scan crosses into a different stereochemical pathway'
                for point in observations[rotor]):
            # Interpolating across another explicitly counted route would
            # overlap its contribution. Retain this rotor's harmonic mode.
            return 'scan crosses into a different stereochemical pathway'
        if rotor >= len(self.hir_status) or len(self.hir_status[rotor]) != self.nrotation:
            return 'no scan recorded'
        status = self.hir_status[rotor]
        if status[0] == 2:
            return 'scan skipped'
        if status[0] == 1:
            return 'reference point failed or scans disabled by the rotor-0 energy test'
        if any(value < 0 for value in status):
            return 'scan incomplete'
        if any(value == 2 for value in status):
            return 'scan contains skipped points'
        if rotor >= len(self.hir_energies) or len(self.hir_energies[rotor]) != self.nrotation:
            return 'no scan energies recorded'
        if any(not np.isfinite(self.hir_energies[rotor][i])
               for i, value in enumerate(status) if value == 0):
            return 'nonfinite successful scan energy'
        if 1 in status and (rotor >= len(self.hir_fourier)
                            or self.hir_fourier[rotor] is None):
            return 'partial scan has no usable Fourier interpolation'
        if rotor < len(self.projection_status):
            if not self.projection_status[rotor]['projected']:
                return self.projection_status[rotor]['reason']
        return None

    def is_valid_rotor(self, rotor):
        """Whether a completed scan has a usable reference point.

        Completion alone does not enable a failed or skipped rotor; such
        rotors stay harmonic oscillators in the frequency set and in MESS.
        """
        return self.invalid_rotor_reason(rotor) is None

    def demoted_rotor_summary(self):
        """One-line account of the rotors MESS will treat as harmonic, or ''."""
        demoted = []
        for rotor, rot in enumerate(self.species.dihed):
            why = self.invalid_rotor_reason(rotor)
            if why is not None:
                demoted.append(f'rotor {rotor} (atoms {rot[1] + 1}-{rot[2] + 1}): {why}')
        if not demoted:
            return ''
        return (f'{len(demoted)} of {len(self.species.dihed)} rotors of '
                f'{self.species.name} will be treated as harmonic oscillators '
                f'in MESS: ' + '; '.join(demoted))

    def check_hir(self, wait=0):
        """
        Return completion, including failed/skipped scans. Use is_valid_rotor
        to decide which completed scans can replace harmonic modes.
        """
        while 1:
            # check if all the calculations are finished
            self.test_hir()
            if len(self.species.dihed) == 0:
                logger.debug(f'No hindered rotors for {self.species.name}.')
            for rotor in range(len(self.species.dihed)):
                status = self.hir_status[rotor]
                if any([st < 0 for st in status]):  # at least one is running
                    continue
                if self.hir_status[rotor][0] in (1, 2):  # failed or skipped
                    continue
                energies = self.hir_energies[rotor]
                if abs(energies[0] - self.species.energy) * constants.AUtoKCAL > 0.1 and self.rotor_0_test:
                    logger.warning(f'\t0 angle rotor for rotor {rotor} has a different energy than '
                                   'the optimized structure for '
                                   f'{self.species.name} ({energies[0]} vs {self.species.energy}).')
                    logger.warning('This might be '
                                   'caused by an SCF convergence issue. '
                                   'Hindered rotors are disabled for this '
                                   'stationary point.')
                    logger.warning(rotor)
                    logger.warning(energies)
                    self.demotion_reason = 'rotor-zero energy differs from the optimized reference'
                    self.statuses_before_demotion = [list(row) for row in self.hir_status]
                    self.hir_status = [[1 for ai in ri] for ri in self.hir_status]

            # if job finishes status set to 0 or 1, if all done then do the following calculation
            if all([all([test >= 0 for test in status]) for status in self.hir_status]):
                for rotor in range(len(self.species.dihed)):
                    if self.hir_status[rotor][0] == 2 or self.hir_status[rotor][0] == 1:  # skipped or corrupted rotor
                        continue
                    if self.species.wellorts:
                        job = self.species.name + '_hir_' + str(rotor)
                    else:
                        job = routing_name(self.species) + '_hir_' + str(rotor)
                    replacements = [point for point in self.scan_jobs[rotor]
                                    if isinstance(point, str) and '_recovery_' in point] \
                                   if rotor < len(self.scan_jobs) else []
                    if replacements:
                        job += '_recovery_' + replacements[0].split('_recovery_')[-1]
                    if self.hir_status[rotor].count(0) < self.nrotation - 2:
                        logger.warning("More than 2 HIR calculations failed for " + job)

                    angles = [i * 2 * np.pi / float(self.nrotation) for i in range(self.nrotation)]
                    # write profile to file
                    self.write_profile(rotor, job)
                    # Check to see if HIR failed, job will continue if failed, but warning will be generated
                    a = self.fourier_fit(job, angles, rotor)
                    if a == 0:
                        logger.warning('HIR fit unavailable for %s: %s', job,
                                       self.invalid_rotor_reason(rotor))
                summary = self.demoted_rotor_summary()
                if summary:
                    logger.warning('\t' + summary)
                return 1
            else:
                if wait:
                    time.sleep(1)
                else:
                    return 0

    def write_profile(self, rotor, job):
        """
        Write a molden-readable file with the
        HIR scan (geometries and energies)
        """
        with open('hir/' + job + '.xyz', 'w') as ff:
            for i in range(self.nrotation):
                geom = np.asarray(self.hir_geoms[rotor][i])
                if geom.shape != (self.species.natom, 3):
                    # Failed calculations may have no geometry to display.
                    continue
                s = str(self.species.natom) + '\n'
                s += 'energy = ' + str(self.hir_energies[rotor][i]) + '\n'
                for j, at in enumerate(self.species.atom):
                    x, y, z = geom[j]
                    s += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
                ff.write(s)
        return

    def fourier_fit(self, job, angles, rotor):
        """
        Create a alternative fourier formulation of a hindered rotor
        profile, the angles are in radians and the energies in
        kcal per mol (Vanspeybroeck et al.)
        """
        energies = self.hir_energies[rotor]
        status = self.hir_status[rotor]
        while len(self.hir_fourier) <= rotor:
            self.hir_fourier.append(None)
        # Retain the optional export slot, without claiming the legacy fit
        # has passed the removed rank/order/coverage diagnostics.
        while len(self.hir_fit_diagnostics) <= rotor:
            self.hir_fit_diagnostics.append(None)
        self.hir_fit_diagnostics[rotor] = None
        while len(self.hir_raw_energies) <= rotor:
            self.hir_raw_energies.append(None)
        if self.hir_raw_energies[rotor] is None:
            self.hir_raw_energies[rotor] = list(energies)

        ang = [angles[i] for i in range(len(status)) if status[i] == 0]
        ens = [(energies[i] - energies[0])*constants.AUtoKCAL for i in range(len(status)) if status[i] == 0]

        X = np.zeros((len(ang), 2 * self.n_terms))
        for i, ai in enumerate(ang):
            for j in range(self.n_terms):
                X[i][j] = (1 - np.cos((j+1) * ai))
                X[i][j+self.n_terms] = np.sin((j+1) * ai)

        if(len(ens) > 0):
            a = 1
            self.A = np.linalg.lstsq(X, np.array(ens), rcond=None)[0]
            self.hir_fourier[rotor] = self.A.copy().tolist()

            for i, si in enumerate(status):
                if si == 1:
                    energies[i] = energies[0] + self.get_fit_value(angles[i])/constants.AUtoKCAL
            if self.plot_hir_profiles:
                # fit the plot to a png file
                plt.plot(ang, ens, 'ro')
                fit_angles = [i * 2. * np.pi / 360 for i in range(360)]
                fit_energies = [self.get_fit_value(ai) for ai in fit_angles]
                plt.plot(fit_angles, fit_energies)
                plt.xlabel('Dihedral angle [radians]')
                plt.ylabel('Energy [kcal/mol]')
                plt.savefig('hir_profiles/{}.png'.format(job))
                plt.clf()
        else:
            self.A = 0
            self.hir_fourier[rotor] = None
            a = 0

        return a

    def get_fit_value(self, ai, rotor=None):
        """
        Get the fitted energy
        """
        coefficients = self.A if rotor is None else self.hir_fourier[rotor]
        if coefficients is None:
            raise ValueError(f'No Fourier fit is available for rotor {rotor}.')
        e = 0.
        for j in range(self.n_terms):
            e += coefficients[j] * (1 - np.cos((j+1) * ai))
            e += coefficients[j+self.n_terms] * np.sin((j+1) * ai)
        return e
