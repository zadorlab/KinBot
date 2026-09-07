import time
import logging
import numpy as np
import matplotlib.pyplot as plt
from kinbot import constants
from kinbot import geometry
from kinbot import zmatrix
from kinbot.frequencies import skip_rotor
from kinbot.stationary_pt import StationaryPoint

logger = logging.getLogger('KinBot')


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
        # number of terms for Fourier
        self.n_terms = 6
        # all the geometries of the HIR scan points
        self.hir_geoms = []

    def generate_hir_geoms(self, cart, rigid):
        """
        Generate the initial geometries of the points along the scans
        """
        # re-initialize the lists in case of a restart of the HIR scans
        self.hir_status = []
        self.hir_energies = []
        self.hir_raw_energies = []
        self.hir_fourier = []
        self.hir_geoms = []

        while len(self.hir_status) < len(self.species.dihed):
            self.hir_status.append([-1 for i in range(self.nrotation)])
            self.hir_energies.append([-1 for i in range(self.nrotation)])
            self.hir_geoms.append([[] for i in range(self.nrotation)])

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
            self.qc.qc_hir(self.species, cart_new, rotor, 0, [fi], rigid)
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
                self.qc.qc_hir(self.species, cart_new, rotor, ai, [fi], rigid)
        return 0

    def test_hir(self):
        for rotor in range(len(self.species.dihed)):
            for ai in range(self.nrotation):
                if ai and self.hir_status[rotor][0] == -1:
                    # A completed point cannot be screened against an unknown
                    # reference energy. Revisit it after point zero finishes.
                    continue
                success = None
                if self.hir_status[rotor][ai] == -1:
                    if self.species.wellorts:
                        job = 'hir/' + self.species.name + '_hir_' + str(rotor) + '_' + str(ai).zfill(2)
                    else:
                        job = 'hir/' + str(self.species.chemid) + '_hir_' + str(rotor) + '_' + str(ai).zfill(2)
                    err, geom = self.qc.get_qc_geom(job, self.species.natom)
                    if err == 1:  # still running
                        continue
                    elif err == -1:  # failed
                        success = -1
                    else:
                        # check if all the bond lenghts are within
                        # 15% of the original bond lengths
                        temp = StationaryPoint('temp',
                                               self.species.charge,
                                               self.species.mult,
                                               atom=self.species.atom,
                                               geom=geom)
                        temp.characterize()
                        if geometry.equal_geom(self.species,
                                               temp,
                                               0.15):
                            err, energy = self.qc.get_qc_energy(job)
                            if err or not np.isfinite(energy):
                                success = -1
                            elif ai == 0:
                                success = 1
                            # cut off barriers above 20 kcal/mol to prevent the Fourier fit to oscillate
                            elif (energy - self.hir_energies[rotor][0]) < 20. / constants.AUtoKCAL:
                                success = 1
                            else:
                                success = -1
                        else:
                            success = -1
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
        the existing Fourier-fill policy and do not make a rotor invalid.
        """
        if rotor >= len(self.hir_status) or len(self.hir_status[rotor]) != self.nrotation:
            return 'no scan recorded'
        status = self.hir_status[rotor]
        if status[0] == 2:
            return 'scan skipped'
        if status[0] == 1:
            return 'reference point failed or scans disabled by the rotor-0 energy test'
        if any(value < 0 for value in status):
            return 'scan incomplete'
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
                    self.hir_status = [[1 for ai in ri] for ri in self.hir_status]

            # if job finishes status set to 0 or 1, if all done then do the following calculation
            if all([all([test >= 0 for test in status]) for status in self.hir_status]):
                for rotor in range(len(self.species.dihed)):
                    if self.hir_status[rotor][0] == 2 or self.hir_status[rotor][0] == 1:  # skipped or corrupted rotor
                        continue
                    if self.species.wellorts:
                        job = self.species.name + '_hir_' + str(rotor)
                    else:
                        job = str(self.species.chemid) + '_hir_' + str(rotor)
                    if self.hir_status[rotor].count(0) < self.nrotation - 2:
                        logger.warning("More than 2 HIR calculations failed for " + job)

                    angles = [i * 2 * np.pi / float(self.nrotation) for i in range(self.nrotation)]
                    # write profile to file
                    self.write_profile(rotor, job)
                    # Check to see if HIR failed, job will continue if failed, but warning will be generated
                    a = self.fourier_fit(job, angles, rotor)
                    if(a == 0):
                        logger.warning("FAILED HIR - empty energy array sent to fourier_fit for " + job)
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
