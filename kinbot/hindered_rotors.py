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

        # -1 (not finished), 0 (successful) or
        # 1 (failed) for each HIR scan point
        self.hir_status = []
        # energies of all the HIR scan points
        self.hir_energies = []
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
                            if ai == 0: 
                                success = 1
                            # cut off barriers above 20 kcal/mol to prevent the Fourier fit to oscillate
                            elif (energy - self.hir_energies[rotor][0]) < 20. / constants.AUtoKCAL:
                                success = 1
                            else:
                                success = -1
                        else:
                            success = -1
                if success == 1:
                    err, energy = self.qc.get_qc_energy(job)
                    self.hir_status[rotor][ai] = 0
                    self.hir_energies[rotor][ai] = energy
                    self.hir_geoms[rotor][ai] = geom
                elif success == -1:
                    logger.warning("Hindered rotor optimization not successful for {}".format(job))
                    self.hir_status[rotor][ai] = 1
                    self.hir_energies[rotor][ai] = -1
                    self.hir_geoms[rotor][ai] = geom

        return 0

    def check_hir(self, wait=0):
        """
        Check for hir calculations and optionally wait for them to finish
        """
        while 1:
            # check if all the calculations are finished
            self.test_hir()
            if len(self.species.dihed) == 0:
                logger.debug(f'No hindered rotors for {self.species.name}.')
            for rotor in range(len(self.species.dihed)):
                status = self.hir_status[rotor]
                if any([st < 0 for st in status]):
                    continue
                energies = self.hir_energies[rotor]
                if abs(energies[0] - self.species.energy) * constants.AUtoKCAL > 0.1:
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
                    return 0
                # energies taken if status = 0, successful geom check or normal gauss termination
                ens = [(energies[i] - energies[0]) * constants.AUtoKCAL 
                       for i in range(len(status)) if status[i] == 0]

            # if job finishes status set to 0 or 1, if all done then do the following calculation
            if all([all([test >= 0 for test in status]) for status in self.hir_status]):
                for rotor in range(len(self.species.dihed)):
                    if self.hir_status[rotor][0] == 2:  # skipped rotor
                        continue
                    if self.species.wellorts:
                        job = self.species.name + '_hir_' + str(rotor)
                    else:
                        job = str(self.species.chemid) + '_hir_' + str(rotor)
                    if len(ens) < self.nrotation - 2:
                        logger.warning("More than 2 HIR calculations failed for " + job)

                    angles = [i * 2 * np.pi / float(self.nrotation) for i in range(self.nrotation)]
                    # write profile to file
                    self.write_profile(rotor, job)
                    # Check to see if HIR failed, job will continue if failed, but warning will be generated
                    a = self.fourier_fit(job, angles, rotor)
                    if(a == 0):
                        logger.warning("FAILED HIR - empty energy array sent to fourier_fit for " + job)
                    else:
                        self.hir_fourier.append(self.fourier_fit(job, angles, rotor))
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
                s = str(self.species.natom) + '\n'
                s += 'energy = ' + str(self.hir_energies[rotor][i]) + '\n'
                for j, at in enumerate(self.species.atom):
                    x, y, z = self.hir_geoms[rotor][i][j]
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
            a = 0

        return a

    def get_fit_value(self, ai):
        """
        Get the fitted energy
        """
        e = 0.
        for j in range(self.n_terms):
            e += self.A[j] * (1 - np.cos((j+1) * ai))
            e += self.A[j+self.n_terms] * np.sin((j+1) * ai)
        return e
