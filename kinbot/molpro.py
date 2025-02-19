import os
import logging
import numpy as np
import re
from dataclasses import dataclass

from kinbot import kb_path
from kinbot import constants
import math

logger = logging.getLogger('KinBot')



class Molpro:
    """
    Class to write and read molpro file and to run molpro
    """
    def __init__(self, species, par):
        self.species = species
        self.par = par

    def create_molpro_input(self, bls=0, name='', shift_vec=None, natom1=None, do_vdW=False, VTS=False, sample=False):
        """
        Create the input for molpro based on the template,
        which is either the one in the system, or provided
        by the user.
        : bls is whether it is a bls rection
        : name is an alternative name for the file
        : shift_vec is for bls to define the direction of shift in prod scan
        : natom1 is the number of atoms in the second fragment
        """
        if bls == 1:
            if shift_vec is None:
                tpl_file = self.par['barrierless_saddle_single_point_template']
            elif shift_vec is not None:
                tpl_file = self.par['barrierless_saddle_prod_single_point_template']
        elif VTS:
            if self.par['vrc_tst_scan_molpro_tpl'] == '':
                tpl_file = f'{kb_path}/tpl/molpro_vts.tpl'
            else:
                tpl_file = self.par['vrc_tst_scan_molpro_tpl']  # currently not really supported
        elif self.par['single_point_template'] == '':
            tpl_file = f'{kb_path}/tpl/molpro.tpl'
        else:
            tpl_file = self.par['single_point_template']

        with open(tpl_file) as f:
            tpl = f.read()

        fname = self.get_name(name, from_name=do_vdW)

        geom = ''
        nelectron = 0
        for i, at in enumerate(self.species.atom):
            x, y, z = self.species.geom[i]
            geom += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
            nelectron += constants.znumber[at]

        geometry_block = 'geometry={ \n'
        geometry_block += f"{self.species.natom};\n{self.species.name};\n{geom}\n"
        geometry_block += "}\n" 

        nelectron -= self.species.charge
        symm = self.molpro_symm()
        spin = self.species.mult - 1

        if not bls and not VTS:
            with open('molpro/' + fname + '.inp', 'w') as outf:
                outf.write(tpl.format(name=fname,
                                      natom=self.species.natom,
                                      geom=geom,
                                      nelectron=nelectron,
                                      symm=symm,
                                      spin=spin,
                                      charge=self.species.charge
                                      ))

        elif bls:
            closed = (nelectron - self.par['barrierless_saddle_nelectron']) / 2
            if closed.is_integer() is not True:
                logger.warning("The number of closed orbitals is not an integer,\n\
                             the CASPT2-like calculation will crash, but\n\
                             KinBot carries on for now. Revise your input,\n\
                             barrierless_saddle_nelectron is incorrect.")
            else:
                closed = int(closed)

            occ = closed + self.par['barrierless_saddle_norbital']

            if shift_vec is None:
                with open('molpro/' + fname + '.inp', 'w') as outf:
                    outf.write(tpl.format(name=fname,
                                          natom=self.species.natom,
                                          geom=geom,
                                          nelectron=nelectron,
                                          symm=symm,
                                          spin=spin,
                                          charge=self.species.charge,
                                          state=self.par['barrierless_saddle_nstate'],
                                          closed=closed,
                                          occ=occ
                                          ))
            else:
                shift_vec = shift_vec / np.linalg.norm(shift_vec) * 0.5  # step of 0.5 A
                geom0 = ''
                scancoo = ''
                scanstart = ''
                shift = ''
                for i, at in enumerate(self.species.atom):
                    x, y, z = self.species.geom[i]
                    if i < self.species.natom - natom1:
                        geom0 += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
                    else:
                        scancoo += '{} s{} s{} s{}\n'.format(at, 3 * i, 3 * i + 1, 3 * i + 2)
                        scanstart += 's{} = {:.8f}\ns{}= {:.8f}\ns{}= {:.8f}\n'.\
                                     format(3 * i, x, 3 * i + 1, y, 3 * i + 2, z)
                        shift += 's{0} = s{0} + {1:.8f}\n'.format(3 * i, shift_vec[0])
                        shift += 's{0} = s{0} + {1:.8f}\n'.format(3 * i + 1, shift_vec[1])
                        shift += 's{0} = s{0} + {1:.8f}\n'.format(3 * i + 2, shift_vec[2])
                with open('molpro/' + fname + '.inp', 'w') as f:
                    f.write(tpl.format(name=fname,
                                       natom=self.species.natom,
                                       geom=geom0,
                                       scanstart=scanstart,
                                       scancoo=scancoo,
                                       shift=shift,
                                       nelectron=nelectron,
                                       symm=symm,
                                       spin=spin,
                                       charge=self.species.charge,
                                       state=self.par['barrierless_saddle_nstate'],
                                       closed=closed,
                                       occ=occ
                                       ))
        elif VTS:
            options = 'memory,3200,M\n' \
                      'GPRINT,ORBITALS,ORBEN,CIVECTOR\n' \
                      'GTHRESH,energy=1.d-7\n' \
                      'angstrom\n' \
                      'orient,noorient\n' \
                      'nosym'

            if sample:
                basis = f"basis = {self.par['vrc_tst_sample_basis']}"
            else:
                basis = f"basis = {self.par['vrc_tst_high_basis']}"
            

            if sample:
                l3_method = self.par["vrc_tst_sample_method"]
            else:
                l3_method = self.par["vrc_tst_high_method"]

            method = ''
            match regex_in(l3_method):
                case r'.*caspt2\([0-9]+,[0-9]+\)':
                    active_electrons = int(l3_method.split("caspt2(")[1].split(",")[0])
                    active_orbitals = int(l3_method.split("caspt2(")[1].split(",")[1][:-1])
                    closed_orbitals = int(math.trunc(nelectron-active_electrons)/2)
                    occ_obitals = closed_orbitals + active_orbitals
                    method += f'{{multi,\nocc,{occ_obitals}\nclosed,{closed_orbitals}\nmaxit,50;}}\n\n'
                    method += f'put, molden, {fname}.mld\n\n'
                    method += '{rs2c, shift=0.3}\n'
                case "uwb97xd":
                    method = " {rhf;wf," + f"{nelectron},{symm},{spin},{self.species.charge}" + "}\n\n"
                    method += " omega=0.2    !range-separation parameter\n"
                    method += " srx=0.222036 !short-range exchange\n"
                    # Same as Gaussian fine grid, important for long range
                    method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                    method += " {int; ERFLERFC,mu=$omega,srfac=$srx}\n"
                    method += " uks,HYB_GGA_XC_WB97X_D\n"
                case "wb97xd":
                    method = " {rhf;wf," + f"{nelectron},{symm},{spin},{self.species.charge}" + "}\n\n"
                    method += " omega=0.2    !range-separation parameter\n"
                    method += " srx=0.222036 !short-range exchange\n"
                    method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                    method += " {int; ERFLERFC,mu=$omega,srfac=$srx}\n"
                    method += " ks,HYB_GGA_XC_WB97X_D\n"
                case "ub3lyp":
                    method = " {rhf;wf," + f"{nelectron},{symm},{spin},{self.species.charge}" + "}\n\n"
                    method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                    method += " uks,HYB_GGA_XC_B3LYP\n"
                case "b3lyp":
                    method = " {rhf;wf," + f"{nelectron},{symm},{spin},{self.species.charge}" + "}\n\n"
                    method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                    method += " ks,HYB_GGA_XC_B3LYP\n"
                case _:
                    method += " rhf\n"
                    method += " {}\n".format(l3_method)

            with open('vrctst/molpro/' + fname + '.inp', 'w') as f:
                    f.write(tpl.format(options=options,
                                       fname=fname,
                                       basis=basis,
                                       geometry_block=geometry_block,
                                       methods=method,
                                       key=self.par['vrc_tst_scan_molpro_key'].upper()))
                                       
        return 0

    def get_molpro_energy(self, key, name='', do_vdW=False, VTS=False):
        """
        Verify if there is a molpro output file and if yes, read the energy
        key is the keyword for the energy we want to read
        returns 1, energy if successful
        returns 0, -1 if the energy or the file was not there
        A non-object-oriented version is used in pes.py
        """
        fname = self.get_name(name, do_vdW)
        if fname == '10000000000000000001':  # proton
            return 1, 0.0
        if not VTS:
            status = os.path.exists('molpro/' + fname + '.out')
            if status:
                molpro_dir = "molpro/"
            else:
                status = os.path.exists('../molpro/' + fname + '.out')
                if status:
                    molpro_dir = "../molpro/"
                else:
                    molpro_dir = f"{fname.split('_')[0]}/molpro/"
                    status = os.path.exists(molpro_dir + fname + '.out')
        else:
            status = os.path.exists('vrctst/molpro/' + fname + '.out')
            if status:
                molpro_dir = 'vrctst/molpro/'
        if status:
            with open(f"{molpro_dir}{fname}.out") as f:
                lines = f.readlines()
            for line in reversed(lines):
                if ('SETTING ' + key) in line:
                    return 1, float(line.split()[3])
        return 0, -1

    def create_molpro_submit(self, name='', do_vdW=False, VTS=False):
        """
        write a pbs or slurm file for the molpro input file
        """
        fname = self.get_name(name, from_name=do_vdW)

        # open the template head and template
        if self.par['q_temp_l3'] == '':
            molpro_head = f'{kb_path}/tpl/{self.par["queuing"]}.tpl'
        else:
            molpro_head = self.par['q_temp_l3']
        with open(molpro_head) as f:
            tpl_head = f.read()
        molpro_tpl = f'{kb_path}/tpl/{self.par["queuing"]}_molpro.tpl'
        with open(molpro_tpl) as f:
            tpl = f.read()
        # substitution
        if not VTS:
            file_string = f'molpro/{fname}.{self.par["queuing"]}'
        else:
            file_string = f'vrctst/molpro/{fname}.{self.par["queuing"]}'
        with open(file_string, 'w') as f:
            if self.par['queuing'] == 'pbs':
                f.write((tpl_head + tpl).format(
                        name=fname,
                        ppn=self.par['single_point_ppn'],
                        queue_name=self.par['queue_name'],
                        errdir='molpro',
                        command=self.par['single_point_command']))
            elif self.par['queuing'] == 'slurm':
                f.write((tpl_head + tpl).format(
                        name=fname,
                        ppn=self.par['single_point_ppn'],
                        queue_name=self.par['queue_name'],
                        errdir='.',
                        command=self.par['single_point_command'],
                        slurm_feature=self.par['slurm_feature']))
        return 0

    def molpro_symm(self):
        if np.array_equal(self.species.atom, ['O']) and self.species.mult == 3:
            return 4
        if np.array_equal(self.species.atom, ['S']) and self.species.mult == 3:
            return 4
        if np.array_equal(self.species.atom, ['O', 'O']) \
                and self.species.mult == 3:
            return 4
        if np.array_equal(sorted(self.species.atom), sorted(['O', 'H'])) \
                and self.species.mult == 2:
            return 2
        return 1

    def get_name(self, name, from_name=False):
        if name != '':
            fname = name
        elif self.species.wellorts or from_name:
            fname = self.species.name
        else:
            fname = str(self.species.chemid)
        return fname

@dataclass
class regex_in:
    string: str

    def __eq__(self, other: str | re.Pattern):
        if isinstance(other, str):
            other = re.compile(other)
        assert isinstance(other, re.Pattern)
        # TODO extend for search and match variants
        return other.fullmatch(self.string) is not None
