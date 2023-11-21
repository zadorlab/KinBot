import os
import logging
import numpy as np

from kinbot import kb_path
from kinbot import constants
from kinbot import zmatrix


logger = logging.getLogger('KinBot')



class Molpro:
    """
    Class to write and read molpro file and to run molpro
    """
    def __init__(self, species, par):
        self.species = species
        self.par = par

    def create_molpro_input(self, bls=0, name='', shift_vec=None, natom1=None, do_vdW=False, VTS=False):
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
            tpl_file = self.par["vrc_tst_scan_parameters"]["molpro_tpl"]
        elif self.par['single_point_template'] == '':
            tpl_file = f'{kb_path}/tpl/molpro.tpl'
        else:
            tpl_file = self.par['single_point_template']

        with open(tpl_file) as f:
            tpl = f.read()

        fname = self.get_name(name, from_name=do_vdW)

        geom_type = "xyz"

        geom = ''
        nelectron = 0
        for i, at in enumerate(self.species.atom):
            x, y, z = self.species.geom[i]
            geom += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
            nelectron += constants.znumber[at]

        geometry_block = "geometry={ \n"
        if geom_type == "xyz":
            geometry_block += f"{self.species.natom};\n{self.species.name};\n{geom}\n"
        #elif geom_type == "zmat":
        #    distances, angles, dihedrals = zmatrix.make_simple_zmat_from_cart(self.species)
        #    for index in range(self.species.natom):
        #        geometry_block += f"{self.species.atom[index]}{index}"
        #        if index > 0:
        #                geometry_block += f", {self.species.atom[index-1]}{index-1}, dist_{index}"
        #                variables_block += f"{'dist_':>10s}{index} = {distances[index-1]:>12.5f};\n"
        #        if index > 1:
        #                geometry_block += f", {self.species.atom[index-2]}{index-2}, angle_{index}"
        #                variables_block += f"{'angle_':>10s}{index} = {angles[index-2]:>12.5f};\n"
        #        if index > 2:
        #                geometry_block += f", {self.species.atom[index-3]}{index-3}, dihed_{index}"
        #                variables_block += f"{'dihed_':>10s}{index} = {dihedrals[index-3]:>12.5f};\n"
        #        geometry_block += "\n"
        geometry_block += "}\n" 

        nelectron -= self.species.charge
        symm = self.molpro_symm()
        spin = self.species.mult - 1

        if bls == 0 and not VTS:
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
            options = "GPRINT,ORBITALS,ORBEN,CIVECTOR \nGTHRESH,energy=1.d-7 \nangstrom \n orient,noorient\n nosym"
            
            if "sample" in fname or "frozen" in fname:
                basis = f"basis = {self.par['vrc_tst_scan_basis']['L3'][0]}"
            else:
                basis = f"basis = {self.par['vrc_tst_scan_basis']['L3'][1]}"
            active_orbitals = 2
            active_electrons = 2
            closed_orbitals = (nelectron-active_electrons)/2
            occ_obitals = closed_orbitals + active_orbitals
            
            method = " {rhf;wf," + f"{nelectron},{symm},{spin},{self.species.charge}" + "}\n\n"

            if "sample" in fname or "frozen" in fname:
                l3_method = self.par["vrc_tst_scan_methods"]["L3"][0]
            else:
                l3_method = self.par["vrc_tst_scan_methods"]["L3"][1]
            match l3_method:
                case "caspt2":
                    method += " {multi,\n" + f" occ,{occ_obitals}\n closed,{closed_orbitals}\n" + " }\n\n"
                    if spin == 0:
                        method += " {rs2c, shift=0.3}\n"
                    else:
                        method += " {rs2, shift=0.3}\n"
                case "ccsd(t)":
                    method += " {ccsd(t)-f12}\n"
                case "uwb97xd":
                    method += " omega=0.2    !range-separation parameter\n"
                    method += " rx=0.222036 !short-range exchange\n"
                    method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                    method += " {int; ERFLERFC,mu=$omega,srfac=$srx}\n"
                    method += " uks,HYB_GGA_XC_WB97X_D\n"
                case "wb97xd":
                    method += " omega=0.2    !range-separation parameter\n"
                    method += " rx=0.222036 !short-range exchange\n"
                    method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                    method += " {int; ERFLERFC,mu=$omega,srfac=$srx}\n"
                    method += " ks,HYB_GGA_XC_WB97X_D\n"
                case "ub3lyp":
                    method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                    method += " uks,HYB_GGA_XC_B3LYP\n"
                case "b3lyp":
                    method += " {grid,wcut=1d-30,min_nr=[175,250,250,250],max_nr=[175,250,250,250],min_L=[974,974,974,974],max_L=[974,974,974,974]}\n"
                    method += " ks,HYB_GGA_XC_B3LYP\n"
                case _:
                    raise NotImplementedError
                

            with open('molpro/' + fname + '.inp', 'w') as f:
                    f.write(tpl.format(options=options,
                                       fname=fname,
                                       basis=basis,
                                       geometry_block=geometry_block,
                                       methods=method,
                                       key=self.par['vrc_tst_scan_parameters']["molpro_key"].upper()))
                                       
        return 0

    def get_molpro_energy(self, key, name='', do_vdW=False):
        """
        Verify if there is a molpro output file and if yes, read the energy
        key is the keyword for the energy we want to read
        returns 1, energy if successful
        returns 0, -1 if the energy or the file was not there
        A non-object-oriented version is used in pes.py
        """
        fname = self.get_name(name, do_vdW)
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
        if fname == '10000000000000000001':  # proton
            return 1, 0.0
        if status:
            with open(f"{molpro_dir}{fname}.out") as f:
                lines = f.readlines()
            for line in reversed(lines):
                if ('SETTING ' + key) in line:
                    return 1, float(line.split()[3])
        return 0, -1

    def create_molpro_submit(self, name='', do_vdW=False):
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
        with open(f'molpro/{fname}.{self.par["queuing"]}', 'w') as f:
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
