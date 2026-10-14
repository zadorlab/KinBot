from kinbot.species_routing import routing_name
import os
import logging
import numpy as np
import re
from dataclasses import dataclass

from kinbot import kb_path
from kinbot.theory import fortran_threshold
from kinbot import constants
import math

logger = logging.getLogger('KinBot')


_VRC_MULTIREFERENCE = re.compile(
    r'(caspt2|mrci)(?:-f12([ab])?)?(\+q)?\((\d+),(\d+)\)')


def vrc_multireference_energy_contract(method):
    """Return the exact Molpro energy selected for a VRC MR method.

    Molpro stores the F12A and F12B results of an unqualified F12 command in
    indices 1 and 2, respectively.  The literature shorthand ``CASPT2-F12``
    and ``MRCI-F12+Q`` is interpreted as F12B here.  Users can still request
    F12A explicitly.  MRCI ``+Q`` always means the relaxed-reference Davidson
    correction (``OPTION,CLUSTER=0``), including the default CABS-singles
    correction in an F12 calculation.
    """
    source = str(method).strip().casefold()
    match = _VRC_MULTIREFERENCE.fullmatch(source)
    if match is None:
        return None
    family, variant, correction, electrons_text, orbitals_text = \
        match.groups()
    if family == 'caspt2' and correction is not None:
        return None
    if family == 'mrci' and correction != '+q':
        return None
    f12 = '-f12' in source
    variant = variant.upper() if f12 and variant else ('B' if f12 else None)
    index = 2 if variant == 'B' else 1
    root = 'CASPT2' if family == 'caspt2' else 'MRCI'
    label = root + (f'-F12{variant}' if f12 else '')
    if family == 'mrci':
        label += '+Q'
    return {
        'method': label,
        'active_electrons': int(electrons_text),
        'active_orbitals': int(orbitals_text),
        'f12_variant': variant,
        'energy_variable': (f'ENERGD({index})' if family == 'mrci'
                            else f'ENERGY({index})'),
        'cabs_singles_included': (
            True if f12 and family == 'mrci' else (False if not f12 else None)),
        'davidson_correction': (
            'relaxed_reference' if family == 'mrci' else None),
        'molpro_cluster_option': 0 if family == 'mrci' else None,
    }


def _vrc_multireference_method(method, nelectron, spin=0, *, shift=None,
                               ipea=None):
    """Render a supported VRC active-space method and its energy variable.

    KinBot uses ``CASPT2(ne,no)`` for the inexpensive ROTD sampling surface
    and can use ``MRCI+Q(ne,no)`` for the trusted one-dimensional correction.
    Molpro obtains the Davidson-corrected MRCI energy from ``ENERGD``; ``MRCI``
    followed by ``energy`` would silently select the variational MRCI value.
    """
    contract = vrc_multireference_energy_contract(method)
    if contract is None:
        return None
    family = 'caspt2' if contract['method'].startswith('CASPT2') else 'mrci'
    active_electrons = contract['active_electrons']
    active_orbitals = contract['active_orbitals']
    if active_electrons < 1 or active_orbitals < 1:
        raise ValueError('VRC active electrons and orbitals must be positive.')
    inactive_electrons = int(nelectron) - active_electrons
    if inactive_electrons < 0 or inactive_electrons % 2:
        raise ValueError(
            f'{method} is incompatible with {nelectron} electrons: the '
            'inactive electron count must be nonnegative and even.')
    closed_orbitals = inactive_electrons // 2
    occupied_orbitals = closed_orbitals + active_orbitals
    multi = ('{multi,\n'
             f'occ,{occupied_orbitals}\n'
             f'closed,{closed_orbitals}\n'
             f'wf,{int(nelectron)},1,{int(spin)}\n'
             'maxit,50;}\n\n')
    f12 = contract['f12_variant'] is not None
    selected = contract['energy_variable'].casefold()
    if family == 'caspt2':
        program = ('rs2-f12' if f12
                   else ('rs2c' if int(spin) == 0 else 'rs2'))
        options = []
        if shift is None:
            shift = 0.3
        if shift is not None:
            options.append(f'shift={float(shift):g}')
        if ipea is not None:
            options.append(f'ipea={float(ipea):g}')
        command = ','.join([program, *options])
        audit = (f'! KinBot energy contract: {contract["method"]} = '
                 f'{contract["energy_variable"]}\n'
                 f'KBVRCSEL={selected}\n')
        return multi + f'{{{command}}}\n' + audit, 'KBVRCSEL'
    program = ('mrci-f12,scalf12=1,cabs_singles=1' if f12 else 'mrci')
    variant_index = 2 if contract['f12_variant'] == 'B' else 1
    audit = (f'! KinBot energy contract: {contract["method"]} = '
             f'{contract["energy_variable"]}; Davidson relaxed reference\n'
             f'KBVRCVAR=energy({variant_index})\n'
             f'KBVRCSEL={selected}\n')
    return (multi + f'{{{program};option,cluster=0;}}\n' + audit,
            'KBVRCSEL')



class Molpro:
    """
    Class to write and read molpro file and to run molpro
    """
    def __init__(self, species, par):
        self.species = species
        self.par = par

    def create_molpro_input(self, bls=0, name='', shift_vec=None, natom1=None,
                            do_vdW=False, VTS=False, sample=False,
                            vrc_level=None):
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
        if vrc_level is not None:
            multiplicity = vrc_level.get('multiplicity', self.species.mult)
            if multiplicity is None:
                multiplicity = self.species.mult
            if (isinstance(multiplicity, bool)
                    or not isinstance(multiplicity, int)
                    or multiplicity < 1):
                raise ValueError('VRC multiplicity must be a positive integer.')
            spin = multiplicity - 1
            if spin > nelectron or (nelectron - spin) % 2:
                raise ValueError('VRC multiplicity is incompatible with the '
                                 'molecular electron count.')

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
            scf_threshold = fortran_threshold(
                self.par.get('vrc_tst_scf_convergence', 1.e-10),
                'vrc_tst_scf_convergence')
            correlation_threshold = fortran_threshold(
                self.par.get('vrc_tst_correlation_convergence', 1.e-8),
                'vrc_tst_correlation_convergence')
            options = 'memory,3200,M\n' \
                      'GPRINT,ORBITALS,ORBEN,CIVECTOR\n' \
                      f'GTHRESH,energy={scf_threshold},' \
                      f'coeff={correlation_threshold},' \
                      f'civec={correlation_threshold}\n' \
                      'angstrom\n' \
                      'orient,noorient\n' \
                      'nosym'

            if vrc_level is None:
                vrc_level = {
                    'method': (self.par['vrc_tst_sample_method'] if sample
                               else self.par['vrc_tst_high_method']),
                    'basis': (self.par['vrc_tst_sample_basis'] if sample
                              else self.par['vrc_tst_high_basis']),
                }
            basis = f"basis = {vrc_level['basis']}"
            l3_method = vrc_level['method']

            method = ''
            energy_expression = 'energy(1)'
            multireference = _vrc_multireference_method(
                l3_method, nelectron, spin=spin,
                shift=vrc_level.get('shift'), ipea=vrc_level.get('ipea'))
            if multireference is not None:
                method, energy_expression = multireference
                # Keep the orbitals used by the correction calculation for
                # direct inspection of an external validation run.
                first_program = method.find('\n\n')
                method = (method[:first_program + 2]
                          + f'put,molden,{fname}.mld\n\n'
                          + method[first_program + 2:])
            else:
                match regex_in(str(l3_method).casefold()):
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
                        method += (" {rhf;wf," +
                                   f"{nelectron},{symm},{spin}" + "}\n")
                        method += "{}\n".format(l3_method)
                        energy_expression = 'energy'

            with open('vrctst/molpro/' + fname + '.inp', 'w') as f:
                    f.write(tpl.format(options=options,
                                       fname=fname,
                                       basis=basis,
                                       geometry_block=geometry_block,
                                       methods=method,
                                       key=self.par['vrc_tst_scan_molpro_key'].upper(),
                                       energy_expression=energy_expression))
                                       
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
            def setting(name):
                pattern = re.compile(
                    rf'^\s*SETTING\s+{re.escape(name)}\s*=\s*'
                    r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[DEde][-+]?\d+)?)')
                for line in reversed(lines):
                    found = pattern.search(line)
                    if found:
                        return float(found.group(1).replace('D', 'E')
                                     .replace('d', 'e'))
                return None

            energy = setting(key)
            if energy is not None:
                # Generated multireference inputs bind the selected
                # ENERGY/ENERGD array element to this audit variable first.
                # Reject a result if the public key does not reproduce it.
                input_path = f"{molpro_dir}{fname}.inp"
                if os.path.isfile(input_path):
                    with open(input_path, encoding='utf-8') as source:
                        audited = 'KBVRCSEL' in source.read()
                    if audited:
                        selected = setting('KBVRCSEL')
                        if (selected is None or not np.isclose(
                                energy, selected, rtol=0., atol=1.e-10)):
                            return 0, -1
                return 1, energy
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
        from kinbot.qc import _ordinary_slurm_walltime, _slurm_exclude_directive
        walltime_directive = (_ordinary_slurm_walltime(
            self.par, self.par['queue_name'])
            if '{walltime_directive}' in tpl_head else '')
        slurm_exclude_directive = _slurm_exclude_directive(self.par)
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
                        slurm_feature=self.par['slurm_feature'],
                        slurm_exclude_directive=slurm_exclude_directive,
                        walltime_directive=walltime_directive))
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
            fname = routing_name(self.species)
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
