import os

from kinbot import kb_path
from kinbot import constants


class Orca:
    """
    Class to write and read orca file and to run orca
    """
    def __init__(self, species, par):
        self.species = species
        self.par = par

    def create_orca_input(self, name=''):
        """
        Create the input for orca based on the template,
        which is provided by the use.
        : bls is whether it is a bls rection
        """
        tpl_file = self.par['single_point_template']

        with open(tpl_file) as f:
            tpl = f.read()

        fname = self.get_name(name)

        geom = ''
        nelectron = 0
        for i, at in enumerate(self.species.atom):
            x, y, z = self.species.geom[i]
            geom += '{} {:.8f} {:.8f} {:.8f}\n'.format(at, x, y, z)
            nelectron += constants.znumber[at]

        with open(f'orca/{fname}.inp', 'w') as outf:
            outf.write(tpl.format(geom=geom,
                                  mult=self.species.mult,
                                  charge=self.species.charge,
                                  ppn=self.par['single_point_ppn'],
                                  ))

        return 0

    def get_orca_energy(self, key, name=''):
        """
        Verify if there is a orca output file and if yes, read the energy
        key is the keyword for the energy we want to read
        returns 1, energy if successful
        returns 0, -1 if the energy or the file was not there
        A non-object-oriented version is used in pes.py
        """
        fname = self.get_name(name)
        status = os.path.exists(f'orca/{fname}_property.txt')
        if status:
            with open(f'orca/{fname}_property.txt') as f:
                lines = f.readlines()
            for index, line in enumerate(reversed(lines)):
                # in log file
                # E(CCSD(T))                                 ...    -75.637732066
                # in property file (used)
                # Total MDCI Energy:                                                -75.6377320661
                if (key) in line:
                    return 1, float(line.split()[-1])
        return 0, -1

    def create_orca_submit(self, name=''):
        """
        write a submission file for the orca input file
        """
        fname = self.get_name(name)

        # open the template head and template
        if self.par['queue_template'] == '':
            orca_head = f'{kb_path}/tpl/{self.par["queuing"]}.tpl'
        else:
            orca_head = self.par['queue_template'] 
        with open(orca_head) as f:
            tpl_head = f.read()
        orca_tpl = f'{kb_path}/tpl/{self.par["queuing"]}_orca.tpl'
        with open(orca_tpl) as f:
            tpl = f.read()
        # substitution
        with open(f"orca/{fname}.{self.par['queuing']}", 'w') as f:
            if self.par['queuing'] == 'pbs':
                f.write((tpl_head + tpl).format(
                        name=fname,
                        ppn=self.par['single_point_ppn'],
                        queue_name=self.par['queue_name'],
                        errdir='orca',
                        command=self.par['single_point_command']))
            elif self.par['queuing'] == 'slurm':
                f.write((tpl_head + tpl).format(
                        name=fname,
                        ppn=self.par['single_point_ppn'],
                        queue_name=self.par['queue_name'],
                        errdir='orca',
                        command=self.par['single_point_command'],
                        slurm_feature=self.par['slurm_feature']))

        return 0

    def get_name(self, name):
        if name != '':
            fname = name
        elif self.species.wellorts:
            fname = self.species.name
        else:
            fname = str(self.species.chemid)
        return fname
