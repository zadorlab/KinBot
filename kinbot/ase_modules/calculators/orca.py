import re

import ase.io.orca as io
from ase.calculators.genericfileio import (
    BaseProfile,
    CalculatorTemplate,
    GenericFileIOCalculator,
)


def get_version_from_orca_header(orca_header):
    match = re.search(r'Program Version (\S+)', orca_header, re.M)
    return match.group(1)


class OrcaProfile(BaseProfile):
    def version(self):
        # XXX Allow MPI in argv; the version call should not be parallel.
        from ase.calculators.genericfileio import read_stdout
        stdout = read_stdout([self.command, "does_not_exist"])
        return get_version_from_orca_header(stdout)

    def get_calculator_command(self, inputfile):
        return [inputfile]


class OrcaTemplate(CalculatorTemplate):
    _label = 'orca'

    def __init__(self):
        super().__init__('orca',
                         implemented_properties=['energy', 'free_energy',
                                                 'forces', 'dipole'])

        self.inputname = f'{self._label}.inp'
        self.outputname = f'{self._label}.out'
        self.errorname = f'{self._label}.err'

    def execute(self, directory, profile) -> None:
        profile.run(directory, self.inputname, self.outputname,
                    errorfile=self.errorname)

    def write_input(self, profile, directory, atoms, parameters, properties):
        parameters = dict(parameters)

        kw = dict(charge=0, mult=1, orcasimpleinput='B3LYP def2-TZVP',
                  orcablocks='%pal nprocs 1 end')
        kw.update(parameters)

        io.write_orca(directory / self.inputname, atoms, kw)

    def read_results(self, directory):
        return io.read_orca_outputs(directory, directory / self.outputname)

    def load_profile(self, cfg, **kwargs):
        return OrcaProfile.from_config(cfg, self.name, **kwargs)

class ORCA(GenericFileIOCalculator):
    """Class for doing ORCA calculations.

    Example:

      calc = ORCA(charge=0, mult=1, orcasimpleinput='B3LYP def2-TZVP',
        orcablocks='%pal nprocs 16 end')
    """

    def __init__(self, *, profile=None, directory='.', **kwargs):
        """Construct ORCA-calculator object.

        Parameters
        ==========
        charge: int

        mult: int

        orcasimpleinput : str

        orcablocks: str


        Examples
        ========
        Use default values:

        >>> from ase.calculators.orca import ORCA
        >>> h = Atoms(
        ...     'H',
        ...     calculator=ORCA(
        ...         charge=0,
        ...         mult=1,
        ...         directory='water',
        ...         orcasimpleinput='B3LYP def2-TZVP',
        ...         orcablocks='%pal nprocs 16 end'))

        """

        super().__init__(template=OrcaTemplate(),
                         profile=profile, directory=directory,
                         parameters=kwargs)
