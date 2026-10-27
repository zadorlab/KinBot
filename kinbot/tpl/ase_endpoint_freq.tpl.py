"""One fixed-geometry endpoint Hessian; no optimization or model fallback."""
import numpy as np
from ase import Atoms
from ase.db import connect
from kinbot.frequencies import frequency_validation_error
from kinbot.utils import iowait

label = '{label}'
mol = Atoms(symbols={atom}, positions={geom})
kwargs = {kwargs}
if '{code}' == 'gauss':
    from kinbot.ase_modules.calculators.gaussian import Gaussian
    from kinbot import reader_gauss as reader
    Gaussian.command = '{qc_command} < PREFIX.com > PREFIX.log'
    mol.calc = Gaussian(**kwargs)
    output, marker = label+'.log', label+'.log'
else:
    from kinbot.ase_modules.calculators.qchem import QChem
    from kinbot import reader_qchem as reader
    calc = QChem(**kwargs)
    calc.command = '{qc_command} -nt {ppn} PREFIX.inp PREFIX.out'
    mol.calc = calc
    output, marker = label+'_freq.out', label+'.out'
try:
    energy = mol.get_potential_energy()
    iowait(output, '{code}')
    frequencies = reader.read_freq(output, {atom})
    zpe = reader.read_zpe(output)
    reason = frequency_validation_error(frequencies, 0, {imagfreq_threshold})
    if reason or not np.isfinite(energy) or not np.isfinite(zpe):
        raise ValueError(reason or 'Nonfinite energy or ZPE.')
    data = {{'status': 'normal', 'energy': energy, 'zpe': zpe,
             'frequencies': np.asarray(frequencies)}}
except Exception as error:
    data = {{'status': 'error', 'error': str(error)}}
connect('{working_dir}/kinbot.db').write(mol, name=label, data=data)
with open(marker, 'a') as handle:
    handle.write('done\n')
