import time
import os
import numpy as np
from numpy.typing import NDArray
import json
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from kinbot.stationary_pt import StationaryPoint


def tail(file_path, lines=10):
    """Returns the specified last number of lines of a file.

    @param file_path: The n to retrieve the last lines from.
    @param lines: The number of lines to be retrieved.
    @return str: The last number of lines
    """
    with open(file_path, 'br') as f:
        total_lines_wanted = lines

        block_size = 1024
        f.seek(0, 2)
        block_end_byte = f.tell()
        lines_to_go = total_lines_wanted
        block_number = -1
        blocks = []
        while lines_to_go > 0 and block_end_byte > 0:
            if block_end_byte - block_size > 0:
                f.seek(block_number * block_size, 2)
                blocks.append(f.read(block_size))
            else:
                f.seek(0, 0)
                blocks.append(f.read(block_end_byte))
            lines_found = blocks[-1].count(b'\n')
            lines_to_go -= lines_found
            block_end_byte -= block_size
            block_number -= 1
        all_read_text = b''.join(reversed(blocks))
    return b'\n'.join(all_read_text.splitlines()[-total_lines_wanted:]).decode()


def iowait(logfile, qc_code):
    """Wait for I/O to be done to a specific file.
    There is a maxtime.
    """
    maxtime = 30  # s
    clock = 0

    if qc_code == 'gauss':
        endings = ['Normal termination', 'Error termination']
    elif qc_code == 'qchem':
        endings = ['Thank you very much for using Q-Chem', 'Please submit a crash report']
    else:
        endings = []
    while True:
        logtail = tail(logfile)
        if any([ee in logtail for ee in endings]):
            break
        clock += 1
        if clock > maxtime:
            break
        time.sleep(1)


def make_dirs(par):
    """Creates the directory hierarchy.

    @param par: Dictionary with the simulation parameters read from the input.
    """
    make_dir('perm')
    make_dir('scratch')
    make_dir(par['single_point_qc'])
    make_dir('me')
    if par['rotor_scan'] == 1:
        make_dir('hir')
        make_dir('hir_profiles')
        make_dir('perm/hir/')
    if par['conformer_search'] == 1:
        make_dir('conf')
        make_dir('perm/conf')
    if par['calc_aie'] == 1:
        make_dir('aie')
        make_dir('perm/aie')
    if par['vrc_tst_scan'] != {}:
        make_dir('vrctst')
        make_dir('vrctst/molpro')
        make_dir('perm/vrctst')
        make_dir('perm/vrctst/molpro')

def make_dir(name):
    '''
    Helper function for make_dirs
    '''
    if not os.path.exists(name):
        os.makedirs(name)
    return

def clean_files():
    """Removes files from jobs that ended up erroneously.
    """
    import logging
    import numpy as np
    from kinbot.ase_modules.io.formats import read
    logger = logging.getLogger('KinBot')
    # delete leftover AM1 calculations
    files = os.listdir()
    com = []
    for ff in files:
        if ff.endswith('.com'):
            com.append(ff)

    for cc in com:
        delfile = False
        with open(cc, 'r') as f:
            if 'am1' in f.read():
                delfile = True
        if delfile:
            ll = cc.split('.')[0] + '.log'
            try:
                os.remove(ll)
                logger.info(f'Stuck AM1 job {ll} is deleted.')
            except FileNotFoundError:
                pass

    # delete empty files
    try:
        conf_files = os.listdir('conf')
        conf_files = [f'conf/{ff}' for ff in conf_files]
    except FileNotFoundError:
        conf_files = []
    try:
        hir_files = os.listdir('hir')
        hir_files = [f'hir/{ff}' for ff in hir_files]
    except FileNotFoundError:
        hir_files = []
    files = files + conf_files + hir_files
    log = []
    for ff in files:
        if len(ff) > 4 and ff[-4:] == '.log':
            log.append(ff)

    for ll in log:
        if not os.path.isfile(ll):
            continue
        if os.path.getsize(ll) < 10:
            os.remove(ll)
            logger.info(f'Empty file {ll} is deleted.')
        else:
            try:
                atoms = read(ll)
            except (StopIteration, ValueError):
                continue
            else:
                if len(atoms.positions) > 1 and np.all(atoms.positions == 0):
                    os.remove(ll)
                    logger.info(f'All coordinates of file {ll} are 0, hence '
                                f'it is deleted.')


def create_matplotlib_graph(x=[0., 1.],
                            data=[[1., 1.]],
                            name="mtpltlb",
                            x_label="x",
                            y_label="y",
                            data_legends=["y0"],
                            comments=[""]):
    """Function that creates the input for a 2D matplotlib plot."""

    if not isinstance(x, list) and not isinstance(data, list):
        return

    scale = 5
    width = 3.236*scale
    height = 2*scale

    x_spln = np.arange(min(x), max(x), 0.1)

    spln = []
    y_spln = []
    min_y = np.inf
    max_y = -np.inf
    for index, y in enumerate(data):
        if min(y) < min_y:
            min_y = min(y)
        if max(y) > max_y:
            max_y = max(y)
        spln.append(make_interp_spline(x, y))
        y_spln.append(spln[index](x_spln))

    fig, ax = plt.subplots()

    markers = ['x', 'o', '+', 'h', '*']
    for index, legend in enumerate(data_legends):
        ax.plot(x, data[index],
                label=legend,
                marker=markers[index],
                linewidth=3.0,
                markersize=15)
        # ax.plot(x_spln, y_spln[index], label=f'spln_{legend}')

    ax.legend(loc=0)
    ax.set_xlabel(r'{}'.format(x_label))
    ax.set_ylabel(r'{}'.format(y_label))
    ax.set_xlim([min(x), max(x)])
    ax.set_ylim([min_y, max_y])

    plt.rcParams.update({'font.size': 30})
    fig.set_figheight(height)
    fig.set_figwidth(width)

    plt.savefig(f'{name}.png',
                transparent=True,
                dpi=400,
                format='png',
                )


def queue_command(qu):
    if qu == 'pbs':
        cmd = 'qsub'
        ext = 'pbs'
    elif qu == 'slurm':
        cmd = 'sbatch'
        ext = 'slurm'
    elif qu == 'local':
        cmd = ''
        ext = ''
        pass
    else:
        raise ValueError(f'Unexpected value for queueing: {par["queuing"]}')
    return cmd, ext


def reorder_coord(mol_A: StationaryPoint,
                  mol_B: StationaryPoint,
                  map_B: NDArray | None = None
                  ) -> None | NDArray:
    """Reorder the coordinates of mol_B to correspond to mol_A.

    Args:
        mol_A (StationaryPoint): Stationary point of molecule A
        mol_B (StationaryPoint): Stationary point of molecule B
        map_B (NDArray): if B is a fragment,
                         change the map to the parent accordingly
    """

    if mol_A.chemid != mol_B.chemid:
        raise AttributeError("Reordering only for identical molecules.")
    new_geom = np.array(mol_B.geom)
    idx_used = []
    new_map = np.zeros(len(mol_B.atom), dtype=int)
    for idxb, aidb in enumerate(mol_B.atomid):
        for idxa, aida in enumerate(mol_A.atomid):
            if aida == aidb and idxa not in idx_used:
                new_geom[idxa] = mol_B.geom[idxb]
                idx_used.append(idxa)
                if map_B is not None:
                    new_map[idxa] = map_B[idxb]
                break
    mol_B.atom = mol_A.atom
    mol_B.geom = new_geom
    mol_B.reset_order()
    if map_B is not None:
        return new_map
    return


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
