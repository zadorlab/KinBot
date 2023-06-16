import time
import os


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
    maxtime = 10  # s
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
    if not os.path.exists('perm'):
        os.makedirs('perm')
    if not os.path.exists('scratch'):
        os.makedirs('scratch')
    if not os.path.exists(par['single_point_qc']):
        os.mkdir(par['single_point_qc'])
    if par['rotor_scan'] == 1:
        if not os.path.exists('hir'):
            os.mkdir('hir')
        if not os.path.exists('hir_profiles'):
            os.mkdir('hir_profiles')
        if not os.path.exists('perm/hir/'):
            os.makedirs('perm/hir/')
    if par['conformer_search'] == 1:
        if not os.path.exists('conf'):
            os.mkdir('conf')
        if not os.path.exists('perm/conf'):
            os.makedirs('perm/conf')
    if par['calc_aie'] == 1:
        if not os.path.exists('aie'):
            os.mkdir('aie')
        if not os.path.exists('perm/aie'):
            os.makedirs('perm/aie')
    if not os.path.exists('me'):
        os.mkdir('me')


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
        if 'com' in ff:
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
            except StopIteration:
                continue
            else:
                if len(atoms.positions) > 1 and np.all(atoms.positions == 0):
                    os.remove(ll)
                    logger.info(f'All coordinates of file {ll} are 0, hence '
                                 f'it is deleted.')
