import numpy as np

def read_imag_mode(job, natom):
    """
    Read the imaginary normal mode displacements from a log file.
    Only for saddle points! It will read the firs normal mode
    for a well, but that's not very useful.
    """

    nmode = np.zeros([natom, 3])
    joblog = '{}_vib/vib.xyz'.format(job)
    with open(joblog) as f:
        lines = f.read().split('\n')

    l = 1
    line = lines[l]
    for n in range(natom):
        mm = lines[l + n + 1].split()
        nmode[n][0] = float(mm[-3])
        nmode[n][1] = float(mm[-2])
        nmode[n][2] = float(mm[-1])

    return nmode

