import copy
import numpy as np

from kinbot import constants


def calc_angle(a, b, c):
    """ Calculate the A - B - C angle in radians"""

    v1 = (b-a) / np.linalg.norm(b-a)
    v2 = (b-c) / np.linalg.norm(b-c)
    return np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))


def is_linear(geom, bond):
    """
    Check if all the atoms are linear and if so, add a dummy atom
    """
    dummy = []  # positions of the dummy atoms
    for i, pos in enumerate(geom):
        # get the number of neighbors:
        neigh = len(np.where(np.asarray(bond[i]) > 0)[0])
        if neigh == 2:
            j = np.where(np.asarray(bond[i]) > 0)[0][0]
            k = np.where(np.asarray(bond[i]) > 0)[0][1]

            alpha = calc_angle(geom[j], pos, geom[k])
            if alpha > 175. * np.pi / 180.:
                # consider an angle as being linear
                # if it is larger than 175 degrees
                v1 = geom[j] - pos
                if v1[0] != 0 or v1[1] != 0:
                    p = [v1[1], -v1[0], 0.]
                elif v1[0] != 0 or v1[2] != 0:
                    p = [v1[2], 0., -v1[0]]
                else:
                    p = [1., 0., 0.]
                dummy.append(pos + p / np.linalg.norm(p))
    return dummy


def calc_out_of_plane_angle(a, b, c, d):
    """
    Calculate the out of plane angle of the A-D vector
    to the A-B-C plane

    Returns the value in radians and a boolean
    telling if b-a-c are near-collinear
    """
    collinear_cutoff = 175./180.
    collinear = 0
    if abs(calc_angle(b, a, c)) > np.pi * collinear_cutoff:
        collinear = 1
    rab = b - a
    rac = c - a
    rad = d - a
    rab /= np.linalg.norm(rab)
    rac /= np.linalg.norm(rac)
    rad /= np.linalg.norm(rad)

    n = np.cross(rab, rac)
    n /= np.linalg.norm(n)

    sin = np.dot(n, rad)
    ang = np.arcsin(sin)

    return ang, collinear


def calc_dihedral(a, b, c, d):
    """
    Calculate the A - B - C - D dihedral angle in radians.
    For collinear or close to collinear structures return a warning.

    Returns the value in degrees
    """
    collinear_cutoff = 175./180.
    collinear = 0
    if (abs(calc_angle(a, b, c)) > np.pi * collinear_cutoff or
            abs(calc_angle(b, c, d)) > np.pi * collinear_cutoff):
        collinear = 1

    b0 = a - b  # reversed on purpose
    b1 = c - b
    b2 = d - c

    # normalize b1 so that it does not influence magnitude of vector
    b1 /= np.linalg.norm(b1)

    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x)), collinear


def new_ring_dihedrals(species, instance, step_nr,
                       total_nr_of_steps, geom=None):
    """
    Calculates the required new dihedrals to create a cyclic TS
    """
    if geom is None:
        geom = species.geom
        frac = (1.+step_nr) / (total_nr_of_steps + 0.)
    else:
        # step starts at zero, not consistent with bond lenghts
        frac = 1. / (total_nr_of_steps - step_nr + 0.)
    if len(instance) > 3:
        if len(instance) == 4:
            fin_dih = 25.
        elif len(instance) == 5:
            fin_dih = 50.
        else:
            fin_dih = 1.
        dihedrals = []
        for i in range(len(instance)-3):
            dihedrals.append(calc_dihedral(geom[instance[i]],
                                           geom[instance[i+1]],
                                           geom[instance[i+2]],
                                           geom[instance[i+3]])[0])
        new_dihedrals = [dij + frac * (fin_dih - dij) for dij in dihedrals]
        return new_dihedrals


def new_bond_length(species, ati, atj, step_nr, total_nr_of_steps,
                    final_val, geom=None):
    """
    Calculates the required new bond lengths to create a TS
    """
    if geom is None:
        geom = species.geom
        frac = (0. + step_nr) / (total_nr_of_steps + 0.)
    else:
        # step starts at 1, not consistent with dihedrals
        frac = 1. / (total_nr_of_steps - step_nr + 1.)

    current_val = np.linalg.norm(geom[ati] - geom[atj])
    new_val = current_val + frac * (final_val - current_val)

    return new_val


def init_ring_dihedral(species, instance, geom=None):
    """
    Calculates the required modifications to a
    structures dihedral to create a cyclic TS
    """
    if geom is None:
        geom = species.geom
    if len(instance) > 3:
        if len(instance) < 6:
            final_dihedral = 15.
        else:
            final_dihedral = 1.
        dihedrals = []
        for i in range(len(instance)-3):
            dihedrals.append(calc_dihedral(geom[instance[i]],
                                           geom[instance[i+1]],
                                           geom[instance[i+2]],
                                           geom[instance[i+3]])[0])
        dihedral_diff = [final_dihedral - dij for dij in dihedrals]

        return dihedral_diff


def translate_and_rotate(cart, atom, i, j):
    """
    translate the molecule as such that the first
    rotor atom i is the center of rotation and the ij
    vector is along the z axis

    The atom argument is deprecated
    """
    # translate the molecule:
    trans = copy.deepcopy(cart[i])
    for ci in cart:
        ci -= trans

    # rotate the molecule
    end_vec = [0., 0., 1.]
    angle = calc_angle(cart[j], cart[i], end_vec)
    if angle != 0:
        axis = np.cross(cart[j], end_vec)
        axis = axis/np.linalg.norm(axis)
        a = np.cos(angle/2)
        b, c, d = -axis*np.sin(angle/2)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        rot_matrix = ([[aa+bb-cc-dd, 2*(bc+ad),   2*(bd-ac)],
                       [2*(bc-ad),   aa+cc-bb-dd, 2*(cd+ab)],
                       [2*(bd+ac),   2*(cd-ab),   aa+dd-bb-cc]])
        for i in range(len(cart)):
            cart[i] = np.dot(rot_matrix, cart[i])

    return cart


def get_center_of_mass(geom, atom):
    list = range(len(atom))

    com = [0., 0., 0.]
    tot_mass = 0.

    for i in list:
        at = atom[i]
        tot_mass += constants.exact_mass[at]
        com += constants.exact_mass[at] * geom[i]

    return com / tot_mass


def rotate_atom(v, n, th):
    """
    Rotate vector v around unit vector n by angle th in 3D.
    """
    w = np.zeros(3)
    w[0] = (v[0] * ms(n[0], th) +
            v[1] * mm(n[0], n[1], n[2], th) +
            v[2] * mp(n[0], n[2], n[1], th))

    w[1] = (v[0] * mp(n[1], n[0], n[2], th) +
            v[1] * ms(n[1], th) +
            v[2] * mm(n[1], n[2], n[0], th))

    w[2] = (v[0] * mm(n[2], n[0], n[1], th) +
            v[1] * mp(n[2], n[1], n[0], th) +
            v[2] * ms(n[2], th))

    v = w

    return v


def ms(x, a):
    """
    Diagonal element of the rotation matrix.
    x: selected coordinate of the unit vector around which
    rotation is done.
    a: angle
    """
    return x * x * (1. - np.cos(a)) + np.cos(a)


def mm(x, y, z, a):
    """
    Off-diagonal element of the rotation matrix with minus sign.
    x, y, x: coordinates of the unit vector around which
    rotation is done. Order matters!
    a: angle
    """
    return x * y * (1. - np.cos(a)) - z * np.sin(a)


def mp(x, y, z, a):
    """
    Off-diagonal element of the rotation matrix with plus sign.
    x, y, x: coordinates of te unit vector around which
    rotation is done. Order matters!
    a: angle
    """
    return x * y * (1. - np.cos(a)) + z * np.sin(a)


def get_moments_of_inertia(geom, atom):
    # translate the molecule to the center of mass
    geom = geom - get_center_of_mass(geom, atom)
    # initialize the inertia tensor
    I = np.zeros((3, 3))
    # add a contribution per atom
    for i, at in enumerate(atom):
        m = constants.exact_mass[at]
        x = geom[i][0]
        y = geom[i][1]
        z = geom[i][2]

        I[0][0] += m * (y**2 + z**2)
        I[1][1] += m * (x**2 + z**2)
        I[2][2] += m * (x**2 + y**2)

        I[0][1] += -m * x * y
        I[1][0] += -m * x * y

        I[1][2] += -m * y * z
        I[2][1] += -m * y * z

        I[0][2] += -m * x * z
        I[2][0] += -m * x * z

    eigvals, eigvecs = np.linalg.eigh(I)

    return eigvals, eigvecs.transpose()


def equal_geom(orig_spec, new_spec, cutoff):
    """
    Test if two geometries are the same based on:
    - bond mx has to be the same
    - bond lenths have to be within cutoff as a percentage change
    Only works for structures with unchanged atom order, e.g.,
    L2 vs L1 or conformers vs base.
    """
    
    max_bond_new = new_spec.bonds[0]
    for b in range(len(new_spec.bonds) - 1):
        max_bond_new = np.maximum(max_bond_new, new_spec.bonds[b + 1])

    max_bond_orig = orig_spec.bonds[0]
    for b in range(len(orig_spec.bonds) - 1):
        max_bond_orig = np.maximum(max_bond_orig, orig_spec.bonds[b + 1])
 
    if max_bond_orig.all() != max_bond_new.all():
        return 0

    for i in range(len(orig_spec.bond[0])-1):
        for j in range(i + 1, len(orig_spec.bond[0])):
            if orig_spec.bond[i][j] > 0:
                orig_dist = np.linalg.norm(orig_spec.geom[i] - orig_spec.geom[j])
                new_dist = np.linalg.norm(new_spec.geom[i] - new_spec.geom[j])
                if np.abs(new_dist - orig_dist) / orig_dist > cutoff:
                    return 0
    return 1


def matrix_corr(p, q):
    """
    Calculates the correlation of two sets of points, p and q,
    where each point in these sets are 3D. The correlation is calculated
    after p is rotated (r) and translated (t) such that it has maximum overlap
    with q in this sense:
    ||(r pi + t) - qi||^2
    """

    # centorids of both point sets
    # pcent = np.zeros(3)
    # for pi in p:
    #     pcent += pi
    #     pcent[0] += pi[0]
    #     pcent[1] += pi[1]
    #     pcent[2] += pi[2]

    # qcent = np.zeros(3)
    # for qi in q:
    #     qcent[0] += qi[0]
    #     qcent[1] += qi[1]
    #     qcent[2] += qi[2]

    # l = len(p)
    # pcent = pcent / l
    # qcent = qcent / l

    pcent = p.mean(axis=0)
    qcent = q.mean(axis=0)

    # shift to centroids
    x = p - pcent
    y = q - qcent

    # 3 by 3 covariance matrix
    s = np.matmul(np.transpose(x), y)
    # s = np.matmul(x.T, y) #  only Python 3.6+ compatible

    # SVD of the covariance mx
    u, _, v = np.linalg.svd(s)

    # to determine the sign (reflection, if needed)
    diag = np.identity(3)
    diag[2][2] = np.linalg.det(np.matmul(np.transpose(v), np.transpose(u)))  # +1 or -1
    # diag[2][2] = np.linalg.det(v.T @ u.T) #  only Python 3.6+ compatible

    # rotation
    r = np.matmul(np.matmul(np.transpose(v), diag), np.transpose(u))
    # translation
    t = qcent - np.matmul(r, pcent)

    # rotated and translated p
    pnew = np.zeros([len(p), 3])  # initialize
    for i, pi in enumerate(p):
        pnew[i] = np.matmul(r, pi) + t

    # return(np.corrcoef(np.matrix.flatten(pnew),np.matrix.flatten(q))[0][1])
    return np.corrcoef(pnew.ravel(), q.ravel())[0][1]
