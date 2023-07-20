import numpy as np
import math

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    u = vector / np.linalg.norm(vector)
    return u

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2':
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def plane_from_points(v0, v1, v2):
    u = np.subtract(v1, v0)
    v = np.subtract(v2, v0)

    normal = np.cross(u, v)
    d= np.dot(normal, v0)

    return normal, d

def dist_point_to_plane(point, plane):
    abc, d = plane
     
    dist = abs((np.dot(abc, point) - d))
    e = np.sqrt(np.sum(np.square(abc)))
    return dist/e

from kinbot.pp_tables import *

v1 = np.array([1, 0, 3])
v2 = np.array([1, 1, 9])
v3 = np.array([1, 2, 9])
ra_pos = np.array([1.05, -1, 9])

n_pp = 2
plane = plane_from_points(v1, v2, v3)
ra_to_plane = dist_point_to_plane(ra_pos, plane)
if abs(ra_to_plane) >= .1:
    #If carbon atom out of plane, only place a single pp on other side of the plane
    n_pp = 1
pp_list = []
plane_direction = unit_vector(np.dot(plane[0], ra_pos))
length = 0.5
for i in range(n_pp):
    pp_orient = np.array(unit_vector(plane[0])*plane_direction*np.power(-1,i), dtype=float)
    pp_vect = length * unit_vector(pp_orient)
    pp_coord = np.add(ra_pos, pp_vect)
    pp_list.append(pp_coord)
