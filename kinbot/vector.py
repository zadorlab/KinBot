###################################################
##                                               ##
## This file is part of the KinBot code v2.0     ##
##                                               ##
## The contents are covered by the terms of the  ##
## BSD 3-clause license included in the LICENSE  ##
## file, found at the root.                      ##
##                                               ##
## Copyright 2018 National Technology &          ##
## Engineering Solutions of Sandia, LLC (NTESS). ##
## Under the terms of Contract DE-NA0003525 with ##
## NTESS, the U.S. Government retains certain    ##
## rights to this software.                      ##
##                                               ##
## Authors:                                      ##
##   Judit Zador                                 ##
##   Ruben Van de Vijver                         ##
##                                               ##
###################################################
import numpy as np



def angle(a, b, c):
    """ Calculate the A - B - C angle in radians"""

    v1 = (b-a) / np.linalg.norm(b-a)
    v2 = (b-c) / np.linalg.norm(b-c) 
    
    return np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))



def dihedral(a, b, c, d):
    """ 
    Calculate the A - B - C - D dihedral angle in radians.
    For collinear or close to collinear structures return a warning.
    """


    collinear = 0

    if abs(abs(angle(a, b, c)) - 180.) < 5. or abs(abs(angle(b, c, d)) - 180.) < 5.:
        collinear = 1        

    b0 = a - b # reversed on purpose
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



def rotate_atom(v, n, th):
    """ 
    Rotate vector v around unit vector n by angle th in 3D.
    """
    
    w = np.zeros(3)
    w[0] = v[0] * ms(n[0], th ) + v[1] * mm(n[0], n[1], n[2], th) + v[2] * mp(n[0], n[2], n[1], th)
    w[1] = v[0] * mp(n[1], n[0], n[2], th) + v[1] * ms(n[1], th) + v[2] * mm(n[1], n[2], n[0], th)
    w[2] = v[0] * mm(n[2], n[0], n[1], th) + v[1] * mp(n[2], n[1], n[0], th) + v[2] * ms(n[2], th)
    
    v = w
    
    return v



def ms(x, a):
    """
    Diagonal element of the rotation matrix. 
    x: selected coordinate of the unit vector around which rotation is done.
    a: angle
    """
    
    
    return x * x * (1. - np.cos(a)) + np.cos(a)



def mm(x, y, z, a):
    """
    Off-diagonal element of the rotation matrix with minus sign.
    x, y, x: coordinates of the unit vector around which rotation is done. Order matters!
    a: angle
    """
    
    
    return x * y * (1. - np.cos(a)) - z * np.sin(a)



def mp(x, y, z, a):
    """
    Off-diagonal element of the rotation matrix with plus sign.
    x, y, x: coordinates of te unit vector around which rotation is done. Order matters!
    a: angle
    """
    
    
    return x * y * (1. - np.cos(a)) + z * np.sin(a)



def main():
    """
    Simple vector algebra.
    """



if __name__ == "__main__":
    main()
