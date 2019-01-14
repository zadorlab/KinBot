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
"""
This class tests the BFGS local optimizer
"""
import unittest
import numpy as np

import kinbot.bfgs as bfgs


class styblinski_tang():
    """
    n-dimensional, continuous, multimodal, non-convex

    global minimum:
    xi = -2.90353403 for xi in x
    local minima where xi can be 2.74680277
    """
    def __init__(self):
        pass

    def eval(self, x):
        sum = 0.
        for xi in x:
            sum += xi**4 - 16 * xi**2 + 5 * xi
        return 0.5*sum

    def gradient(self, x):
        grad = np.zeros(len(x))
        for i, xi in enumerate(x):
            grad[i] += 2. * xi**3 - 16. * xi + 2.5
        return grad


class TestBFGS(unittest.TestCase):
    def setUp(self):
        pass

    def test(self):
        f = styblinski_tang()
        # use a 3 dimensional function initiated
        # at random numbers (with seed)
        np.random.seed(1)
        x = np.random.uniform(low=-5, high=5, size=(3, ))
        opt = bfgs.BFGS()
        xk, x_i, g_i = opt.optimize(f, x)
        x_expected = [-2.90353403, 2.74680277, 2.74680277]
        for i, xki in enumerate(xk):
            self.assertAlmostEqual(x_expected[i], xki, places=6)


if __name__ == '__main__':
    unittest.main()
