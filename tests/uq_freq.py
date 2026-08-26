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
This class tests the frequency scaling used in the UQ analysis.

The UQ factor is drawn log-uniformly, so a factor and its reciprocal are
equally likely and have to produce reciprocal changes to every frequency.
Otherwise the perturbed frequencies carry a systematic bias, which is
largest for the low-frequency modes the scaling exists to emphasize.
"""
import unittest

from kinbot.mess import MESS

FREQS = [10., 25., 50., 100., 500., 1000., 3000.]
FACTORS = [1.05, 1.1, 1.2, 1.5, 2.0]


class StubMESS:
    """
    make_freq and scale_freq only need the parameters, not the templates
    and species that MESS.__init__ reads.
    """

    def __init__(self, ref=100., max_exp=4.):
        self.par = {'freq_uq_ref': ref, 'freq_uq_max_exp': max_exp}

    scale_freq = MESS.scale_freq
    make_freq = MESS.make_freq


class TestUQFreq(unittest.TestCase):

    def setUp(self):
        self.mess = StubMESS()

    def testReciprocal(self):
        """A factor and its reciprocal have to cancel exactly."""
        for fr in FREQS:
            for factor in FACTORS:
                up = self.mess.scale_freq(fr, factor) / fr
                down = self.mess.scale_freq(fr, 1. / factor) / fr
                warn = 'Factor {} and its reciprocal do not cancel at '.format(factor)
                warn += '{} cm-1: product is {}'.format(fr, up * down)
                self.assertAlmostEqual(up * down, 1., places=12, msg=warn)

    def testReference(self):
        """At the reference frequency the factor is applied as is."""
        ref = self.mess.par['freq_uq_ref']
        for factor in FACTORS:
            for f in (factor, 1. / factor):
                warn = 'Factor {} is not applied as is at {} cm-1'.format(f, ref)
                self.assertAlmostEqual(self.mess.scale_freq(ref, f) / ref, f,
                                       places=12, msg=warn)

    def testMonotonic(self):
        """
        The perturbation is damped as the frequency grows. Frequencies in the
        capped region share an exponent, so ties are allowed.
        """
        tol = 1e-12
        for factor in FACTORS:
            up = [self.mess.scale_freq(fr, factor) / fr for fr in FREQS]
            down = [self.mess.scale_freq(fr, 1. / factor) / fr for fr in FREQS]
            for i in range(len(FREQS) - 1):
                warn = 'Perturbation with factor {} is not '.format(factor)
                warn += 'decreasing between {} and {} cm-1'.format(FREQS[i], FREQS[i + 1])
                self.assertGreaterEqual(up[i], up[i + 1] - tol, warn)
                self.assertLessEqual(down[i], down[i + 1] + tol, warn)

    def testPositive(self):
        """Frequencies stay positive for any factor, however small."""
        for fr in FREQS:
            for factor in [0.1, 0.5, 1. / 1.2, 1., 1.2, 10.]:
                scaled = self.mess.scale_freq(fr, factor)
                warn = 'Frequency {} cm-1 is not positive after '.format(fr)
                warn += 'scaling with {}: {}'.format(factor, scaled)
                self.assertGreater(scaled, 0., warn)

    def testCap(self):
        """The amplification is frozen below ref / max_exp."""
        ref = self.mess.par['freq_uq_ref']
        max_exp = self.mess.par['freq_uq_max_exp']
        cutoff = ref / max_exp
        for factor in FACTORS:
            expected = factor ** max_exp
            for fr in (cutoff, cutoff / 2., cutoff / 10.):
                warn = 'Amplification is not capped at {} cm-1 '.format(fr)
                warn += 'for factor {}'.format(factor)
                self.assertAlmostEqual(self.mess.scale_freq(fr, factor) / fr,
                                       expected, places=12, msg=warn)

    def testLinearLimit(self):
        """For factors near 1 the old constant shift is recovered."""
        ref = self.mess.par['freq_uq_ref']
        factor = 1.001
        for fr in [fr for fr in FREQS if fr >= ref]:
            old = fr + ref * (factor - 1.)
            new = self.mess.scale_freq(fr, factor)
            warn = 'New and old scaling differ by more than 0.01% at '
            warn += '{} cm-1: {} vs {}'.format(fr, old, new)
            self.assertLess(abs(new - old) / old, 1e-4, warn)

    def testMakeFreq(self):
        """The formatted output keeps one frequency per input, three per line."""
        freqs = [100., 200., 300., 400.]
        out = self.mess.make_freq(freqs, 1.2, 0)
        self.assertEqual(len(out.split()), len(freqs))
        self.assertEqual(out.split()[0], '120.0')
        # for a saddle point the imaginary mode is dropped
        out = self.mess.make_freq(freqs, 1.2, 1)
        self.assertEqual(len(out.split()), len(freqs) - 1)


if __name__ == "__main__":
    unittest.main()
