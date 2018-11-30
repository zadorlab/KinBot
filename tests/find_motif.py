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
This class tests the resonance structure generation algorithm

The data is a dictionary of which the keys are the smiles
and the values are the expected number of resonance isomers
"""
import unittest

from kinbot.stationary_pt import StationaryPoint
import kinbot.find_motif as find_motif


class TestFindMotif(unittest.TestCase):
    def setUp(self):
        pass

    def testNumberOfHits(self):
        """
        Test the molecular formula creator, which is
        basically a RDKit feature, but used in KinBot
        for postprocessing
        """
        # the data contains the smiles as key and as value:
        # 1. the multiplicity
        # 2. the motif which will be searched for
        # 3. the expected number of hits of the motif
        data = {'CCCO[O]': [2, ['H', 'C', 'C', 'C', 'O', 'O'], 3],
                'C=C': [1, ['X', 'X', 'X'], 12],
                'S=S': [1, ['S', 'S'], 2],
                }

        for smi in data:
            mult = data[smi][0]
            motif = data[smi][1]
            exp = data[smi][2]
            st_pt = StationaryPoint(smi, 0, mult, smiles=smi)
            st_pt.characterize()
            bond = st_pt.bond
            natom = st_pt.natom
            atom = st_pt.atom
            # do not use any equivalencies
            eqv = [[k] for k in range(natom)]
            hits = find_motif.start_motif(motif, natom, bond, atom, -1, eqv)
            cal = len(hits)
            warn = 'Unexpected number of motif hits for '
            warn += '{}, expected {}, calculated {}'.format(smi, exp, cal)
            self.assertEqual(exp, cal, warn)

    def testNumberOfHitsWithStartAtom(self):
        """
        Test the molecular formula creator, which is
        basically a RDKit feature, but used in KinBot
        for postprocessing
        """
        # the data contains the smiles as key and as value:
        # 1. the multiplicity
        # 2. start atom
        # 3. the motif which will be searched for
        # 4. the expected number of hits of the motif
        data = {'CCCO[O]': [2, 0, ['C', 'C', 'C', 'O', 'O'], 1],
                'C=C': [1, 0, ['X', 'X', 'X'], 2],
                'S=S': [1, 0, ['S', 'S'], 1],
                }

        for smi in data:
            mult = data[smi][0]
            start = data[smi][1]
            motif = data[smi][2]
            exp = data[smi][3]
            st_pt = StationaryPoint(smi, 0, mult, smiles=smi)
            st_pt.characterize()
            bond = st_pt.bond
            natom = st_pt.natom
            atom = st_pt.atom
            # do not use any equivalencies
            eqv = [[k] for k in range(natom)]
            hits = find_motif.start_motif(motif, natom, bond, atom, start, eqv)
            cal = len(hits)
            warn = 'Unexpected number of motif hits for '
            warn += '{}, expected {}, calculated {}'.format(smi, exp, cal)
            self.assertEqual(exp, cal, warn)

    def testNumberOfHitsWithEquivalences(self):
        """
        Test the molecular formula creator, which is
        basically a RDKit feature, but used in KinBot
        for postprocessing
        """
        # the data contains the smiles as key and as value:
        # 1. the multiplicity
        # 2. start atom
        # 3. the motif which will be searched for
        # 4. the expected number of hits of the motif
        data = {'CCCO[O]': [2, ['H', 'C', 'C', 'C', 'O', 'O'], 1],
                'CC': [1, ['C', 'C'], 1],
                }

        for smi in data:
            mult = data[smi][0]
            motif = data[smi][1]
            exp = data[smi][2]
            st_pt = StationaryPoint(smi, 0, mult, smiles=smi)
            st_pt.characterize()
            bond = st_pt.bond
            natom = st_pt.natom
            atom = st_pt.atom
            # do not use any equivalencies
            eqv = st_pt.atom_eqv
            hits = find_motif.start_motif(motif, natom, bond, atom, -1, eqv)
            cal = len(hits)
            warn = 'Unexpected number of motif hits for '
            warn += '{}, expected {}, calculated {}'.format(smi, exp, cal)
            self.assertEqual(exp, cal, warn)

    def testBondFilter(self):
        """
        Test the molecular formula creator, which is
        basically a RDKit feature, but used in KinBot
        for postprocessing
        """
        smi = 'CCC=CCC'
        motif = ['C', 'C', 'C', 'C']
        bondpattern = [2, 'X', 'X']
        exp = 2
        st_pt = StationaryPoint(smi, 0, 1, smiles=smi)
        st_pt.characterize()
        bond = st_pt.bond
        natom = st_pt.natom
        atom = st_pt.atom
        # do not use any equivalencies
        eqv = [[k] for k in range(natom)]
        hits = find_motif.start_motif(motif, natom, bond, atom, -1, eqv)
        count = 0
        for hit in hits:
            if find_motif.bondfilter(hit, bond, bondpattern) == 0:
                count += 1
        warn = 'Unexpected number of motif hits for '
        warn += '{}, expected {}, calculated {}'.format(smi, exp, count)
        self.assertEqual(exp, count, warn)


if __name__ == "__main__":
    unittest.main()
