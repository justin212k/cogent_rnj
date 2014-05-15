#! /usr/bin/env python
import unittest, os

from cogent.phylo.distance import *
from rnj import rnj
from cogent import LoadTree

__author__ = "Justin Kuczynski"
__copyright__ = ""
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

class TreeReconstructionTests(unittest.TestCase):
    def setUp(self):
        pass
        
    def _test_tree(self, method, treestring):
        t = LoadTree(treestring=treestring)
        t_distances = t.getDistances()
        reconstructed = method(t_distances)
        distances = reconstructed.getDistances()
        for key in t_distances:
            self.assertAlmostEqual(t_distances[key], distances[key])

    def _test_phylo_method(self, method):
        """testing (well, exercising at least), rnj"""
        self._test_tree(method, '((a:3,b:4):20,(c:6,d:7):30,e:5)')
        self._test_tree(method, '((a:3,b:4):0,(c:6,d:7):30,e:5)')
        self._test_tree(method, '((a:3,b:4,c:6,d:7):30,e:5)')

    def test_rnj(self):
        """testing (well, exercising at least), rnj"""
        self._test_phylo_method(rnj)

if __name__ == '__main__':
    unittest.main()
