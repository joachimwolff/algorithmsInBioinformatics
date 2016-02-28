#! /usr/bin/python
# Copyright 2015 Joachim Wolff
# Programming Course: Algorithms in Bioinformatics
# Tutors: Robert Kleinkauf, Omer Alkhnbashi
# Winter semester 2014/2015
#
# Chair of Bioinformatics
# Department of Computer Science
# Faculty of Engineering
# Albert-Ludwig-University Freiburg im Breisgau
#
# UPGMA/WPGMA test class
import unittest
import os, sys

lib_path = os.path.abspath('../../')
sys.path.append(lib_path)

from multiple import UpgmaWpgma


class UpgmaWpgmaTestClass(unittest.TestCase):
    """Test class to test the correct computation of the UPGMA/WPGMA algorithm."""

    def test_computeMinimalDistance(self):
        distanceDictionary = {"0 1": 1, "0 2": 2, "0 3": 3, "1 2": 2, "1 3": 3, "1 4": 3}
        upgma = UpgmaWpgma(distanceDictionary, 4)
        expectedValue = ["0 1", 1]
        self.assertEqual(expectedValue, upgma.compute_minimal_distance())

    def test_computeClustering(self):
        distanceDictionary = {"0 1": 1, "0 2": 2, "0 3": 3, "1 2": 2, "1 3": 3, "2 3": 3}
        upgma = UpgmaWpgma(distanceDictionary, 4)
        expectedValue = {"0 1": 4, "2 4": 5, "3 5": 6}
        upgma.compute_clustering()
        print upgma.get_newick_tree()
        self.assertEqual(expectedValue, upgma.mapping)

        print upgma.get_newick_tree(with_edge_weights=True)
        distanceDictionary = {"0 1": 6, "0 2": 10, "0 3": 10, "0 4": 10, "1 2": 10, "1 3": 10, "1 4": 10, "2 3": 2,
                              "2 4": 6, "3 4": 6}
        upgma2 = UpgmaWpgma(distanceDictionary, 5)
        expectedValue = {"2 3": 5, "0 1": 7, "4 5": 6, "6 7": 8}
        upgma2.compute_clustering()
        print upgma2.get_newick_tree(with_edge_weights=False)
        self.assertEqual(expectedValue, upgma2.mapping)
        print upgma2.get_newick_tree(with_edge_weights=True)




    def test_getNewickTree(self):
        mapping = {'1 3': 5, '4 6': 7, '5 7': 8, '0 2': 6}
        distanceDictionary = {}
        upgma = UpgmaWpgma(distanceDictionary, 5)
        upgma.mapping = mapping
        upgma.get_newick_tree()


if __name__ == "__main__":
    unittest.main()  # run all tests
