#! /usr/bin/python
# Copyright 2014 Joachim Wolff
# Programming Course: Algorithms in Bioinformatics
# Tutors: Robert Kleinkauf, Omer Alkhnbashi
# Winter semester 2014/2015
#
# Chair of Bioinformatics
# Department of Computer Science
# Faculty of Engineering
# Albert-Ludwig-University Freiburg im Breisgau
#
# Feng-Doolittle test class
import unittest
import os, sys
lib_path = os.path.abspath('../../')
sys.path.append(lib_path)

from multiple import FengDoolittle

class FengDoolittleTestClass(unittest.TestCase):
    """Test class to test the correct computation of the Needleman-Wunsch n=3 algorithm."""
    def test_computeAlignments(self):
    	sequences = ["ACTG", "AT", "ACG"]
        expectedAlignments = [[["ACTG", "A-T-"],0,1], [["ACTG", "AC-G"],0,2], [["AT-", "ACG"],1,2]]
        fd = FengDoolittle(sequences, "weightFunctionDifference", "pam250")
        fd.computeAlignments()
        self.assertEqual(expectedAlignments, fd.alignments)
    def test_computeDistanceDictionary(self):
        sequences = ["ACCCAT", "ACGGAT", "AACCT"]
        expectedAlignments = [["AC-CAT", "ACGGAT"], ["ACGGAT", "AACCAT"], ["-ACCAT", "AACCAT"]]
        fd = FengDoolittle(sequences, "weightFunctionDifference", "pam250")
        fd.computeAlignments()
        fd.computeDistanceDictionary()
    def test_computeOrderOfSequencesToAlign(self):
        sequences = ["ACTG", "AT", "ACG"]
        fd = FengDoolittle(sequences, "weightFunctionDifference", "pam250")
        fd.computeAlignments()
        fd.computeDistanceDictionary()
        fd.buildTree()
        # print "NewickTree: ",fd.newickTree
        expectedResult = [[[0],[2]], [[0,2],[1]]]
        # print "asd"
        fd.computeOrderOfSequencesToAlign()
        # print "asd"
        self.assertEqual(expectedResult, fd.orderToAlign)
    def test_computeMultipleAlignment(self):
        sequences = ["ACTG", "AT", "ACG"]
        expectedResult = {0: 'ACTG', 1: 'AXTX', 2: 'ACXG'}
        fd = FengDoolittle(sequences, "weightFunctionDifference", "pam250")
        fd.computeMultipleAlignment()
        self.assertEqual(expectedResult, fd.computeMultipleAlignment())

        sequences = ["ACCAT", "ACGGAT", "AACCAT"]
        expectedResult = {0: 'AXCCXAT', 1: 'AXCGGAT', 2: 'AACCXAT'}
        fd2 = FengDoolittle(sequences, "weightFunctionDifference", "pam250")
        fd2.computeMultipleAlignment()
        self.assertEqual(expectedResult, fd2.computeMultipleAlignment())
if __name__ == "__main__":
    unittest.main() # run all tests     
