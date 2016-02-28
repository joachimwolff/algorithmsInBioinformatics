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
# Needleman-Wunsch with n=3 test class
import unittest
import os, sys
lib_path = os.path.abspath('../../')
sys.path.append(lib_path)

from multiple import NeedlemanWunschN3 as nw
from helper import MultipleAlignmentHelper as mah
from helper import MathHelper as mathHelper
class NeelemanWunschN3TestClass(unittest.TestCase):
    """Test class to test the correct computation of the Needleman-Wunsch n=3 algorithm."""
    def test_computeMatrix(self):
        sequenceA = "AC"
        sequenceB = "AGT"
        sequenceC = "AGT"
        expectedMatrix = [[[0 for i in range(len(sequenceC)+1) ] for j in range(len(sequenceB)+1)] for k in range(len(sequenceA)+1 )]
        for i in range(1, len(sequenceA)+1):
            expectedMatrix[i][0][0] = expectedMatrix[i-1][0][0] + mah().weightFunctionDifference("", "", sequenceA[i-1])
        for i in range(1, len(sequenceB)+1):
            expectedMatrix[0][i][0] = expectedMatrix[0][i-1][0] + mah().weightFunctionDifference("", "", sequenceB[i-1])
        for i in range(1, len(sequenceC)+1):
            expectedMatrix[0][0][i] = expectedMatrix[0][0][i-1] + mah().weightFunctionDifference("", "", sequenceC[i-1])
        for i in range(1, len(sequenceA)+1):
            for j in range(1, len(sequenceB)+1):
                expectedMatrix[i][j][0] = expectedMatrix[i-1][j-1][0] + mah().weightFunctionDifference(sequenceA[i-1], sequenceB[j-1], "")
        for i in range(1, len(sequenceA)+1):
            for k in range(1, len(sequenceC)+1):
                expectedMatrix[i][0][k] = expectedMatrix[i-1][0][k-1] + mah().weightFunctionDifference(sequenceA[i-1], "", sequenceC[k-1])
        for j in range(1, len(sequenceB)+1):
            for k in range(1, len(sequenceC)+1):
                expectedMatrix[0][j][k] = expectedMatrix[0][j-1][k-1] + mah().weightFunctionDifference("", sequenceB[j-1], sequenceC[k-1])

        assertEqual()

if __name__ == "__main__":
    unittest.main() # run all tests     
