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

# Test class for the Needleman-Wunsch algorithm
# All test cases are written with PyUnit: http://pyunit.sourceforge.net/

import unittest
import os, sys
lib_path = os.path.abspath('../../')
sys.path.append(lib_path)

from pairwise import NeedlemanWunsch as nw
from helper import PairwiseAlignmentHelper as pah

class NeedlemanWunschTestClass(unittest.TestCase):
    """Class to test the correctness of the computation for the class NeedlemanWunsch."""
    def test_computeMatrix(self):
        """Test of the computation of the matrix."""
        a = "AGC"
        b = "AC"
        computedMatrix = [[0 for i in range(len(b)+1) ] for j in range(len(a)+1)]
       
        # initalize matrix
        for i in range(1, len(a)+1):
            computedMatrix[i][0] = computedMatrix[i-1][0] + pah().weightFunctionDifference("", a[i-1])
        for i in range(1, len(b)+1):
            computedMatrix[0][i] = computedMatrix[0][i-1] + pah().weightFunctionDifference("", b[i-1])
       
        # define values that should be computed by Needleman-Wunsch algorithm
        computedMatrix[1][1] = 0
        computedMatrix[2][1] = 1
        computedMatrix[3][1] = 2
        
        computedMatrix[1][2] = 1
        computedMatrix[2][2] = 1
        computedMatrix[3][2] = 1

        # check if the values computed by Needleman-Wunsch are correct
        self.assertEqual(computedMatrix, nw().compute_matrix(a, b, pah().weightFunctionDifference))
        
    def test_traceback(self):
        """Test of the traceback computation."""
        # test case with a single traceback
        a = "AGC"
        b = "AC"
        computedAlignment = [["AGC", "A-C"]]
        computedMatrix = nw().compute_matrix(a, b, pah().weightFunctionDifference)
        self.assertEqual(computedAlignment, 
                nw().traceback(a, b, computedMatrix,pah().weightFunctionDifference))

        # test case with a multiple traceback
        a = "AT"
        b = "AAGT"
        computedMatrix = nw().compute_matrix(a, b, pah().weightFunctionDifference)
        computedAlignment = [["A--T","AAGT"], ["-A-T","AAGT"]]
        self.assertEqual(computedAlignment, 
                nw().traceback(a, b, computedMatrix,pah().weightFunctionDifference))

if __name__ == "__main__":
    unittest.main() # run all tests     
