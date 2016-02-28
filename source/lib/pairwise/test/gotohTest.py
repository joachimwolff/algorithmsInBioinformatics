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
# Gotoh test class
import unittest
import os, sys
lib_path = os.path.abspath('../../')
sys.path.append(lib_path)

from pairwise import Gotoh 
from helper import PairwiseAlignmentHelper as pah
from helper import MathHelper as mathHelper
class GotohTestClass(unittest.TestCase):
    """Test class to test the correct computation of the Gotoh algorithm."""
    def test_computeMatrix(self):
        """Test method to test the correct computation of the matrix."""
        a = "AGC"
        b = "AC"
        computedMatrixD = [[0 for i in range(len(b)+1) ] for j in range(len(a)+1)]
       	computedMatrixP = [[0 for i in range(len(b)+1) ] for j in range(len(a)+1)]
       	computedMatrixQ = [[0 for i in range(len(b)+1) ] for j in range(len(a)+1)]
        # initalize matrix
        for i in range(1, len(a)+1):
            computedMatrixD[i][0] = pah().gapCost(i)
            computedMatrixP[i][0] = mathHelper.NaN
            computedMatrixQ[i][0] = mathHelper.Inf
        for i in range(1, len(b)+1):
            computedMatrixD[0][i] = pah().gapCost(i)
            computedMatrixP[0][i] = mathHelper.Inf
            computedMatrixQ[0][i] = mathHelper.NaN

       
        # define values that should be computed by Gotoh algorithm
        # matrix D
        computedMatrixD[1][1] = 0
        computedMatrixD[2][1] = 3
        computedMatrixD[3][1] = 4
        
        computedMatrixD[1][2] = 3
        computedMatrixD[2][2] = 1
        computedMatrixD[3][2] = 3

        # matrix P
        computedMatrixP[1][1] = 6
        computedMatrixP[2][1] = 3
        computedMatrixP[3][1] = 4
        
        computedMatrixP[1][2] = 7
        computedMatrixP[2][2] = 6
        computedMatrixP[3][2] = 4

        # matrix Q
        computedMatrixQ[1][1] = 6
        computedMatrixQ[2][1] = 7
        computedMatrixQ[3][1] = 8
        
        computedMatrixQ[1][2] = 3
        computedMatrixQ[2][2] = 6
        computedMatrixQ[3][2] = 7

        computedMatrix = [computedMatrixD, computedMatrixP, computedMatrixQ]
        # print "test: ", computedMatrix
        # check if the values computed by Gotoh are correct
        gotoh = Gotoh(a, b, "weightFunctionDifference", "gapCost")
        gotoh.compute_matrix()
        # print gotoh.computationMatrix
        self.assertEqual(computedMatrix, gotoh.computationMatrix)
        
    def test_traceback(self):
        """Test method to test the correct computation of the traceback."""
        #test case with a single traceback
        a = "AGC"
        b = "AC"
        gotoh = Gotoh(a, b, "weightFunctionDifference", "gapCost")
        computedAlignment = [["AGC", "A-C"]]
        gotoh.compute_matrix()
        gotoh.traceback()
        self.assertEqual(computedAlignment, gotoh.computedAlignment)

        # test case with a multiple traceback
        a = "CC"
        b = "ACCT"
        gotoh2 = Gotoh(a, b, "weightFunctionDifference", "gapCost")
        gotoh2.compute_matrix()
        computedAlignment = [["--CC","ACCT"], ["CC--","ACCT"]]
        gotoh2.traceback()
        self.assertEqual(computedAlignment, gotoh2.computedAlignment)

if __name__ == "__main__":
    unittest.main() # run all tests     
