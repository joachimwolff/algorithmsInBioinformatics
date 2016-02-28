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

from structurePrediction import Nussinov

class NussinovTestClass(unittest.TestCase):
	def test_computeMatrix(self):
		# example for the slides of Prof. Backofen
		expectedMatrix = [[0,0,1,1,1,2,2,2,3], [0,0,0,0,0,1,1,1,2], [0,0,0,0,0,1,1,1,2], [0,0,0,0,0,1,1,1,2], [0,0,0,0,0,0,0,1,1], [0,0,0,0,0,0,0,0,1], [0,0,0,0,0,0,0,0,1],[0,0,0,0,0,0,0,0,0]]
		rnaSequence = "GCACGACG"
		nussinov = Nussinov(rnaSequence)
		nussinov.compute_matrix()
		self.assertEqual(expectedMatrix, nussinov.computationMatrix)

	def test_traceback(self):
		# example for the slides of Prof. Backofen
		expectedMatrix = {1:2, 4:8, 5:7}
		rnaSequence = "GCACGACG"
		nussinov = Nussinov(rnaSequence)
		nussinov.compute_matrix()
		nussinov.traceback(0, len(rnaSequence))
		self.assertEqual(expectedMatrix, nussinov.pairedBases)

if __name__ == "__main__":
    unittest.main() # run all tests     