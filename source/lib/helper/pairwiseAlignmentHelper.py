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

class PairwiseAlignmentHelper():
    """Class to support the pairwise alignment algorithms Needleman-Wunsch and Gotoh."""
    diagonalD = 0
    dotQ = 1
    dotP = 2
    upD = 3
    upP = 4
    leftD = 5
    leftQ = 6
    matrixIndexD = 0
    matrixIndexP = 1
    matrixIndexQ = 2
    left = 0
    up = 1
    diagonal = 2

    def weightFunctionDifference(self, a, b):
        """Weight function with 0 if a==b and 1 else."""
        if a == b:
            return 0
        elif a != b:
            return 1

    def gapCost(self, x):
        """Returns a gap cost of g(x) = 2 + k."""
        return 2 + x

    def pam250(self, a, b):
        """Returns the value of an amino acid given a pam250 matrix. If it is a gap, 1 is returned.
            Source: http://www.icp.ucl.ac.be/~opperd/private/pam250.html"""
        pam250 = [[13, 6, 9, 9, 5, 8, 9, 12, 6, 8, 6, 7, 7, 4, 11, 11, 11, 2, 4, 9]
            , [3, 17, 4, 3, 2, 5, 3, 2, 6, 3, 2, 9, 4, 1, 4, 4, 3, 7, 2, 2]
            , [4, 4, 6, 7, 2, 5, 6, 4, 6, 3, 2, 5, 3, 2, 4, 5, 4, 2, 3, 3]
            , [5, 4, 8, 11, 1, 7, 10, 5, 6, 3, 2, 5, 3, 1, 4, 5, 5, 1, 2, 3]
            , [2, 1, 1, 1, 52, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 3, 2, 1, 4, 2]
        , [3, 5, 5, 6, 1, 10, 7, 3, 7, 2, 3, 5, 3, 1, 4, 3, 3, 1, 2, 3]
        , [5, 4, 7, 11, 1, 9, 12, 5, 6, 3, 2, 5, 3, 1, 4, 5, 5, 1, 2, 3]
        , [12, 5, 10, 10, 4, 7, 9, 27, 5, 5, 4, 6, 5, 3, 8, 11, 9, 2, 3, 7]
        , [2, 5, 5, 4, 2, 7, 4, 2, 15, 2, 2, 3, 2, 2, 3, 3, 2, 2, 3, 2]
        , [3, 2, 2, 2, 2, 2, 2, 2, 2, 10, 6, 2, 6, 5, 2, 3, 4, 1, 3, 9]
        , [6, 4, 4, 3, 2, 6, 4, 3, 5, 15, 34, 4, 20, 13, 5, 4, 6, 6, 7, 13]
        , [6, 18, 10, 8, 2, 10, 8, 5, 8, 5, 4, 24, 9, 2, 6, 8, 8, 4, 3, 5]
        , [1, 1, 1, 1, 0, 1, 1, 1, 1, 2, 3, 2, 6, 2, 1, 1, 1, 1, 1, 2]
        , [2, 1, 2, 1, 1, 1, 1, 1, 3, 5, 6, 1, 4, 32, 1, 2, 2, 4, 20, 3]
        , [7, 5, 5, 4, 3, 5, 4, 5, 5, 3, 3, 4, 3, 2, 20, 6, 5, 1, 2, 4]
        , [9, 6, 8, 7, 7, 6, 7, 9, 6, 5, 4, 7, 5, 3, 9, 10, 9, 4, 4, 6]
        , [8, 5, 6, 6, 4, 5, 5, 6, 4, 6, 4, 6, 5, 3, 6, 8, 11, 2, 3, 6]
        , [0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 55, 1, 0]
        , [1, 1, 2, 1, 3, 1, 1, 1, 3, 2, 2, 1, 2, 15, 1, 2, 2, 3, 31, 2]
        , [7, 4, 4, 4, 4, 4, 4, 4, 5, 4, 15, 10, 4, 10, 5, 5, 5, 72, 4, 17]]

        pamdict = {"A":0, "R":1, "N":2, "D":3, "C":4, "Q":5, "E":6, "G":7, "H":8, "I":9, "L":10, "K":11, "M":12,
                   "F":13, "P":14, "S":15, "T":16, "W":17, "Y":18, "V":19}

        if a in pamdict and b in pamdict:
            return pam250[pamdict[a]][pamdict[b]]
        else:
            return 1



