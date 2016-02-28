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
# Nussinov algorithm

class Nussinov():
    """The algorithm of Nussinov is a RNA secondary structure folding algorithm. It was developed by Ruth Nussinov et al.
    and was published in 1978:
            Nussinov, Ruth, et al. "Algorithms for loop matchings."
            SIAM Journal on Applied mathematics 35.1 (1978): 68-82.
            http://rci.rutgers.edu/~piecze/GriggsNussinovKleitmanPieczenik.pdf
        """
    def __init__(self, rnaSequence):
        """rnaSequence: The RNA sequence for which the folding should be computed."""
        self.sequence = rnaSequence
        self.pairedBases = {}
        self.computationMatrix = [[]]

    def computeMatrix(self):
        """This function computes the matrix which the Nussinov-algorithm is based on."""
        self.computationMatrix = [[0 for i in range(len(self.sequence)+1) ] for j in range(len(self.sequence))]
        i = 2
        while i <= len(self.sequence):
            k = i
            j = 0
            while j <= (len(self.sequence)-2) and k <= (len(self.sequence)):
                self.computeMatrixCell(j, k)
                j += 1
                k += 1
            i += 1

    def computeMatrixCell(self, i, j):
        """This function computes the value for every cell of the matrix for the Nussinov-algorithm.
            i:  First index of cell of the Nussinov-matrix
            j:  Second index of cell of the Nussinov-matrix
        Every cell is the maximum of:
                            |       N_(i, j-1)
            N_(i,j) = max   |max i <= k < j N_(i, k-1) + N_(k+1, j-1) + 1
                            |       S_k and S_j are complementary
        """
        self.computationMatrix[i][j-1]
        maximumValue = [0,0,0]
        k = i
        while i <= k and k < j:
            if self.complementary(self.sequence[k], self.sequence[j-1]):
                pairingValue = self.computationMatrix[i][k-1] + self.computationMatrix[k+1][j-1] + 1
                if maximumValue[2] < pairingValue:
                    maximumValue[0] = k
                    maximumValue[1] = j
                    maximumValue[2] = pairingValue
            k += 1
        self.computationMatrix[i][j] = max(self.computationMatrix[i][j-1], maximumValue[2])

    def complementary(self, characterA, characterB):
        """Returns True if two RNA nucleotides are complementary, False otherwise.
        Nucleotides are complemetary if there are "A" and "U" or "C" and "G".
            characterA: First nucleotide
            characterB: Second nucleotide"""
        if characterA == "A" and characterB == "U":
            return True
        elif characterA == "U" and characterB == "A":
            return True
        elif characterA == "C" and characterB == "G":
            return True
        elif characterA == "G" and characterB == "C":
            return True
        return False

    def traceback(self, i, j):
        """Computes the traceback for the Nussinov-algorithm.
            i:  First index of cell of the Nussinov-matrix
            j:  Second index of cell of the Nussinov-matrix
            """
        if j <= i:
            return
        elif self.computationMatrix[i][j] == self.computationMatrix[i][j-1]:
            self.traceback(i, j-1)
            return
        else:
            k = i
            while i <= k and k < j:
                if self.complementary(self.sequence[k-1], self.sequence[j-1]):

                    if self.computationMatrix[i][j] == self.computationMatrix[i][k-1] + self.computationMatrix[k][j-1] + 1:
                        self.pairedBases[k] = j
                        self.traceback(i, k-1)
                        self.traceback(k, j -1)
                        return
                k += 1

    def execute(self):
        """To compute the Nussinov-algorithm execute this method. It returns a dictionary with the paired bases."""
        self.computeMatrix()
        self.traceback(0, len(self.sequence))
        print self.pairedBases
        print len(self.pairedBases)
        return self.pairedBases