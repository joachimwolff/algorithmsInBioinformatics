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
# Gotoh algorithm
from helper import PairwiseAlignmentHelper as pah
from helper import MathHelper as mathHelper
import sys

class Gotoh():
    """This class holds methods which are needed to compute the pairwise 
    alignment algorithm from Osamu Gotoh, published in 1982:
        Osamu Gotoh (1982). "An improved algorithm for matching biological sequences".
        Journal of molecular biology 162: 705.
        https://www.cs.umd.edu/class/spring2003/cmsc838t/papers/gotoh1982.pdf
    """
    def __init__(self, sequenceA, sequenceB, scoreFunction, costFunction):
        """Initalize all variables and methods needed to compute the Gotoh algorithm.
            sequenceA:      A string with the first DNA sequence.
            sequenceB:      A string with the second DNA sequence.
            scoreFunction:  The name of a weight function as a String which is defined 
                            in the pairwiseAlignmentHelper-class.
            costFunction:   The name of a gap cost function as a String which is defined 
                            in the pairwiseAlignmentHelper-class.
        """
        if scoreFunction in dir(pah) and callable(getattr(pah, scoreFunction)):
            scoreFunctionObj = eval('pah().' + scoreFunction)
        else:
            print "Score function not found!"
            sys.exit()       
        if costFunction in dir(pah) and callable(getattr(pah, costFunction)):
            costFunctionObj = eval('pah().' + costFunction)
        else:
            print "Gap cost function not found!"
            sys.exit() 
        self.computationMatrix = [[],[],[]]
        self.sequenceA = sequenceA
        self.sequenceB = sequenceB
        self.scoreFunction = scoreFunctionObj
        self.costFunction = costFunctionObj
        self.beta = self.costFunction(1) - self.costFunction(0)
        self.i = 0
        self.j = 0
        self.tracebackStack = [[]]
        self.tracebackStackIndex = 0
        self.indiciesStack = [[]]
        self.computedAlignment = []


    def computeMatrix(self):
        """Initalize the three matricies needed for the Gotoh-Algorithm.
        The sequences A and B, the weight function and the gap costs have to be defined 
        by the creation of the object of this class."""
        computationMatrixD = [[0 for i in range(len(self.sequenceB)+1) ] for j in range(len(self.sequenceA)+1)]
        computationMatrixP = [[0 for i in range(len(self.sequenceB)+1) ] for j in range(len(self.sequenceA)+1)]
        computationMatrixQ = [[0 for i in range(len(self.sequenceB)+1) ] for j in range(len(self.sequenceA)+1)]
        # initalize matrix
        for i in range(1, len(self.sequenceA)+1):
            computationMatrixD[i][0] = self.costFunction(i)
            computationMatrixP[i][0] = mathHelper.NaN
            computationMatrixQ[i][0] = mathHelper.Inf
        for i in range(1, len(self.sequenceB)+1):
            computationMatrixD[0][i] = self.costFunction(i)
            computationMatrixP[0][i] = mathHelper.Inf
            computationMatrixQ[0][i] = mathHelper.NaN

        for i in range(1, len(self.sequenceA)+1):
            for j in range(1, len(self.sequenceB)+1):
                computationMatrixP[i][j] = self.computeP(computationMatrixD[i-1][j], computationMatrixP[i-1][j], self.costFunction, self.beta)  
                computationMatrixQ[i][j] = self.computeQ(computationMatrixD[i][j-1], computationMatrixQ[i][j-1], self.costFunction, self.beta)
                computationMatrixD[i][j] = self.computeD(computationMatrixD[i-1][j-1], computationMatrixP[i][j], computationMatrixQ[i][j], self.sequenceA[i-1], self.sequenceB[j-1], self.scoreFunction)
        self.computationMatrix = [computationMatrixD, computationMatrixP, computationMatrixQ]

    def computeP(self, valueOfD, valueOfP, costFunction, beta):
        """Compute the values for matrix P.
            This is the minimum value of:
                matrix D of cell (i-1, j) + gap costs
                and 
                matrix P of cell (i-1, j) + 1
            valueOfD:       The value from matrix D of cell i-1, j.
            valueOfP:       The value from matrix P of cell i-1, j.
            costFunction:   The gap cost function defined at the object creation.
            beta:           The beta value from the gap costs."""
        return min(valueOfD + costFunction(1), valueOfP + beta)
      
    def computeQ(self, valueOfD, valueOfQ, costFunction, beta):
        """Compute the values for matrix Q.
            This is the minimum value of:
                matrix D of cell (i, j-1) + gap costs
                and 
                matrix Q of cell (i, j-1) + 1
            valueOfD:       The value from matrix D of cell i, j-1.
            valueOfQ:       The value from matrix Q of cell i, j-1.
            costFunction:   The gap cost function defined at the object creation.
            beta:           The beta value from the gap costs."""
        return min(valueOfD + costFunction(1), valueOfQ + beta)

    def computeD(self, valueOfD, valueOfP, valueOfQ, characterA, characterB, scoreFunction):
        """Compute the values for matrix D.
            This is the minimum value of:
                matrix D of cell (i-1, j-1) + w(a,b)
                and 
                matrix P of cell (i, j)
                and 
                matrix Q of cell (i, j)
            valueOfD:       The value from matrix D of cell i-1, j-1.
            valueOfP:       The value from matrix P of cell i, j.
            valueOfQ:       The value from matrix Q of cell i, j.
            characterA:     The character in sequence A at position i.
            characterB:     The character in sequence B at position j.
            scoreFunction:  The weight cost function defined at the object creation."""
        return min(valueOfP, min(valueOfQ, valueOfD + scoreFunction(characterA, characterB)))

    def traceback(self):
        """Computes the traceback for the Gotoh algorithm."""
        self.j = len(self.computationMatrix[0][0]) - 1
        self.i = len(self.computationMatrix[0]) - 1
        self.tracebackStackIndex = 0
        self.indiciesStack[self.tracebackStackIndex] = [self.i, self.j, pah.matrixIndexD]
        tracebackDone = False
        while not tracebackDone:
            while self.i > 0 or self.j > 0:
                if self.indiciesStack[self.tracebackStackIndex][2] == pah.matrixIndexD:
                    self.tracebackD()
                elif self.indiciesStack[self.tracebackStackIndex][2] == pah.matrixIndexP:
                    self.tracebackP()
                elif self.indiciesStack[self.tracebackStackIndex][2] == pah.matrixIndexQ:
                    self.tracebackQ()
                self.i = self.indiciesStack[self.tracebackStackIndex][0]
                self.j = self.indiciesStack[self.tracebackStackIndex][1]
            tracebackDone = True
            for i in range(0, len(self.indiciesStack)):
                if self.indiciesStack[i][0] > 0 or self.indiciesStack[i][1] > 0:
                    self.tracebackStackIndex = i
                    tracebackDone = False
                    break
            self.i = self.indiciesStack[self.tracebackStackIndex][0]
            self.j = self.indiciesStack[self.tracebackStackIndex][1]
        for i in range(0, len(self.tracebackStack)):
            self.computedAlignment.append(self.buildAlignment(self.tracebackStack[i]))

    
    def tracebackD(self):
        """Computes the traceback for a cell of the matrix D."""
        a = self.sequenceA[self.i - 1]
        b = self.sequenceB[self.j - 1]
        split = 0
        pathVariableI = self.i
        pathVariableJ = self.j
        if self.j > 0 and self.i > 0:
            if self.computationMatrix[pah.matrixIndexD][self.i][self.j] == self.computationMatrix[pah.matrixIndexD][self.i-1][self.j-1] + self.scoreFunction(a,b):
                self.tracebackStack[self.tracebackStackIndex].append(pah.diagonalD)
                pathVariableI -= 1
                pathVariableJ -= 1
                split = 1
            if self.computationMatrix[pah.matrixIndexD][self.i][self.j] == self.computationMatrix[pah.matrixIndexQ][self.i][self.j]:
                if split == 0:
                    self.tracebackStack[self.tracebackStackIndex].append(pah.dotQ)
                    self.indiciesStack[self.tracebackStackIndex][2] = pah.matrixIndexQ
                    split = 1
                else:
                    self.tracebackStack.append(self.tracebackStack[self.tracebackStackIndex][0:-1])
                    self.tracebackStack[len(self.tracebackStack)-1].append(pah.dotQ)
                    self.indiciesStack.append([self.i,self.j, pah.matrixIndexQ])
            if self.computationMatrix[pah.matrixIndexD][self.i][self.j] == self.computationMatrix[pah.matrixIndexP][self.i][self.j]:
                if split == 0: 
                    self.tracebackStack[self.tracebackStackIndex].append(pah.dotP)
                    self.indiciesStack[self.tracebackStackIndex][2] = pah.matrixIndexP
                else:
                    self.tracebackStack.append(self.tracebackStack[self.tracebackStackIndex][0:-1])
                    self.tracebackStack[len(self.tracebackStack)-1].append(pah.dotP)
                    self.indiciesStack.append([self.i, self.j, pah.matrixIndexP])

        if self.i == 0:
            self.tracebackStack[self.tracebackStackIndex].append(pah.leftD)
            pathVariableJ -= 1
        if self.j == 0:
            self.tracebackStack[self.tracebackStackIndex].append(pah.upD)
            pathVariableI -= 1
        if self.i <= 0 or pathVariableI <= 0:
            pathVariableI = 0
        if self.j <= 0 or pathVariableJ <= 0:
            pathVariableJ = 0
        self.indiciesStack[self.tracebackStackIndex][0] = pathVariableI
        self.indiciesStack[self.tracebackStackIndex][1] = pathVariableJ


    def tracebackP(self):
        """Computes the traceback for a cell of the matrix P"""
        split = False
        if self.i > 0:
            if self.computationMatrix[pah.matrixIndexP][self.i][self.j] == self.computationMatrix[pah.matrixIndexD][self.i-1][self.j] + self.costFunction(1):
                self.tracebackStack[self.tracebackStackIndex].append(pah.upD)
                self.indiciesStack[self.tracebackStackIndex][0] -= 1
                self.indiciesStack[self.tracebackStackIndex][2] = pah.matrixIndexD
                split = True
            if self.computationMatrix[pah.matrixIndexP][self.i][self.j] == self.computationMatrix[pah.matrixIndexP][self.i-1][self.j] + self.beta:
                if split:
                    self.tracebackStack.append(self.tracebackStack[self.tracebackStackIndex][0:-1])
                    self.tracebackStack[len(self.tracebackStack)-1].append(pah.upP)
                    self.indiciesStack.append([self.i - 1, self.j, pah.matrixIndexP])
                else:
                    self.tracebackStack[self.tracebackStackIndex].append(pah.upP)
                    self.indiciesStack[self.tracebackStackIndex][0] -= 1
                    self.indiciesStack[self.tracebackStackIndex][2] = pah.matrixIndexP

    def tracebackQ(self):
        """Computes the traceback for a cell of the matrix Q"""
        split = False
        if self.j > 0:
            if self.computationMatrix[pah.matrixIndexQ][self.i][self.j] == self.computationMatrix[pah.matrixIndexD][self.i][self.j-1] + self.costFunction(1):
                self.tracebackStack[self.tracebackStackIndex].append(pah.leftD)
                self.indiciesStack[self.tracebackStackIndex][1] -= 1
                self.indiciesStack[self.tracebackStackIndex][2] = pah.matrixIndexD
                split = True

            if self.computationMatrix[pah.matrixIndexQ][self.i][self.j] == self.computationMatrix[pah.matrixIndexQ][self.i][self.j-1] + self.beta:
                if split:
                    self.tracebackStack.append(self.tracebackStack[self.tracebackStackIndex][0:-1])
                    self.tracebackStack[len(self.tracebackStack)-1].append(pah.leftQ)
                    self.indiciesStack.append([self.i , self.j - 1, pah.matrixIndexQ])
                else:
                    self.tracebackStack[self.tracebackStackIndex].append(pah.leftQ)
                    self.indiciesStack[self.tracebackStackIndex][1] -= 1
                    self.indiciesStack[self.tracebackStackIndex][2] = pah.matrixIndexQ

    def buildAlignment(self, tracebackStack):
        """A method to compute the alignment of a given traceback of the Gotoh algorithm.
            tracebackStack: The computed traceback path for one alignment as a list."""
        i = 0
        j = 0
        k = len(tracebackStack)-1
        alignmentOfA = ""
        alignmentOfB = ""
        while len(tracebackStack) > 0:
            try:
                tracebackElement = tracebackStack.pop(k)
                if pah.leftQ == tracebackElement or pah.leftD == tracebackElement:
                    alignmentOfA += "-"
                    alignmentOfB += self.sequenceB[j]
                    j += 1
                elif pah.upP == tracebackElement or pah.upD == tracebackElement:
                    alignmentOfA += self.sequenceA[i]
                    alignmentOfB += "-"
                    i += 1
                elif pah.diagonalD == tracebackElement:
                    alignmentOfA += self.sequenceA[i]
                    alignmentOfB += self.sequenceB[j]
                    i += 1
                    j += 1
                k -= 1

            except:
                print "An error occured."
                sys.exit()

        while i < len(self.sequenceA):
            alignmentOfA += self.sequenceA[i]
            i += 1
        while j < len(self.sequenceB):
            alignmentOfB += self.sequenceB[j]
            j += 1   
        alignment = [alignmentOfA, alignmentOfB]
        return alignment

    def compute(self):
        """Method to start the computation of the Gotoh algorithm."""
        self.computeMatrix()
        self.traceback()
        return self.computedAlignment

