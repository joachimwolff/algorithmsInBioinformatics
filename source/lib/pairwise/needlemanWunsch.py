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
# Needleman-Wunsch algorithm
import sys
from helper import PairwiseAlignmentHelper as pah


class NeedlemanWunsch():
    """This class holds methods which are needed to compute the pairwise
    alignment algorithm from Saul Needleman and Christian Wunsch, published in 1970:
        Needleman, Saul B.; and Wunsch, Christian D. (1070).
        A general method applicable to search for similarities in the aminoacid
        sequence of two proteins. Journal of Molecular Biology 48 (3): 443-53
        http://www.cise.ufl.edu/class/cis4930sp09rab/00052.pdf"""

    def computeMatrix(self, sequenceA, sequenceB, scoreFunction):
        """Initalize and computes the values for the Needleman-Wunsch matrix.
            sequenceA:      A string with the first DNA sequence.
            sequenceB:      A string with the second DNA sequence.
            scoreFunction:  The name of a weight function as a String which is defined 
                            in the pairwiseAlignmentHelper-class."""
        computationMatrix = [[0 for i in range(len(sequenceB) + 1)] for j in range(len(sequenceA) + 1)]

        # initalize matrix
        for i in range(1, len(sequenceA) + 1):
            computationMatrix[i][0] = computationMatrix[i - 1][0] + scoreFunction("", sequenceA[i - 1])
        for i in range(1, len(sequenceB) + 1):
            computationMatrix[0][i] = computationMatrix[0][i - 1] + scoreFunction("", sequenceB[i - 1])

        for i in range(1, len(sequenceA) + 1):
            for j in range(1, len(sequenceB) + 1):
                computationMatrix[i][j] = self.computeMinimum(sequenceA[i - 1], sequenceB[j - 1],
                                                              computationMatrix[i][j - 1], computationMatrix[i - 1][j],
                                                              computationMatrix[i - 1][j - 1], scoreFunction)
        return computationMatrix

    def computeMinimum(self, characterOfA, characterOfB, predecessorLeft, predecessorUp, predecessorDiagonal,
                       scoreFunction):
        """Computes the minimum of a given cell for the Needleman-Wunsch matrix.
                characterA:             The character in sequence A at position i.
                characterB:             The character in sequence B at position j.
                predecessorLeft:        The value i, j-1 in the matrix.
                predecessorUp:          The value i-1, j in the matrix.
                predecessorDiagonal:    The value i-1, j-1 in the matrix.
                scoreFunction:          The weight function defined in 
                                        class pairwiseAlignmentHelper."""
        costUp = predecessorUp + scoreFunction(characterOfA, "")
        costDiagonal = predecessorDiagonal + scoreFunction(characterOfA, characterOfB)
        costLeft = predecessorLeft + scoreFunction("", characterOfB)
        return min(costUp, costDiagonal, costLeft)

    def traceback(self, sequenceA, sequenceB, computationMatrix, scoreFunction, maxOptimalSolutions=-1):
        """Computes the traceback for the Needleman-Wunsch matrix.
                sequenceA:          A string with the first DNA sequence.
                sequenceB:          A string with the second DNA sequence.
                computationMatrix:  The computed matrix for the two sequences.
                scoreFunction:      The name of a weight function as a String which is defined 
                                    in the pairwiseAlignmentHelper-class.
            """
        tracebackStack = [[]]
        indiciesStack = [[len(computationMatrix) - 1, len(computationMatrix[0]) - 1]]
        tracebackCount = 0
        tracebackDone = False
        optimalSolutionsCount = 0
        l = 0
        allTracebacksComputed = 0
        appendTracebackStack = tracebackStack.append
        appendIndices = indiciesStack.append
        while not tracebackDone:

            optimalSolutionsCount += 1
            i = indiciesStack[tracebackCount][0]
            j = indiciesStack[tracebackCount][1]
            split = False
            appendTraceback = tracebackStack[tracebackCount].append

            while i > 0 or j > 0:
                pathVariableI = i
                pathVariableJ = j
                # left arrow
                if j > 0:
                    if computationMatrix[i][j] == computationMatrix[i][j - 1] + scoreFunction("", sequenceB[j - 1]):
                        # tracebackStack[tracebackCount].append(pah.left)
                        appendTraceback(pah.left)
                        pathVariableJ -= 1  # change j
                        split = True

                # up arrow
                if i > 0:
                    if computationMatrix[i][j] == computationMatrix[i - 1][j] + scoreFunction(sequenceA[i - 1], ""):
                        if split == False:
                            appendTraceback(pah.up)
                            # tracebackStack[tracebackCount].append(pah.up)
                            pathVariableI -= 1
                            split = True
                        else:
                            appendTracebackStack(tracebackStack[tracebackCount][0:-1])
                            tracebackStack[len(tracebackStack) - 1].append(pah.up)
                            appendIndices([i - 1, j])

                # diagonal arrow
                if i > 0 and j > 0:
                    if computationMatrix[i][j] == computationMatrix[i - 1][j - 1] + scoreFunction(sequenceA[i - 1],
                                                                                                  sequenceB[j - 1]):
                        if split == False:
                            appendTraceback(pah.diagonal)
                            # tracebackStack[tracebackCount].append(pah.diagonal)
                            pathVariableI -= 1
                            pathVariableJ -= 1
                        elif split == True:
                            appendTracebackStack(tracebackStack[tracebackCount][0:-1])
                            tracebackStack[len(tracebackStack) - 1].append(pah.diagonal)
                            appendIndices([i - 1, j - 1])
                split = 0
                i = pathVariableI
                j = pathVariableJ

            indiciesStack[tracebackCount][0] = i
            indiciesStack[tracebackCount][1] = j
            l = tracebackCount
            while l < len(indiciesStack):
                if indiciesStack[l][0] == 0 and indiciesStack[l][1] == 0:
                    allTracebacksComputed += 1
                else:
                    tracebackCount = l
                    l = len(indiciesStack)
                l += 1
            if allTracebacksComputed >= len(indiciesStack):
                tracebackDone = True
            if maxOptimalSolutions != -1 and optimalSolutionsCount >= maxOptimalSolutions:
                tracebackDone = True
                # allTracebacksComputed = 0

        computedAlignment = []
        if maxOptimalSolutions == -1:
            for i in range(0, len(tracebackStack)):
                computedAlignment.append(self.buildAlignment(tracebackStack[i], sequenceA, sequenceB))
        else:
            for i in range(0, maxOptimalSolutions):
                computedAlignment.append(self.buildAlignment(tracebackStack[i], sequenceA, sequenceB))
        return computedAlignment

    def buildAlignment(self, tracebackStack, sequenceA, sequenceB):
        """Builds the alignment for one traceback path.
                tracebackStack: The computed tracebackpath as a list = []
                sequenceA:      A string with the first DNA sequence.
                sequenceB:      A string with the second DNA sequence.
            """
        i = 0
        j = 0
        k = len(tracebackStack) - 1
        alignmentOfA = ""
        alignmentOfB = ""

        while len(tracebackStack) > 0:
            try:
                tracebackElement = tracebackStack.pop(k)
                if pah.left == tracebackElement:
                    alignmentOfA += "-"
                    alignmentOfB += sequenceB[j]
                    j += 1
                elif pah.up == tracebackElement:
                    alignmentOfA += sequenceA[i]
                    alignmentOfB += "-"
                    i += 1
                elif pah.diagonal == tracebackElement:
                    alignmentOfA += sequenceA[i]
                    alignmentOfB += sequenceB[j]
                    i += 1
                    j += 1
                k -= 1
            except:
                print "An error occured."
                sys.exit()
        while i < len(sequenceA):
            alignmentOfA += sequenceA[i]
            i += 1
        while j < len(sequenceB):
            alignmentOfB += sequenceB[j]
            j += 1
        alignment = [alignmentOfA, alignmentOfB]
        return alignment

    def compute(self, sequences, scoreFunction, maxOptimalSolutions=-1):
        """Method to execute the Needleman-Wunsch algorithm.
            sequences:      A list with two strings which represents the DNA sequences.
            scoreFunction:  The name of the weight function defined in 
                            class pairwiseAlignmentHelper."""
        if scoreFunction in dir(pah) and callable(getattr(pah, scoreFunction)):
            scoreFunctionObj = eval('pah().' + scoreFunction)
        else:
            print "Score function not found!"
            sys.exit()
        if maxOptimalSolutions == -1:
            return self.traceback(sequences[0], sequences[1],
                                  self.computeMatrix(sequences[0], sequences[1], scoreFunctionObj), scoreFunctionObj)
        else:
            return self.traceback(sequences[0], sequences[1],
                                  self.computeMatrix(sequences[0], sequences[1], scoreFunctionObj), scoreFunctionObj,
                                  maxOptimalSolutions)

