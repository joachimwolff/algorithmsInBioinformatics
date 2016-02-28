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
from helper import PairwiseAlignmentHelper as pah
class MultipleAlignmentHelper():
  
    noGap = 0
    gapA = 1
    gapB = 2
    gapC = 3
    gapAB = 4
    gapBC = 5
    gapAC = 6
    
    def weightFunctionDifference(self, a, b, c):
        """Weight function with 0 if a==b==c, 1 if a==b, a==c or b==c, 2 else."""
        if a == b and b == c:
            return 0
        elif a == b:
            return 1
        elif b == c:
            return 1
        elif a ==c :
            return 1
        else:
            return 2
    def createDataForUpgmaWpgma(self, sequences):
        """Preprocessing of the sequences for the upgm/wpgm algorithm."""
        differenceDictionary = {}
        sequenceToIdMapping = {}
        sequenceToLengthMapping = {}
        mappingCount = 0
        for i in sequences:
            sequenceToIdMapping[i] = mappingCount
            sequenceToLengthMapping[mappingCount] = len(i)
            mappingCount += 1

        differenceScore = 0
        for i in range(0, len(sequences)):
            for j in range(i+1, len(sequences)):
                for k in range(0, max(len(sequences[i]), len(sequences[j]))):
                    if k < len(sequences[j]) and k < len(sequences[i]):
                        differenceScore += pah().weightFunctionDifference(sequences[i][k], sequences[j][k])
                    elif k < len(sequences[i]):
                        differenceScore += pah().weightFunctionDifference(sequences[i][k], "-", )
                    elif k < len(sequences[j]):
                        differenceScore += pah().weightFunctionDifference("-", sequences[j][k])
                key = str(i) + " " + str(j)
                differenceDictionary[key] = differenceScore
                differenceScore = 0
        return [differenceDictionary, sequenceToIdMapping, sequenceToLengthMapping]                       





