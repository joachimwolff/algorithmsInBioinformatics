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
import sys
from pairwise import NeedlemanWunsch
from math import log10
from helper import PairwiseAlignmentHelper as pah
from multiple import UpgmaWpgma
import random


class FengDoolittle():
    """This class computes the Feng-Doolittle algorithm by Da-Fei Feng and Russell F. Doolittle:
            Feng, Da-Fei, and Russell F. Doolittle.
            "Progressive sequence alignment as a prerequisitetto correct phylogenetic trees."
            Journal of molecular evolution 25.4 (1987): 351-360.
            http://dna.bio.puc.cl/cardex/papersbio252/Grupo06-2013.pdf"""
    def __init__(self, sequences, weightFunction, similarityScore):
        """To initialize an object of class FengDoolittle you have to define:
                sequences:          A list of sequences for the multiple alignment.
                weightFunction:     A string containing the name of your preferred weight function.
                                    The weight function have to be defined in package helper, class PairwiseAlignmentHelper.
                similarityScore:    A string containing the name of your preferred similarity score like pam250.
                                    The similarity score have to be defined in package helper, class PairwiseAlignmentHelper."""
        self.sequences = sequences
        self.alignments = []
        self.weightFunction = weightFunction
        if similarityScore in dir(pah) and callable(getattr(pah, similarityScore)):
            similarityScoreObj = eval('pah().' + similarityScore)
        else:
            print "Score function not found!"
            sys.exit()
        self.similarityScore = similarityScoreObj
        self.alignmentToIndexMapping = {}
        self.sequenceToIndexMapping = {}
        self.distanceDictionary = {}
        self.newickTree = ""
        self.orderToAlign = []

    def computeAlignments(self):
        """This function computes all pairwise alignments between every sequence with the Needleman-Wunsch algorithm."""
        nw = NeedlemanWunsch()
        alignmentsAppend = self.alignments.append
        for i in range(0, len(self.sequences)):
            for j in range(i + 1, len(self.sequences)):
                alignmentsAppend([nw.compute([self.sequences[i], self.sequences[j]], self.weightFunction,1)[0], i, j])

    def computeDistanceDictionary(self):
        """This function computes the distance between every alignment. The distances are used to generate a phylogenetic tree."""
        for i in range(0, len(self.alignments)):
            index = str(self.alignments[i][1]) + " " + str(self.alignments[i][2])
            self.distanceDictionary[index] = self.similarityToDistance(self.alignments[i][0])

    def similarityToDistance(self, alignment):
        """Computes from the given similarity the distance measure."""
        sMax = self.similarity(alignment[0], alignment[0]) + self.similarity(alignment[1], alignment[1])
        sMax /= 2
        alignmentAsList = list(alignment[0])
        alignmentAsList1 = list(alignment[1])
        random.shuffle(alignmentAsList)
        random.shuffle(alignmentAsList1)
        alignmentShuffel0 = "".join(alignmentAsList)
        alignmentShuffel1 = "".join(alignmentAsList1)

        sRand = self.similarity(alignmentShuffel0, alignmentShuffel1)
        if sMax == sRand:
            sRand = sRand - 0.0001
        else:
            sEff = (self.similarity(alignment[0], alignment[1]) - sRand) / float(sMax - sRand)
        if sEff <= 0.0:
            return 1
        distance = -log10(sEff)
        return distance

    def similarity(self, a, b):
        """Returns the similarity of two sequences a and b with the similarity score defined at the initialization."""
        similarity = 0
        for i in range(0, len(a)):
            similarity += self.similarityScore(a[i], b[i])
        return similarity

    def buildTree(self):
        """This function computes the phylogenetic tree with UPGMA and stores it in the Newick-Tree format."""
        upgma = UpgmaWpgma(self.distanceDictionary, len(self.sequences))
        upgma.compute_clustering()
        self.newickTree = upgma.get_newick_tree()

    def buildMultipleAlignment(self, group0, group1):
        """This function returns which is the best pairwise alignment out of all alignments of group0 and group1."""
        highestScore = 0
        optimalAlignment = []
        for i in group0:
            for j in group1:
                nw = NeedlemanWunsch()
                alignment = nw.compute([i[0], j[0]], self.weightFunction, 1)
                score = self.similarity(alignment[0][0], alignment[0][1])
                if highestScore < score:
                    highestScore = score
                    optimalAlignment = [alignment[0][0], alignment[0][1], i[1], j[1]]
        return optimalAlignment


    def computeOrderOfSequencesToAlign(self):
        """This function computes out of the phylogenetic tree in which order the sequences are aligned."""
        indexBegin = 0
        indexEnd = len(self.newickTree)
        while indexEnd != -1:
            indexBegin = self.newickTree.rfind("(", indexBegin, indexEnd)
            if indexBegin == -1:
                break
            i = indexBegin + 1
            stack = 0
            while stack >= 0 and i < len(self.newickTree):
                if self.newickTree[i] == "(":
                    stack += 1
                elif self.newickTree[i] == ")":
                    stack -= 1
                i += 1
            indexEnd = i

            group0 = ""
            group1 = ""
            substring = self.newickTree[indexBegin:indexEnd]
            if substring[1] != "(":
                indexGroup0 = substring.find(",")
                group0 = substring[0:indexGroup0].strip(",")
                group1 = substring[indexGroup0:-1].strip(",")
            else:
                k = 1
                stack = 0
                while k < len(substring):
                    if substring[k] == "(":
                        stack += 1
                    elif substring[k] == ")":
                        stack -= 1
                    k += 1
                    if stack <= 0:
                        break
                group0 = substring[0:k].strip(",")
                group1 = substring[k:-1].strip(",")
            group0List = group0.split(",")
            group1List = group1.split(",")
            list0 = []
            list1 = []
            for j in group0List:
                list0.append(int(j.strip("(").strip(")").strip(",")))
            for j in group1List:
                list1.append(int(j.strip("(").strip(")").strip(",")))

            self.orderToAlign.append(sorted([sorted(list0), sorted(list1)]))
            indexEnd = indexBegin
            indexBegin = 0

    def computeMultipleAlignment(self):
        """This function returns the multiple sequence alignment."""
        self.computeAlignments()
        self.computeDistanceDictionary()
        self.buildTree()
        self.computeOrderOfSequencesToAlign()
        i = 0
        indexAlignments = {}
        # create index to algnment realation
        while i < len(self.orderToAlign):
            if len(self.orderToAlign[i][0]) == 1 and len(self.orderToAlign[i][1]):
                for j in self.alignments:
                    if (j[1] == self.orderToAlign[i][0][0] and j[2] == self.orderToAlign[i][1][0]):
                        indexAlignments[self.orderToAlign[i][0][0]] = j[0][0]
                        indexAlignments[self.orderToAlign[i][1][0]] = j[0][1]
                        break
                    elif(j[1] == self.orderToAlign[i][1][0] and j[2] == self.orderToAlign[i][0][0]):
                        indexAlignments[self.orderToAlign[i][0][0]] = j[0][1]
                        indexAlignments[self.orderToAlign[i][1][0]] = j[0][0]
                        break
            elif len(self.orderToAlign[i][0]) == 1:
                indexAlignments[self.orderToAlign[i][0][0]] = self.sequences[self.orderToAlign[i][0][0]]
            elif len(self.orderToAlign[i][1]) == 1:
                try:
                    indexAlignments[self.orderToAlign[i][1][0]] = self.sequences[self.orderToAlign[i][1][0]]
                except:
                    print "Exception!"
                    print "i: ", i
                    print "OrderToAlign: ", self.orderToAlign
                    print "orderAlign:", self.orderToAlign[i][1][0]
                    print self.sequences
            i += 1

        for i in self.orderToAlign:
            # one sequence with one sequence
            if len(i[0]) == 1 and len(i[1]):
                indexAlignments[i[0][0]] = indexAlignments[i[0][0]].replace("-", "X")
                indexAlignments[i[1][0]] = indexAlignments[i[1][0]].replace("-", "X")
            # one sequence with one group
            # two groups
            else:
                group0 = []
                group1 = []
                for j in i[0]:
                    group0.append([indexAlignments[j], j])
                for j in i[1]:
                    group1.append([indexAlignments[j],j])
                pairwiseAlignment = self.buildMultipleAlignment(group0, group1)
                indexAlignments[pairwiseAlignment[2]] = pairwiseAlignment[0].replace("-", "X")
                indexAlignments[pairwiseAlignment[3]] = pairwiseAlignment[1].replace("-", "X")

                for j in i[0]:
                    nw = NeedlemanWunsch()
                    alignment = nw.compute([pairwiseAlignment[0], indexAlignments[j]], self.weightFunction, 1)
                    indexAlignments[j] = alignment[0][1]
                for j in i[1]:
                    nw = NeedlemanWunsch()
                    alignment = nw.compute([pairwiseAlignment[1], indexAlignments[j]], self.weightFunction, 1)
                    indexAlignments[j] = alignment[0][1]
                for j in indexAlignments:
                    indexAlignments[j] = indexAlignments[j].replace("-", "X")
        return indexAlignments
