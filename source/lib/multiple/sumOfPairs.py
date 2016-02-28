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
# Sum of pairs algorithm
from helper import PairwiseAlignmentHelper as pah
import sys


class SumOfPairs():
    """This class computes the Sum-of-pairs algorithm by Carillo and Lipman:
            Carrillo, Humberto, and David Lipman.
            "The multiple sequence alignment problem in biology."
            SIAM Journal on Applied Mathematics 48.5 (1988): 1073-1082.
            http://www.academia.edu/download/30855770/Articulo03.pdf"""
    def __init__(self, sequences, similarity_score):
        """To initialize a object of the SumOfPairs class please define a list with the multiple sequence alignment and
         a similarity score method which is defined in class PairwiseAlignmentHelper of package helper.
         sequences: The multiple alginment as a list.
         similarity_score: The scoring functions name as a string."""
        self.sequences = sequences
        if similarity_score in dir(pah) and callable(getattr(pah, similarity_score)):
            similarity_score_obj = eval('pah().' + similarity_score)
        else:
            print "Score function not found!"
            sys.exit()
        self.score_function = similarity_score_obj

    def execute(self):
        """Run this method to compute the sum of pairs scoring for multiple alignment."""
        score_value = 0
        for i in range(0, len(self.sequences)):
            for j in range(i+1, len(self.sequences)):
                score_value += self.score(self.sequences[i], self.sequences[j])
        return score_value

    def score(self, sequence_a, sequence_b):
        """Returns the pairwise alignment for sequence_a and sequence_b."""
        score_value = 0
        for i in range(0, max(len(sequence_a), len(sequence_b))):
            if i < len(sequence_a) and i < len(sequence_b):
                score_value += self.score_function(sequence_a[i], sequence_b[i])
            elif i < len(sequence_a):
                score_value += self.score_function(sequence_a[i], "")
            elif i < len(sequence_b):
                score_value += self.score_function("", sequence_b[i])
            i += 1
        return score_value
