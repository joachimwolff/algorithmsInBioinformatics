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
from helper import MultipleAlignmentHelper as mah


class NeedlemanWunschN3():
    """This class computes the Needleman-Wunsch algorithm with three sequences."""

    def __init__(self, sequence_a, sequence_b, sequence_c, score_function):
        """Initalize all variables and methods needed to compute the Needleman-Wunsch algorithm with three sequences.
            sequenceA:      A string with the first DNA sequence.
            sequenceB:      A string with the second DNA sequence.
            sequenceC:      A string with the third DNA sequence.
            scoreFunction:  The name of a weight function as a String which is defined 
                            in the pairwiseAlignmentHelper-class.
        """
        if score_function in dir(mah) and callable(getattr(mah, score_function)):
            score_function_obj = eval('mah().' + score_function)
        else:
            print "Score function not found!"
            sys.exit()

        self.computation_matrix = [[[]]]
        self.sequence_a = sequence_a
        self.sequence_b = sequence_b
        self.sequence_c = sequence_c
        self.score_function = score_function_obj
        self.i = 0
        self.j = 0
        self.traceback_stack = [[]]
        self.traceback_stack_index = 0
        self.indices_stack = [[]]
        self.computed_alignment = []

    def compute_matrix(self):
        """Computes the matrix which is needed by the Needleman-Wunsch algorithm for three sequences."""
        self.computation_matrix = [
            [[0 for i in range(len(self.sequence_c) + 1)] for j in range(len(self.sequence_b) + 1)] \
            for k in range(len(self.sequence_a) + 1)]
        # initalize matrix
        for i in range(1, len(self.sequence_a) + 1):
            self.computation_matrix[i][0][0] = self.computation_matrix[i - 1][0][0] \
                                               + self.score_function("", "", self.sequence_a[i - 1])
        for i in range(1, len(self.sequence_b) + 1):
            self.computation_matrix[0][i][0] = self.computation_matrix[0][i - 1][0] \
                                               + self.score_function("", "", self.sequence_b[i - 1])
        for i in range(1, len(self.sequence_c) + 1):
            self.computation_matrix[0][0][i] = self.computation_matrix[0][0][i - 1] \
                                               + self.score_function("", "", self.sequence_c[i - 1])
        for i in range(1, len(self.sequence_a) + 1):
            for j in range(1, len(self.sequence_b) + 1):
                self.computation_matrix[i][j][0] = self.computation_matrix[i - 1][j - 1][0] \
                                                   + self.score_function(self.sequence_a[i - 1], self.sequence_b[j - 1],
                                                                         "")
        for i in range(1, len(self.sequence_a) + 1):
            for k in range(1, len(self.sequence_c) + 1):
                self.computation_matrix[i][0][k] = self.computation_matrix[i - 1][0][k - 1] \
                                                   + self.score_function(self.sequence_a[i - 1], "",
                                                                         self.sequence_c[k - 1])
        for j in range(1, len(self.sequence_b) + 1):
            for k in range(1, len(self.sequence_c) + 1):
                self.computation_matrix[0][j][k] = self.computation_matrix[0][j - 1][k - 1] \
                                                   + self.score_function("", self.sequence_b[j - 1],
                                                                         self.sequence_c[k - 1])

        for i in range(1, len(self.sequence_a) + 1):
            for j in range(1, len(self.sequence_b) + 1):
                for k in range(1, len(self.sequence_c) + 1):
                    self.computation_matrix[i][j][k] = self.compute_minimum(i, j, k)

    def compute_minimum(self, i, j, k):
        """Compute the minimal value for a given cell of the matrix.
            The minimum is choosen of the following values:
                D(i-1, j-1, k-1) + w(a_i-1, b_j-1, c_k-1)
                D(i, j-1, k-1) + w(a_i, b_j-1, c_k-1)
                D(i-1, j, k-1) + w(a_i-1, b_j, c_k-1)
                D(i-1, j-1, k) + w(a_i-1, b_j-1, c_k)
                D(i, j, k-1) + w(a_i, b_j, c_k-1)
                D(i-1, j, k) + w(a_i-1, b_j, c_k)
                D(i, j-1, k) + w(a_i, b_j-1, c_k)
            i: index of sequence A
            j: index of sequence B
            k: index of sequence C
        """
        # no gap
        no_gap = self.computation_matrix[i - 1][j - 1][k - 1] \
                 + self.score_function(self.sequence_a[i - 1], self.sequence_b[j - 1], self.sequence_c[k - 1])
        # one gap
        gap_a = self.computation_matrix[i][j - 1][k - 1] \
                + self.score_function("", self.sequence_b[j - 1], self.sequence_c[k - 1])
        gap_b = self.computation_matrix[i - 1][j][k - 1] \
                + self.score_function(self.sequence_a[i - 1], "", self.sequence_c[k - 1])
        gap_c = self.computation_matrix[i - 1][j - 1][k] \
                + self.score_function(self.sequence_a[i - 1], self.sequence_b[j - 1], "")
        # two gaps
        gap_ab = self.computation_matrix[i][j][k - 1] + self.score_function("", "", self.sequence_c[k - 1])
        gap_bc = self.computation_matrix[i - 1][j][k] + self.score_function(self.sequence_a[i - 1], "", "")
        gap_ac = self.computation_matrix[i][j - 1][k] + self.score_function("", self.sequence_b[j - 1], "")
        possible_values = [no_gap, gap_a, gap_b, gap_c, gap_ab, gap_bc, gap_ac]
        return min(possible_values)

    def traceback(self, maximal_optimal_solutions=-1):
        """Computes the traceback for the Needleman-Wunsch n=3 matrix."""
        self.traceback_stack = [[]]
        self.indices_stack = [[len(self.computation_matrix) - 1, len(self.computation_matrix[0]) - 1,
                               len(self.computation_matrix[0][0]) - 1]]
        self.traceback_stack_index = 0
        traceback_done = False
        optimal_solutions_count = 0
        while not traceback_done:

            i = self.indices_stack[self.traceback_stack_index][0]
            j = self.indices_stack[self.traceback_stack_index][1]
            k = self.indices_stack[self.traceback_stack_index][2]
            optimal_solutions_count += 1
            split = False
            while i > 0 or j > 0 or k > 0:
                path_variable_i = i
                path_variable_j = j
                path_variable_k = k
                # no gap
                if i > 0 and j > 0 and k > 0:
                    if self.computation_matrix[i][j][k] == self.computation_matrix[i - 1][j - 1][k - 1] \
                            + self.score_function(self.sequence_a[i - 1], self.sequence_b[j - 1],
                                                  self.sequence_c[k - 1]):
                        self.traceback_stack[self.traceback_stack_index].append(mah.noGap)
                        path_variable_i -= 1  # change i
                        path_variable_j -= 1  # change j
                        path_variable_k -= 1  # change k
                        split = True

                # a gap in sequence a
                if j > 0 and k > 0:
                    if self.computation_matrix[i][j][k] == self.computation_matrix[i][j - 1][k - 1] \
                            + self.score_function("", self.sequence_b[j - 1], self.sequence_c[k - 1]):
                        if split == False:
                            self.traceback_stack[self.traceback_stack_index].append(mah.gapA)
                            path_variable_j -= 1
                            path_variable_k -= 1
                            split = True
                        else:
                            self.split([i, j - 1, k - 1], mah.gapA)
                # a gap in sequence b
                if i > 0 and k > 0:
                    if self.computation_matrix[i][j][k] == self.computation_matrix[i - 1][j][k - 1] \
                            + self.score_function(self.sequence_a[i - 1], "", self.sequence_c[k - 1]):
                        if split == False:
                            self.traceback_stack[self.traceback_stack_index].append(mah.gapB)
                            path_variable_i -= 1
                            path_variable_k -= 1
                        elif split == True:
                            self.split([i - 1, j, k - 1], mah.gapB)
                # a gap in sequence c
                if i > 0 and j > 0:
                    if self.computation_matrix[i][j][k] == self.computation_matrix[i - 1][j - 1][k] \
                            + self.score_function(self.sequence_a[i - 1], self.sequence_b[j - 1], ""):
                        if split == False:
                            self.traceback_stack[self.traceback_stack_index].append(mah.gapC)
                            path_variable_i -= 1
                            path_variable_j -= 1
                        elif split == True:
                            self.split([i - 1, j - 1, k], mah.gapC)
                # a gap in sequence a and b
                if k > 0:
                    if self.computation_matrix[i][j][k] == self.computation_matrix[i][j][k - 1] \
                            + self.score_function("", "", self.sequence_c[k - 1]):
                        if split == False:
                            self.traceback_stack[self.traceback_stack_index].append(mah.gapAB)
                            path_variable_k -= 1
                        elif split == True:
                            self.split([i, j, k - 1], mah.gapAB)
                # a gap in sequence a and c
                if j > 0:
                    if self.computation_matrix[i][j][k] == self.computation_matrix[i][j - 1][k] \
                            + self.score_function("", self.sequence_b[j - 1], ""):
                        if split == False:
                            self.traceback_stack[self.traceback_stack_index].append(mah.gapAC)
                            path_variable_j -= 1
                        elif split == True:
                            self.split([i, j - 1, k], mah.gapAC)
                # a gap in sequence b and c
                if i > 0:
                    if self.computation_matrix[i][j][k] == self.computation_matrix[i - 1][j][k] \
                            + self.score_function(self.sequence_a[i - 1], "", ""):
                        if split == False:
                            self.traceback_stack[self.traceback_stack_index].append(mah.gapBC)
                            path_variable_i -= 1
                        elif split == True:
                            self.split([i - 1, j, k], mah.gapBC)
                split = False
                i = path_variable_i
                j = path_variable_j
                k = path_variable_k
            if maximal_optimal_solutions != -1 and optimal_solutions_count >= maximal_optimal_solutions:
                break
            self.indices_stack[self.traceback_stack_index][0] = i
            self.indices_stack[self.traceback_stack_index][1] = j
            self.indices_stack[self.traceback_stack_index][2] = k
            l = 0
            all_tracebacks_computed = 0
            while l < len(self.indices_stack):
                if self.indices_stack[l][0] == 0 and self.indices_stack[l][1] == 0 and self.indices_stack[l][2] == 0:
                    all_tracebacks_computed += 1
                else:
                    self.traceback_stack_index = l
                    l = len(self.indices_stack)
                l += 1
            if all_tracebacks_computed >= len(self.indices_stack):
                traceback_done = True
                # all_tracebacks_computed = 0
        if maximal_optimal_solutions != -1 and optimal_solutions_count >= maximal_optimal_solutions:
            for i in range(0, maximal_optimal_solutions):
                self.computed_alignment.append(self.build_alignment(self.traceback_stack[i]))
        else:
            for i in range(0, len(self.traceback_stack)):
                self.computed_alignment.append(self.build_alignment(self.traceback_stack[i]))

    def split(self, index, gapSymbol):
        """Splits the actual traceback path into two paths.
            index:      The index values for the next cell of the path.
            gapSymbol:  A symbol for the computed step for the path."""
        self.traceback_stack.append(self.traceback_stack[self.traceback_stack_index][0:-1])
        self.traceback_stack[len(self.traceback_stack) - 1].append(gapSymbol)
        self.indices_stack.append(index)

    def build_alignment(self, tracebackStack):
        """Builds the alignment for one traceback path.
                tracebackStack: The computed tracebackpath as a list = []
            """
        i = 0
        j = 0
        k = 0
        l = len(tracebackStack) - 1
        alignment_of_a = ""
        alignment_of_b = ""
        alignment_of_c = ""

        while len(tracebackStack) > 0:
            try:
                tracebackElement = tracebackStack.pop(l)
                if mah.noGap == tracebackElement:
                    alignment_of_a += self.sequence_a[i]
                    alignment_of_b += self.sequence_b[j]
                    alignment_of_c += self.sequence_c[k]
                    i += 1
                    j += 1
                    k += 1
                elif mah.gapA == tracebackElement:
                    alignment_of_a += "-"
                    alignment_of_b += self.sequence_b[j]
                    alignment_of_c += self.sequence_c[k]
                    j += 1
                    k += 1
                elif mah.gapB == tracebackElement:
                    alignment_of_a += self.sequence_a[i]
                    alignment_of_b += "-"
                    alignment_of_c += self.sequence_c[k]
                    i += 1
                    k += 1
                elif mah.gapC == tracebackElement:
                    alignment_of_a += self.sequence_a[i]
                    alignment_of_b += self.sequence_b[j]
                    alignment_of_c += "-"
                    i += 1
                    j += 1
                elif mah.gapAB == tracebackElement:
                    alignment_of_a += "-"
                    alignment_of_b += "-"
                    alignment_of_c += self.sequence_c[k]
                    k += 1
                elif mah.gapAC == tracebackElement:
                    alignment_of_a += "-"
                    alignment_of_b += self.sequence_b[j]
                    alignment_of_c += "-"
                    j += 1
                elif mah.gapBC == tracebackElement:
                    alignment_of_a += self.sequence_a[i]
                    alignment_of_b += "-"
                    alignment_of_c += "-"
                    i += 1
                l -= 1
            except:
                print "An error occured."
                sys.exit()
        while i < len(self.sequence_a):
            alignment_of_a += self.sequence_a[i]
            i += 1
        while j < len(self.sequence_b):
            alignment_of_b += self.sequence_b[j]
            j += 1
        while k < len(self.sequence_c):
            alignment_of_b += self.sequence_c[k]
            k += 1
        alignment = [alignment_of_a, alignment_of_b, alignment_of_c]
        return alignment

    def execute(self, maximalOptimalSolutions=-1):
        """Method to start the computation of the Needleman-Wunsch algorithm with three sequences. It returns the computed alignment.
        [maximalOptimalSolutions]: Define how many optimal solutions should be computed. If not defined, all optimal solutions are computed."""
        self.compute_matrix()
        if maximalOptimalSolutions == -1:
            self.traceback()
        else:
            self.traceback(maximalOptimalSolutions)
        return self.computed_alignment
