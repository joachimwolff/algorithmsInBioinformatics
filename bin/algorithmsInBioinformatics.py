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
# main class

import argparse
import os, sys
if os.name == "posix":
    lib_path = os.path.abspath('../lib')
elif os.name == "nt":
    lib_path = os.path.abspath('..\lib')
sys.path.append(lib_path)

from helper import IOHelper as io
from helper import MultipleAlignmentHelper as mah
from pairwise import NeedlemanWunsch as nw
from pairwise import Gotoh
from multiple import NeedlemanWunschN3 as NW3
from multiple import UpgmaWpgma
from multiple import FengDoolittle
from multiple import SumOfPairs
from structurePrediction import Nussinov
def main():
    """Method to parse the arguments and start the defined algorithms."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--algorithm",
                        choices=["nw", "gotoh", "nw3", "fengDoolittle", "sumOfPairs","upgma", "wpgma", "nussinov"],
                        required=True,
        help="Define which algorithm should be executed. "
             "\nOptions are: 'nw' for the algorithm of Needleman and Wunsch,\n"
             "'gotoh' for the algorithm of Osamu Gotoh, \n"
             "'nw3' for the Needleman-Wunsch algorithm with three sequences, \n"
            "'fengDoolittle' for the heuristic multiple sequence alignment algorithm by Da-Fei Feng and Russell F. Doolittle,"
            "'sumOfPairs' for the scoring of a multiple sequence alignment by Humberto Carrillo and David Lipman."
             "'upgma' or 'wpgma' is a clustering method to generate phylogenetic trees, \n"
             "'nussinov' for the RNA secondary structure prediction algorithm by Ruth Nussinov.")
    parser.add_argument("-f", "--inputFile", dest="inputFile",
                        help="Define the file in which the input sequences are defined. It have to be in fasta-format.")
    parser.add_argument("-o", "--outputFile", help="Define in which file the output should be written. "
                                        "If not defined, it is written to \"outputFile.fas\" in the local directory.")
    parser.add_argument("-gc", "--gapCost", dest="gapCost",
                        help="Name of a gap function definde in class PairwiseAligmentHelper.")
    parser.add_argument("--numberOfSolutions", dest="numberOfSolutions",
                        help="Define the number of optimal solutions the Needleman-Wunsch algorithm should compute.")
    parser.add_argument("--outputFormat", dest="outputFormat", choices=["graphML", "newickTree"],
                        help="Define the output format of the output file. "
                             "This function is only parsed if you choose 'upgma' or 'wpgma' as an algorithm. Default is"
                             " Newick tree")
    parser.add_argument("--scoring", dest="similarityScore",
                        help="Name of a similarity score defined in class PairwiseAligmentHelper. Per default "
                             "\"pam\" and \"blosum\" (pam250 and blosum62) are implemented. Feel free to extend, you can find the "
                             "file \"PairwiseAligmentHelper.py\" in lib/helper. If this option is not defined, the pam250 matrix is choosen.")
    parser.add_argument("--gapPenalty", dest="gapPenalty", help="Define a gap penalty. Default for pam is 8 and blosum 6.")
    args = parser.parse_args()

    outputFile = ""
    weightFunction = ""
    if args.outputFile:
        outputFile = args.outputFile
    if args.similarityScore:
        weightFunction = args.similarityScore

    sequences = getSequencesFromFile(args.inputFile)
    if len(sequences) > 1:

        # pairwise alignment
        if args.algorithm == "nw":
            if outputFile == "":
                outputFile = "needlemanWunsch.fas"
            if weightFunction == "":
                weightFunction = "pam"
            numberOfSolutions = -1
            if args.numberOfSolutions:
                numberOfSolutions = args.numberOfSolutions
            needlemanWunsch(sequences[0:2], scoreFunction = weightFunction, outputFile = outputFile, numberOfSolutions=numberOfSolutions)

        elif args.algorithm == "gotoh":
            if outputFile == "":
                outputFile = "gotoh.fas"
            if weightFunction == "":
                weightFunction = "pam"
            gapCost = "gapCost"
            if args.gapCost:
                gapCost = args.gapCost
            gotoh(sequences[0:2], scoreFunction = weightFunction, costFunction =  gapCost, outputFile = outputFile)

        # multiple alignment

        elif args.algorithm == "upgma" or args.algorithm == "wpgma":
            newickTree = True
            if args.outputFormat == "graphML":
                newickTree = False
            if outputFile == "":
                if args.algorithm == "upgma":
                    outputFile = "upgma"
                else:
                    outputFile = "wpgma"
            upgmaWpgma(args.algorithm == "upgma", sequences, outputFile, newickTree)

        elif args.algorithm == "fengDoolittle":
            if outputFile == "":
                outputFile = "fengDoolittle.fas"
            if weightFunction == "":
                weightFunction = "pam"
            similarityScore = "pam"
            if args.similarityScore:
                similarityScore = args.similarityScore
            fengDoolittle(sequences, weightFunction, similarityScore, outputFile)
        elif args.algorithm == "sumOfPairs":
            similarityScore = "pam"
            if args.similarityScore:
                similarityScore = args.similarityScore
            if args.gapPenalty:
                sumOfPairs(sequences, similarityScore, args.gapPenalty)
            else:
                sumOfPairs(sequences, similarityScore)

        elif args.algorithm == "nw3":
            if not (len(sequences) == 3):
                print "Wrong number of input sequences. Needleman-Wunsch n=3 needs exactly three sequences; ", \
                    len(sequences) , " sequences are given."
                sys.exit()
            if weightFunction == "":
                weightFunction = "pam"
            if outputFile == "":
               outputFile = "nw3.fas"
            needlemanWunschN3(sequences[0:3], weightFunction = weightFunction, outputFile = outputFile)

    # multiple alignment

    elif len(sequences) == 1:
        # structure prediction
        if args.algorithm == "nussinov":
            if outputFile == "":
                outputFile = "nussinov.dotBracket"
            nussinov(sequences[0:1], outputFile)
        else:
            print "You have defined only one input sequence, but your defined algorithm \'",\
                args.algorithm, "\' needs at least two sequences."
    else:
        print "No sequences in input file defined."
        sys.exit(0)
def getSequencesFromFile(inputFile):
    """Parse the input file to get the sequences. Returns the sequences as an array.
            inputFile:  A fasta format file with the input sequences."""
    sequences = io().readFastaFile(inputFile)
    return sequences
def needlemanWunsch(sequences, scoreFunction, outputFile, numberOfSolutions):
    """Executes the Needleman-Wunsch algorithm with a default score function defined as: a == b -> 0 and a !=b --> 1.\n
    Stores the alignments per default in file needlemanWunsch.fas.
    To change the score function define a function in class PairwiseAligmentHelper and define the name as an input paramter.
        scoreFunction:      The name of the weigh function which is defined in class PairwiseAligmentHelper.
        outputFile:         The path to the output file.
        numberOfSolutions:  Maximal number of optimal solutions which should be computed."""
    print "\nThe following sequences are given:"
    for i in sequences:
        print i
    print "\nComputing solution...\n\n"
    result = nw().compute(sequences, scoreFunction, int(numberOfSolutions), scoringValue=True)
    print "\nScore: ", result[1]
    print "Number of optimal solutions: ", len(result[0])
    print "\nOne solution is:\n", result[0][0][0], "\n", result[0][0][1]
    print "\nFor more solutions look in the file \"needlemanWunsch.fas\" in the bin directory.\n"
    io().writeFastaFile(result[0], outputFile)
def gotoh(sequences, scoreFunction="weightFunctionDifference", costFunction="gapCost", outputFile="gotoh.fas"):
    """Executes the Gotoh algorithm with a default score function defined as: a == b -> 0 and a !=b --> 1 and a cost function defined as: g(x) = 2 + k.\n
    Stores the alignments per default in file gotoh.fas.
    To change the score or cost function define a function in class PairwiseAligmentHelper and define the name as an input paramter.
        scoreFunction:  The name of the weigh function which is defined in class PairwiseAligmentHelper.
        costFunction:   The name of the gap cost function which is defined in class PairwiseAligmentHelper.
        outputFile:     The path to the output file.
        """
    print "The following sequences are given:"
    for i in sequences:
        print i
    print "Computing solution..."
    gotoh = Gotoh(sequences[0], sequences[1], scoreFunction, costFunction)
    result = gotoh.compute()
    io().writeFastaFile(result, outputFile)
    print "Number of solutions: ", len(result)
    print "Score:", max(gotoh.computationMatrix[0][-1][-1], max(gotoh.computationMatrix[1][-1][-1], gotoh.computationMatrix[2][-1][-1]))
    print "One solution is:\n", result[0][0], "\n", result[0][1]
    print "For more solutions look in the file \"gotoh.fas\" in the bin directory."

def needlemanWunschN3(sequences, weightFunction="weightFunctionDifference", outputFile="nw3.fas"):
    """Executes the Needleman-Wunsch algorithm with three sequences"""
    print "\nThe following sequences are given:"
    for i in sequences:
        print i
    print "\nComputing solution...\n\n"
    nw3 = NW3(sequences[0], sequences[1], sequences[2], weightFunction)
    result = nw3.execute()

    io().writeFastaFile(result, outputFile)
    print "\nScore: ", nw3.computation_matrix[-1][-1][-1]
    print "Number of optimal solutions: ", len(result)
    print "\nOne solution is:\n", result[0][0], "\n", result[0][1], "\n", result[0][2]
    print "\nFor more solutions look in the file \"nw3.fas\" in the bin directory.\n"

def upgmaWpgma(upgmaWpgma, sequences, outputFile, fileFormat):
    """Executes the a phylogenetic clustering with a upgm or wpgm weighting.
        sequences:  All defined input sequences as a list.
        outputFile: The name of the output file
        fileFormat: The file format of the output file"""
    #create
    print "The following sequences are given:"
    for i in sequences:
        print i
    print "Computing clustering..."
    data = mah().createDataForUpgmaWpgma(sequences)
    if upgmaWpgma:
        upgma = UpgmaWpgma(data[0], len(data[1]))
        upgma.compute_clustering()
        if not fileFormat:
            outputFile += ".graphML"
            io().writeGraphMLFile(upgma.mapping, outputFile)
            print "Clustering written as graphML file: ", os.path.abspath(outputFile)
        else:
            outputFile += ".newickTree"
            cluster = upgma.get_newick_tree(with_edge_weights=True)
            io().writeNewickTree(cluster, outputFile)
            print "Computed upgma cluster: ", cluster
            print "The clustering was also written to: ", os.path.abspath(outputFile)
    else:
        wpgma = UpgmaWpgma(data[0], len(data[1]), False, data[2])
        wpgma.compute_clustering()
        if not fileFormat:
            outputFile += ".graphML"
            io().writeGraphMLFile(wpgma.mapping, outputFile)
            print "Clustering written as graphML file: ", os.path.abspath(outputFile)
        else:
            outputFile += ".newickTree"
            cluster = wpgma.get_newick_tree(with_edge_weights=True)
            io().writeNewickTree(cluster, outputFile)
            print "Computed wpgma cluster: ", cluster
            print "The clustering was also written to: ", os.path.abspath(outputFile)


def nussinov(sequence, outputFile):
    """Executes the RNA-folding algorithm from Nussinov.
        sequence:   The RNA-sequnce as a list.
        outputFile: The name of the output file."""
    print "\nThe following sequence is given:"
    print sequence[0]
    print "\n"
    nussinov = Nussinov(sequence[0])
    nussinov.execute()
    print "\nDot-bracket: "
    io().writeRnaDotBracketNotation(sequence[0], nussinov.pairedBases, outputFile)
    print "The result was also written to: ", os.path.abspath(outputFile)

def sumOfPairs(sequences, scoringFunction, gapPenalty=-1):
    """This method scores a multiple sequence alignment with the sum of pairs algorithm.
        sequences:          The multiple sequence alignment.
        scoringFunction:    Name of a similarity score defined in class PairwiseAligmentHelper."""
    print "The following sequences are given:"
    for i in sequences:
        print i
    if gapPenalty == -1:
        sof = SumOfPairs(sequences, scoringFunction)
    else:
        sof = SumOfPairs(sequences, scoringFunction, gapPenalty)
    print "Sum-of-pairs scoring: ", sof.execute()
def fengDoolittle(sequences, weightFunction, similarityScore, outputFile):
    """Executes the heuristic multiple sequence alignment by Feng and Doolittle.
        sequences:          All input sequnces to align.
        weightFunction:     The weight function defined in class PairwiseAlignmentHelper for the Needleman-Wunsch algorithm to compute the optimal local alignment.
        similarityScore:    Name of a similarity score defined in class PairwiseAligmentHelper.
        outputFile:         The output file name."""
    fd = FengDoolittle(sequences, weightFunction, similarityScore)
    alignmentDict = fd.computeMultipleAlignment()
    alignment = [[]]
    for i in alignmentDict:
        alignment[0].append(alignmentDict[i])
    io().writeFastaFile(alignment, outputFile)
    print "Input sequences:\n"
    for i in sequences:
        print i
    print "\nAlignment:"
    for i in alignmentDict:
        print alignmentDict[i]
    print sumOfPairs(alignment[0], weightFunction)


if __name__ == "__main__":
    # try:
    main()
    # except:
    #     "You discovered a bug! Please write an email to wolffj@informatik.uni-freiburg.de with your input parameters and I try to fix it."