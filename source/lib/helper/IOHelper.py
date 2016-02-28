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
import os
class IOHelper():
    """Helper class for reading an writing files in different formats."""
    def readFastaFile(self, inputFileName):
        """Reads a given fasta file and returns it as a array. 

            inputFileName:              The path (relative or absolut) to the input fasta file."""
        sequence = []
        if not os.path.exists(inputFileName):
            return sequence

        fileToRead = open(inputFileName, "r")
        i = 0
        for line in fileToRead.readlines():
            if line.startswith(">"):
                continue
            sequence.append(line.strip("\n"))
            i += 1
        fileToRead.close()
        return sequence

    def writeFastaFile(self, sequences, outputFileName):
        """Writes a the given sequences to a file in the fasta format.
            sequences:      All computed alignemnts.
                            A list of lists with two elements: [[,],...,[,]].
            outputFileName: The path (relative or absolut) and the output file name.
                            e.g.: "/path/to/file" or "file" to write it in the local directory."""
        if not outputFileName.endswith(".fas"):
            outputFileName += str(".fas")
        fileToWrite = open(outputFileName, "w")
        i = 0
        while i < len(sequences):
            for sequence in sequences[i]:
                fileToWrite.write('>Alignment '+ str(i) +' sequence ' + str(sequences[i].index(sequence)) + '\n')
                fileToWrite.write(sequence + '\n')
            i += 1
        fileToWrite.close()

    def writeGraphMLFile(self, clusteredNodesDictionary, outputFileName):
        """Writes a tree computed by the UpgmaWpgma class in graphML-format to specified outputFileName."""
        if not outputFileName.endswith(".graphml"):
            outputFileName += str(".graphml")
        fileToWrite = open(outputFileName, "w")
        fileToWrite.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
            +"\n<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\""
            +"\n\t\txmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
            +"\t\txsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n"
            +"\t\thttp://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">"
            +"\n\t<graph id=\"G\" edgedefault=\"undirected\">\n")
        for i in clusteredNodesDictionary:
            nodes = i.split(" ")
            fileToWrite.write("\t\t<node id=\"" + nodes[0] + "\"/>\n")
            fileToWrite.write("\t\t<node id=\"" + nodes[1] + "\"/>\n")
            fileToWrite.write("\t\t<node id=\"" + str(clusteredNodesDictionary[i]) + "\"/>\n")
        j = 0
        for i in clusteredNodesDictionary:
            nodes = i.split(" ")
            fileToWrite.write("\t\t<edge id=\"" + str(j) + "\" source=\"" + nodes[0] + "\" target=\""+ str(clusteredNodesDictionary[i]) + "\"/>\n")
            j += 1
            fileToWrite.write("\t\t<edge id=\"" + str(j) + "\" source=\"" + nodes[1] + "\" target=\""+ str(clusteredNodesDictionary[i]) + "\"/>\n")
            j += 1

        fileToWrite.write("\t</graph>\n</graphml>")
        fileToWrite.close()
    def writeRnaDotBracketNotation(self, sequence, pairedBases, outputFileName):
        """Writes a given RNA sequence and the computed matching bases in dot-bracket notation to the file outputFileName."""
        stack = {}
        for i in range (0, len(sequence)):
            if i in pairedBases:
                stack[i] = "("
                stack[pairedBases[i]] = ")"
            else:
                if not i in stack:
                    stack[i] = "."
        fileToWrite = open(outputFileName, "w")
        fileToWrite.write(sequence+"\n")
        for i in sorted(stack):
            fileToWrite.write(stack[i])
    def writeNewickTree(self, newickTree, outputFileName):
        fileToWrite = open(outputFileName, "w")
        fileToWrite.write(newickTree)
