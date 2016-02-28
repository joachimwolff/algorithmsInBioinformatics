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

import unittest
import os, sys
lib_path = os.path.abspath('../../')
sys.path.append(lib_path)
from helper import IOHelper as io

class IOHelperTestClass(unittest.TestCase):
    """Test class to check the correctness of the methods in IOHelper."""
    def test_readFastaFile(self):
        """Test method to test the correct reading of a fasta file."""
        if os.path.exists("testReadFasta.fas"):
            os.remove("testReadFasta.fas")

        # first test case: two sequences
        sequenceToWrite = [["ACGT", "ACGTAATTA"]]
        expectedSequence = ["ACGT", "ACGTAATTA"]
        io().writeFastaFile(sequenceToWrite, "testReadFasta.fas")
        readSequence = io().readFastaFile("testReadFasta.fas", multipleSequenceAlignment=False)
        self.assertEqual(expectedSequence, readSequence)

        # second test case: two sequences but there are multilpe ones
        sequenceToWrite = [["ACGT", "ACGTAATTA", "AGTTG"]]
        expectedSequence = ["ACGT", "ACGTAATTA", "AGTTG"]
        io().writeFastaFile(sequenceToWrite, "testReadFasta.fas")
        readSequence = io().readFastaFile("testReadFasta.fas", multipleSequenceAlignment=False)
        self.assertNotEqual(expectedSequence, readSequence)

        # third test case: multiple sequences
        readSequence = io().readFastaFile("testReadFasta.fas", multipleSequenceAlignment=True)
        self.assertEqual(expectedSequence, readSequence)

        os.remove("testReadFasta.fas")

    def test_writeFastaFile(self):
        """Test method to test the correct writing of a fasta file."""
        if os.path.exists("testWriteFasta.fas"):
            os.remove("testWriteFasta.fas")
        sequence = [["ACGT", "ACGTAATTA"]]
        expectedReadSequence = [">Alignment 0 sequence 0", "ACGT", ">Alignment 0 sequence 1", "ACGTAATTA"]
        readInputSequence = []

        # first test case, filename with extension
        io().writeFastaFile(sequence, "testWriteFasta.fas")     
        testInputFile = open("testWriteFasta.fas")
        for line in testInputFile.readlines():
            readInputSequence.append(line.strip("\n"))
        self.assertEqual(expectedReadSequence, readInputSequence)
        testInputFile.close()
        os.remove("testWriteFasta.fas")

        # second test case, filename without extension
        readInputSequence = []
        io().writeFastaFile(sequence, "testWriteFasta")
        testInputFile = open("testWriteFasta.fas")
        for line in testInputFile.readlines():
            readInputSequence.append(line.strip("\n"))
        self.assertEqual(expectedReadSequence, readInputSequence)
        testInputFile.close()
        os.remove("testWriteFasta.fas")

if __name__ == "__main__":
    unittest.main() # run all tests  