"""
Amplicon Ratios
"""

import os
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

class AmpliconRatios(object):
    """
    Handles one paired end read
    """

    def __init__(self, read_1, read_2):
        """
        The Amplicon ratio object
        :param read_1: the R1 fastq file
        :param read_2: the R2 fastq file
        """
        self.forward_read = read_1
        self.reverse_read = read_2


    def parse_paired_read(self):
        """
        Takes 2 .fastq files and iterates though them as a pair
        Useful for pairing up forward and reverse reads as one
        :return: Nothing
        """
        index = 0
        print("Counting", self.forward_read, "and", self.reverse_read, "...")
        for record1, record2 in zip(
                SeqIO.parse(self.forward_read, "fastq"),
                SeqIO.parse(self.reverse_read, "fastq")):
            sequence1 = record1.seq
            sequence2 = record2.seq
            print("index %i Read 1: %s \n Read 2: %s" % (index, sequence1, sequence2))
            index += 1
        print("Number of sequences:", index)