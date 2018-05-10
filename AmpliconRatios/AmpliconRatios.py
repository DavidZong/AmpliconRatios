"""
Amplicon Ratios
"""

import os
import csv
from Bio import SeqIO
from Bio import pairwise2
import numpy as np


class AmpliconRatios(object):
    """
    Handles one paired end read from my protocol for NGS
    """

    def __init__(self, read_1, read_2):
        """
        The Amplicon ratio object
        :param read_1: the R1 fastq file
        :param read_2: the R2 fastq file
        """
        self.forward_read = read_1
        self.reverse_read = read_2
        self.forward_indices = ["ATACC", "TTCTT", "TCTGT",
                                "CATAA", "TCGAA", "CTTGC"]
        self.reverse_indices = ["CAATG", "CTGTG", "ACGCC",
                                "TATGG", "CCTAC", "TTAGC",
                                "TTCAC", "CGACC", "AATTG"]
        self.barcodes = ["TTCTCGGTCGGGTCATATCTAAGGT",
                         "GTGTGTCCGATGAACGCGACGTGAT",
                         "CTTGTTACGGGACCTAGTATCCCTA"]

    def parse_paired_read(self, output_folder, output_file_raw, output_file_processed):
        """
        The workhorse function that takes a pair of fastq files and returns
        location and barcode, saving as a .csv file (the raw one), also
        calculates the frequency of each well and the barcode in that well
        :return:
        """
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        destination = os.path.join(output_folder, output_file_raw)
        with open(destination, "w") as text_file:
            index = 0
            omitted = 0
            wr = csv.writer(text_file, quoting=csv.QUOTE_ALL)
            for record1, record2 in zip(
                    SeqIO.parse(self.forward_read, "fastq"),
                    SeqIO.parse(self.reverse_read, "fastq")):
                output = []
                indexing_result = self.match_indices(record1, record2)
                barcode_result = self.find_barcode(record1, record2)
                row = calculate_row(indexing_result)
                output.append(index)
                output = output + indexing_result + row + barcode_result
                print(output)
                if -1 in output:
                    omitted += 1
                else:
                    wr.writerow(output)
                index += 1
            print("Completed")
            print("Number of sequences:", index)
            print("Number of omissions:", omitted)
            accepted_percent = (index - omitted) / index
            print("Percent good:", accepted_percent)

    # # For debugging only, delete later
    # def parse_first_pair(self):
    #     """
    #     Gets the 1st record from a fastq file
    #     :param file: name of the fastq file
    #     :return: the record as a Biopython seq object
    #     """
    #     index = 0
    #     output = []
    #     records_f = SeqIO.parse(self.forward_read, "fastq")
    #     record_f = next(records_f)
    #     records_r = SeqIO.parse(self.reverse_read, "fastq")
    #     record_r = next(records_r)
    #     indexing_result = self.match_indices(record_f, record_r)
    #     barcode_result = self.find_barcode(record_f, record_r)
    #     output.append(index)
    #     output = output + indexing_result + barcode_result
    #     print(output)

    def match_indices(self, forward_sequence, reverse_sequence):
        """
        Takes in the forward and reverse sequence and tries to match it up to
        an index. Returns the index (as an integer), -1 if no match is found
        :param forward_sequence: the forward sequence as a Biopython Seqrecord
        :param reverse_sequence: the reverse sequence as a Biopython Seqrecord
        :return: an array where the first number is the forward's index, and
        the second number is the reverse's index, and the third is a string
        of the wellplate position (e.g. B2)
        """
        result = [None] * 2
        first5_f = forward_sequence[0:5].seq
        first5_r = reverse_sequence[0:5].seq
        result[0] = index_iterate(first5_f, self.forward_indices)
        result[1] = index_iterate(first5_r, self.reverse_indices)
        return result

    # # Ideally this should be used for finding barcodes to reduce rejected sequences
    # def align_barcode(self, forward_sequence, reverse_sequence, sensitivity):
    #     """
    #     Takes in the forward and reverse sequence and tries aligning the
    #     three different barcodes to both sequences.
    #     :param forward_sequence: the forward sequence as a Biopython Seqrecord
    #     :param reverse_sequence: the reverse sequence as a Biopython Seqrecord
    #     :param sensitivity: an integer that sets the minimum acceptable score
    #     :return: -1 if no alignments found, otherwise the index of the
    #     alignment along with the name of the barcode as a list of length 2
    #     """
    #     seq1 = forward_sequence.seq
    #     seq2 = reverse_sequence.seq
    #     seq2_rc = seq2.reverse_complement()
    #     for index, barcode in enumerate(self.barcodes):
    #         alignment_f = pairwise2.align.globalxx(seq1, barcode)
    #         alignment_r = pairwise2.align.globalxx(seq2_rc, barcode)

    def find_barcode(self, forward_sequence, reverse_sequence):
        """
        Takes in the forward and reverse sequence and tries finding the
        three different barcodes to both sequences.
        :param forward_sequence: the forward sequence as a Biopython Seqrecord
        :param reverse_sequence: the reverse sequence as a Biopython Seqrecord
        :return: -1 if no alignments found, otherwise the index of the
        alignment of the forward and reverse
        """
        output = [None] * 3
        seq1 = forward_sequence.seq
        seq2 = reverse_sequence.seq
        seq2_rc = seq2.reverse_complement()
        output[0] = index_iterate(seq1, self.barcodes)
        output[1] = index_iterate(seq2_rc, self.barcodes)
        if output[0] == output[1]:
            barcode_names = ['X', 'Y', 'Z']
            output[2] = barcode_names[output[0]]
        return output


def index_iterate(target, index_list):
    """
    A helper function to iterate through a list and find a target sequence
    Returns the index of the found sequence, or if not found -1
    :param target:
    :param index_list:
    :return: index if found, -1 if not found
    """
    result = -1
    for i, index in enumerate(index_list):
        if target.find(index) != -1:
            result = i
    return result

def calculate_row(result):
    """
    Given a pair of indicies, calculate which well it came from
    :return:
    """
    output = [None] * 2
    if (result[0] != -1) & (result[1] != -1):
        rows = ['B', 'C', 'D', 'E', 'F', 'G']
        cols = ['2', '3', '4', '5', '6', '7', '8', '9', '10']
        row = rows[result[0]]
        col = cols[result[1]]
        output[0] = (result[0] + 1) + (result[1] * 6)
        output[1] = row + col
    else:
        output[0] = -1
        output[1] = -1
    return output


def initialize_freq_table():
    """
    Sets up an empty frequency table
    :return: an empty frequency table
    """
    w, h = 3, 54;
    table = [[0 for x in range(w)] for y in range(h)]
    return table