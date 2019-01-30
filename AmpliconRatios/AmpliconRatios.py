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

    def parse_paired_read(self, output_folder, output_file_raw,
                          output_file_processed):
        """
        The workhorse function that takes a pair of fastq files and returns
        location on plate and barcode, saving as a .csv file (the raw one),
        calculates the frequency of each well and the barcode in that well
        :return:
        """
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        destination = os.path.join(output_folder, output_file_raw)
        with open(destination, "w") as text_file:
            index = 0
            omitted = 0
            wr = csv.writer(text_file, quoting=csv.QUOTE_ALL , lineterminator='\n')
            freq_table = initialize_freq_table()
            for record1, record2 in zip(
                    SeqIO.parse(self.forward_read, "fastq"),
                    SeqIO.parse(self.reverse_read, "fastq")):
                output = []
                indexing_result = self.match_indices(record1, record2)
                barcode_result = self.find_barcode(record1, record2)
                if -1 in barcode_result:
                    barcode_result = self.align_barcode(record1, record2, 20)
                row = calculate_row(indexing_result)
                output.append(index)
                output = output + indexing_result + row + barcode_result
                print(output)
                wr.writerow(output)
                if -1 in output:
                    omitted += 1
                else:
                    row = output[3]
                    actual_barcode = output[-2]
                    freq_table[row-1][actual_barcode] += 1
                index += 1
        print("Completed")
        print("Number of sequences:", index)
        print("Number of omissions:", omitted)
        good_sequences = index - omitted
        accepted_percent = good_sequences / index
        print("Percent good:", accepted_percent)
        print("Saving frequencies table...")
        # for j in range(54):
        #     sum_of_well = freq_table[j][0] + freq_table[j][1] + freq_table[j][2]
        #     freq_table[j][3:6] = [quantity / sum_of_well for quantity in freq_table[j][0:3]]
        #     freq_table[6] = sum_of_well / good_sequences
        destination = os.path.join(output_folder, output_file_processed)
        with open(destination, "w") as text_file:
            wr = csv.writer(text_file, quoting=csv.QUOTE_ALL, lineterminator='\n')
            wr.writerows(freq_table)
        print("Done.")

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

    # Ideally this should be used for finding barcodes to reduce rejected sequences
    def align_barcode(self, forward_sequence, reverse_sequence, sensitivity):
        """
        Takes in the forward and reverse sequence and tries aligning the
        three different barcodes to both sequences.
        :param forward_sequence: the forward sequence as a Biopython Seqrecord
        :param reverse_sequence: the reverse sequence as a Biopython Seqrecord
        :param sensitivity: an integer that sets the minimum acceptable score, perfect is 25
        :return: -1 if no alignments found, otherwise the index of the
        alignment along with the name of the barcode as a list of length 2
        """
        score_f = [None] * 3
        score_r = [None] * 3
        seq1 = forward_sequence.seq
        seq2 = reverse_sequence.seq
        seq2_rc = seq2.reverse_complement()
        for index, barcode in enumerate(self.barcodes):
            score_f[index] = pairwise2.align.localms(seq1, barcode, 1, -1, -5, -1, score_only=True)
            score_r[index] = pairwise2.align.localms(seq2_rc, barcode, 1, -1, -5, -1, score_only=True)
        max_f = max(score_f)
        max_r = max(score_r)
        if max(max_f, max_r) > sensitivity:
            output = [None] * 3
            output[0] = score_f.index(max_f)
            output[1] = score_r.index(max_r)
            barcode_names = ['X', 'Y', 'Z']
            output[2] = barcode_names[output[0]]
            return output
        else:
            return [-1, -1, None]

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
