from Bio import SeqIO
import os
import time

print("Working from: ", os.getcwd())

def count_fastq(filename):
    """
    Takes a .fastq and counts the number of sequences in the file
    Also measures the time it takes to do the counting
    :param filename: name of .fastq file to be counted
    :return: nothing
    """
    counts = 0
    print("Counting", filename, "...")
    start = time.time()
    for record in SeqIO.parse(filename, "fastq"):
        counts += 1
    end = time.time()
    elapsed = end - start
    print("Number of sequences in", filename, ":", counts)
    print("Operation took", elapsed, "seconds")

def list_fastq(filename):
    """
    Takes a .fastq and prints the sequences to console
    :param filename:
    :return:
    """
    print("The sequences from", filename)
    for index, record in enumerate(SeqIO.parse(filename, "fastq")):
        sequence = record.seq
        print("index %i Sequence: %s" % (index, sequence))

# count_fastq("testdata/DZ01-2_R1_001.fastq")
# count_fastq("testdata/DZ01-2_R2_001.fastq")
list_fastq("testdata/DZ01-2_R2_001.fastq")