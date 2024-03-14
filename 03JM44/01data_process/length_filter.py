# TE annotation analysis
import numpy as np
import os

# while annotation result have different length, we can use this function to filter the sequence with length provided
def checkSeqLength(fileName, length):
    with open(fileName) as seq:
        for line in seq:
            if len(line) <= length:
                print(title + line)
            title = line

fileName = os.path.join(os.getcwd(), "../02result/teSeq.fasta")
checkSeqLength(fileName, 5)
