
# coding: utf-8

from Bio import SeqIO
from pyspark import SparkContext
import numpy as np
import time
import scipy.spatial as spt

def k_mer(read):
    pass

def dist(v1, v2):
    return spt.distance.cosine(v1,v2)

def calculate_distance(query, db):
    dim = [len(query), len(db)]
    D = np.empty(dim)
    for i in range(len(query)):
        for j in range(len(db)):
            D[i,j] = dist(query[i], db[j])
    return D

def convert(read):
    seq = read.seq
    rd = np.zeros(len(seq), dtype=np.float)
    for i in range(len(seq)):
        if seq[i] == 'A' or seq[i] == 'a':
            rd[i] = 0;
            continue
        if seq[i] == 'C' or seq[i] == 'c':
            rd[i] = 1;
            continue
        if seq[i] == 'G' or seq[i] == 'g':
            rd[i] = 2;
            continue
        if seq[i] == 'T' or seq[i] == 't':
            rd[i] = 3;
            continue
        if seq[i] == 'N' or seq[i] == 'n':
            rd[i] = -1;
            continue
    return rd

def load_sequences(file_content):

    reads = list();
    for read in file_content:
        reads.append(convert(read))

    return np.array(reads)

if __name__ == "__main__":
    gst = time.time()

    query_file = "SRX5702553_250seqs.fasta"
    db_file = "SRX5702553.fasta"

    #sc = SparkContext("local[4]", "genenmatch")
    #sc.setSystemProperty('spark.executor.memory', '8g')
    #sc.setSystemProperty('spark.driver.memory', '8g')
    #reads = sc.parallelize(list(SeqIO.parse(filename, "fasta")))
    #conv = reads.map(convert).collect()

    pst = time.time()
    query_content = list(SeqIO.parse(query_file, "fasta"))
    db_content = list(SeqIO.parse(db_file, "fasta"))
    pend = time.time()
    print("Read time: ", pend - pst)

    pst = time.time()
    db_reads = load_sequences(db_content)
    query_reads = load_sequences(query_content)
    del(db_content)
    del(query_content)
    pend = time.time()
    print("Convertion time: ", pend - pst)

    pst = time.time()
    D = calculate_distance(query_reads, db_reads)
    print(D)
    pend = time.time()
    print("Distance time: ", pend - pst)

    gend = time.time()
    print("Walltime: ", gend - gst)
