
# coding: utf-8

from Bio import SeqIO
import pyspark
import numpy as np
import time
#from sompy.sompy import SOMFactory
from matplotlib import pyplot as plt
import gc
import pickle
import skfuzzy as fuzz
from scipy.sparse import csr_matrix


def k_mer_freq_count(read, k):
    uniqueValues, occurCount = np.unique(read, return_counts=True)
    count = np.zeros(pow(4,k))
    for i in range(len(uniqueValues)):
        count[uniqueValues[i]] = occurCount[i]
    return count

# Transforma os k-mers em índices na base 10
def k_mer(read, k):
    mer = list()
    for i in range(0, len(read)-k):
        index = 0
        for j, n in enumerate(read[i:i+k]):
            if n == -1:
                index = -1
                continue
            index = index + (n * pow(4, ((k-1)-j)))
        mer.append(index)
    mer = np.array(mer, dtype = np.int)
    return mer


# Converte A, C, G e T para 0, 1, 2 e 3 respectivamente
def convert(read):
    seq = read.seq
    rd = np.zeros(len(seq), dtype=np.int)
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


def kmer_count(filename, k, context):
    conf = pyspark.SparkConf().setAppName(context)
    conf = (conf.setMaster('local[*]')\
           .set('spark.executor.memory', '8G')\
           .set('spark.driver.memory', '8G')\
           .set('spark.driver.maxResultSize', '8G'))
    sc = pyspark.SparkContext(conf=conf)

    content = sc.parallelize(list(SeqIO.parse(filename, "fasta")))

    kmers = np.array(content.map(convert).map(lambda read: k_mer(read, k)).map(lambda read: k_mer_freq_count(read, k)).collect())
    sc.stop()
    return kmers
