
# coding: utf-8

import pyspark
import numpy as np
from scipy.sparse import lil_matrix

def k_mer_freq_count(read, k):
    uniqueValues, occurCount = np.unique(read, return_counts=True)
    return uniqueValues, occurCount

# Transforma os k-mers em índices na base 10
def k_mer(read, k):
    mer = list()
    for i in range(0, len(read)-k):
        index = 0
        for j, n in enumerate(read[i:i+k]):
            if n == -1 or n == -2:
                index = -1
                continue
            index = index + (n * pow(4, ((k-1)-j)))
        if index != -1:
            mer.append(index)
    mer = np.array(mer, dtype=np.int)
    return mer


# Converte A, C, G e T para 0, 1, 2 e 3 respectivamente
def convert(read):
    seq = read['seq']
    rd = np.zeros(len(seq))
    for i in range(len(seq)):
        if seq[i] == "\n":
            print("Erro!")
            continue
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
        rd[i] = -2
    return rd


def parse_fasta(filename):
    read = {'header': '', 'seq': ''}
    reads = list()
    lines =  open(filename, 'r').readlines()
    for i in range(len(lines)):
        if lines[i][0] == '>':
            read['header'] = lines[i][1:].split('\n')[0]
        else:
            if len(read['seq']) != 0:
                read['seq'] = read['seq'] + lines[i].split('\n')[0]
            else:
                read['seq'] = lines[i].split('\n')[0]
        if i < len(lines)-1:
            if lines[i+1][0] == '>':
                reads.append(read)
                read = {'header': '', 'seq': ''}
    return reads


def kmer_count(filename, k, context):
    conf = pyspark.SparkConf().setAppName(context)
    conf = (conf.setMaster('local[48]')\
           .set('spark.executor.memory', '250G')\
           .set('spark.driver.memory', '250G')\
           .set('spark.driver.maxResultSize', '250G'))
    sc = pyspark.SparkContext(conf=conf)

    content = sc.parallelize(parse_fasta(filename))

    kmers_list = content.map(convert).map(lambda read: k_mer(read, k)).map(lambda read: k_mer_freq_count(read, k)).collect()
    sc.stop()

    kmers = lil_matrix((len(kmers_list), 4**k))
    for i in range(0,len(kmers_list)):
        kmers[i, kmers_list[i][0]] = kmers_list[i][1]
#    print(kmers)

    return kmers
