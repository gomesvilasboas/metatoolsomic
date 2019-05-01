
# coding: utf-8

from Bio import SeqIO
import pyspark
import numpy as np
import time
import scipy.spatial as spt
import gc


def prepare_data(kmers):
    dim = [kmers.shape[0], len(max(kmers, key=len))]
    print(dim)
    data = np.zeros(dim, dtype=np.int)
    for i in range(dim[0]):
        data[i][0:len(kmers[i])] = kmers[i][:]
    return data

# Transforma os k-mers em índices na base 10
def k_mer(read):
    k = 31
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


if __name__ == "__main__":
    gst = time.time()

    filename = "/home/fabricio/ncbi/public/sra/SRR8945190.fasta"
    #filename = "/home/fabricio/ncbi/public/sra/SRR6515506.fasta"

    conf = pyspark.SparkConf().setAppName("MetaToolsOmic")
    conf = (conf.setMaster('local[*]')\
           .set('spark.executor.memory', '187G')\
           .set('spark.driver.memory', '187G')\
           .set('spark.driver.maxResultSize', '187G'))
    sc = pyspark.SparkContext(conf=conf)

    pst = time.time()

    content = sc.parallelize(list(SeqIO.parse(filename, "fasta")))

    pend = time.time()
    print("Reading time: ", pend - pst)
   
    kmers = np.array(content.map(convert).map(k_mer).collect())

    pend = time.time()
    print("k-mer processing time: ", pend - pst)

    pst = time.time()

    data = prepare_data(kmers)
    print(data)

    pend = time.time()
    print("Prepare data time: ", pend - pst)

    gend = time.time()
    print("Walltime: ", gend - gst)
