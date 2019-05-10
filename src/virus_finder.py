
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


def fuzzy_train(tain_data):
    fpcs = list()
    for k in range(2,10,1):
        print("Training: k = ", k)
        cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(data, k, 2, error=0.005, maxiter=1000)
        fpcs.append((k,fpc))

    print(fpcs)

#def som(data):
#    sm = SOMFactory().build(data, normalization = 'var', initialization='random')
#    sm.train(n_job=48, verbose=False, train_rough_len=2, train_finetune_len=5)
#    pickle.dump(sm, open( "sompy_model.pickle", "wb" ))
#    #topographic_error = sm.calculate_topographic_error()
#    #quantization_error = np.mean(sm._bmu[1])
#    #print ("Topographic error = %s; Quantization error = %s" % (topographic_error, quantization_error))


def prepare_data(kmers):
    dim = [kmers.shape[0], len(max(kmers, key=len))]
    print(dim)
    data = np.zeros(dim, dtype=np.int)
    for i in range(dim[0]):
        data[i][0:len(kmers[i])] = kmers[i][:]
    return data

# Transforma os k-mers em índices na base 10
def k_mer(read):
    k = 16
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

    filename = "/home/fabricio/Documents/projects/mackenzie/data_mining/SRX5784792.fasta"

    conf = pyspark.SparkConf().setAppName("MetaToolsOmic")
    conf = (conf.setMaster('local[*]')\
           .set('spark.executor.memory', '4G')\
           .set('spark.driver.memory', '4G')\
           .set('spark.driver.maxResultSize', '4G'))
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
#    print(data)

    pend = time.time()
    print("Prepare data time: ", pend - pst)

#    pst = time.time()
#
#    som(data)
#
#    pend = time.time()
#    print("SOM time: ", pend - pst)

    pst = time.time()

    fuzzy_train(data)

    pend = time.time()
    print("Fuzzy time: ", pend - pst)

    gend = time.time()
    print("Walltime: ", gend - gst)
