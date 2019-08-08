
# coding: utf-8

import sys
import time
import csv
import gc
import pyspark as sp
import numpy as np
import multiprocessing as mp
from pyspark.ml.feature import CountVectorizer
from math import floor
from pyspark.ml.clustering import KMeans
from pyspark.sql.functions import col, udf
from pyspark.ml.evaluation import ClusteringEvaluator


def k_mer(read, k):
    mer = list()
    for i in range(0, len(read['seq'])-k):
        mer.append(read['seq'][i:i+k])
    kmer = {'header': read['header'], 'seq': read['seq'], 'kmer': mer}
    return kmer


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


if __name__ == '__main__':

    filename = sys.argv[1]
    outdir = sys.argv[2]
    k = 10
    ng = 100

    conf = sp.SparkConf().setAppName('VirusFinder')
    conf = (conf.setMaster('local[4]')\
           .set('spark.executor.memory', '8G')\
           .set('spark.driver.memory', '8G')\
           .set('spark.driver.maxResultSize', '8G')\
           .set('spark.network.timeout','100000s')\
           .set('spark.executor.heartbeatInterval', '10000s'))

    print('Reading file')
    st = time.time()
    data = parse_fasta(filename)
    nd = time.time()
    print('Reading time: ', (nd - st))

    st = time.time()
    sc = sp.SparkContext(conf=conf)
    sql = sp.sql.SQLContext(sc)
    partitions = floor((sys.getsizeof(data) / 1024.) / 10.)
    rdd = sc.parallelize(data, partitions)
    kmers = rdd.map(lambda read: k_mer(read, k)).toDF()
    kmers_freq= CountVectorizer(inputCol='kmer', outputCol='features').fit(kmers).transform(kmers)
    print('kmeans')
    kmeans = KMeans().setK(ng).setSeed(ng).fit(kmers_freq).transform(kmers_freq).select('prediction', 'header', 'seq', 'features')
#    print('saving')
#    kmeans.write.parquet(outdir)
#    print('saved')
    nd = time.time()
    print('kmeans processing time: ', (nd-st))
