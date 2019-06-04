
# coding: utf-8

import numpy as np
from math import floor
from multiprocessing import Process, Manager
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
    rd = np.zeros(len(seq), dtype=np.int)
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


def run_rest(reads, k, kmers_list, stt):
    for i in range(len(reads)):
        idx = stt + i
        kmers_list[idx] = k_mer_freq_count(k_mer((convert(reads[i])), k), k)
        
def run(reads, k, tid, offset, kmers_list):
    for i in range(len(reads)):
        idx = (tid * offset) + i
        kmers_list[idx] = k_mer_freq_count(k_mer((convert(reads[i])), k), k)
   

def kmer_count(filename, k, context):
    nt = 24
    manager = Manager()

    print('Reading file')
    reads = parse_fasta(filename)

    print('K-mer processing')
    offset = floor(len(reads)/nt)
    rest = len(reads) - (offset * nt)
    kmers_list = manager.list([None] * len(reads))
    p = [None] * nt
    for i in range(nt):
        stt = i * offset
        stp = stt + offset
        p[i] = Process(target=run, args=(reads[stt:stp], k, i, offset, kmers_list, ))
        p[i].start()

    stt = offset * nt
    run_rest(reads[stt:], k, kmers_list, stt)

#    kmers_list = list(kmers_list)

    print('Sparse matrix composition')
    kmers = lil_matrix((len(kmers_list), 4**k), dtype=np.int)
    for i in range(0,len(kmers_list)):
        kmers[i, kmers_list[i][0]] = kmers_list[i][1]
#        print(kmers_list[i][0], kmers_list[i][1])

    for i in range(nt):
        p[i].join()

    return kmers
