
# coding: utf-8

from Bio import SeqIO
import pyspark
import numpy as np
import time
from sompy.sompy import SOMFactory
from sompy.visualization.mapview import View2DPacked
from sompy.visualization.umatrix import UMatrixView
from matplotlib import pyplot as plt
import gc
import pickle
import skfuzzy as fuzz
from sklearn.cluster import KMeans
from sklearn import decomposition
from kmer_count import kmer_count


def kmeans(data):
    k = 3
    return KMeans(n_clusters=k, random_state=0).fit(data)


def fuzzy(data):
    #fpcs = list()
    k = 3
#    for k in range(2,10,1):
#        print("Training: k = ", k)
    cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(data, k, 2, error=0.0005, maxiter=100000)
#        fpcs.append((k,fpc))
    return cntr, u, u0, d, jm, p, fpc


def som(data):
    mapsize = [15,10]
    som = SOMFactory().build(data, mapsize, normalization = 'var', initialization='random')
    som.train(n_job=4, verbose=False, train_rough_len=2, train_finetune_len=5)
    pickle_out =  open('som.dump', 'wb')
    pickle.dump(som, pickle_out, pickle.HIGHEST_PROTOCOL)
    pickle_out.close()
    #map_labels = som.cluster(n_clusters=4)
    #data_labels = np.array([map_labels[int(k)] for k in som._bmu[0]])
    #print(data_labels.shape)
#    u = UMatrixView(50, 50, 'umatrix', show_axis=True, text_size=8, show_text=True)
#    UMAT  = u.build_u_matrix(som, distance=1, row_normalized=False)
#    UMAT = u.show(som, distance2=1, row_normalized=False, show_data=True, contooor=True, blob=False)
#    v = View2DPacked(2, 2, 'k-means',text_size=8)
#    cl = som.cluster(n_clusters=3)
#    v.show(som, what='cluster')


if __name__ == "__main__":
    gst = time.time()

    query = "/home/fabricio/Documents/projects/mackenzie/sequences/SRX5784792.fasta"
    #db = "/home/fabricio/Documents/projects/mackenzie/sequences/AF033819.3.fasta"

    pst = time.time()

    k = 3
    kmer_query_freq = kmer_count(query, k, "K-merQuery")
    #kmer_db_freq = kmer_count(db, k, "K-merDB")

    pend = time.time()
    print("K-mer freq count time: ", (pend-pst))

    pst = time.time()

    som(kmer_query_freq)

    pend = time.time()
    print("SOM time: ", (pend-pst))

#    pst = time.time()
#
#    cntr, u, u0, d, jm, p, fpc = fuzzy(kmer_query_freq)
#    print(cntr)
#    print(u[:1])
#    print(u0[:1])
#    print(d[:1])
#    print(jm)
#    print(p)
#    print(fpc)
#
#    pend = time.time()
#    print("Fuzzy time: ", (pend-pst))

#    u_2d = decomposition.PCA(n_components=2).fit_transform(u)
#    freq_2d = decomposition.PCA(n_components=2).fit_transform(kmer_query_freq)
#    cntr_2d = decomposition.PCA(n_components=2).fit_transform(cntr)
#    cluster_membership = np.argmax(u_2d, axis=0)  # Hardening for visualization
#    fig, ax = plt.subplots()
#    ax.set_title('Fuzzy clustering')
#    for j in range(3):
#        ax.plot(freq_2d[cluster_membership == j, 0],
#                freq_2d[cluster_membership == j, 1], 'o',
#                label='series ' + str(j))
#    ax.legend()
#    plt.show()

#    pst = time.time()
#    kmeans_result = kmeans(kmer_query_freq)
#    freq_2d = decomposition.PCA(n_components=2).fit_transform(kmer_query_freq)
#    cntr_2d = decomposition.PCA(n_components=2).fit_transform(kmeans_result.cluster_centers_)
#    plt.scatter(freq_2d[:,0], freq_2d[:,1], c= kmeans_result.labels_.astype(float), s=50, alpha=0.5)
#    plt.scatter(cntr_2d[:,0], cntr_2d[:,0], c='red', s=50)
#    plt.show()

    gend = time.time()
    print("Walltime: ", gend - gst)
