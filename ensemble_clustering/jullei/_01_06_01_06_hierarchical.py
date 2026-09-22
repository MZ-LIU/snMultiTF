import os
os.environ["OMP_NUM_THREADS"] = "2"

import matplotlib.pyplot as plt
import numpy as np
import os,itertools
from pathlib import Path 
import time
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering






def cengcicsp(data,n_cluster,linkage_type) :

    hc = AgglomerativeClustering(n_clusters=n_cluster, linkage=linkage_type)
    labels_pred = hc.fit_predict(data)

    return labels_pred


if __name__ == '__main__' :
    a = '_10_5/2_csp/' + Path(r'data_set\s2.txt').stem
   
    dat = read_data('Dataset/data_s/s-set2.txt')[:,:-1]
    time0 = time.time()
    cc = k_meanscsp(dat,15)
    print(cc)
    cc = cengcicsp(dat,15,'single')
    print(cc)
    time1 = time.time()
    print(time1 - time0)
    # print(cc[1:])

