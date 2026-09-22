import os
os.environ["OMP_NUM_THREADS"] = "2"

import matplotlib.pyplot as plt
import numpy as np

import os,itertools
from pathlib import Path 
import time
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering



from sklearn.metrics import (
    davies_bouldin_score, silhouette_score, calinski_harabasz_score,
    normalized_mutual_info_score, adjusted_rand_score, fowlkes_mallows_score
)

def k_meanscsp(data, n_cluster, random_state=42, input_type='embed'):
    """
    input_type:
        'embed' - Original behavior, data is the (n_cells, n_features) embedding matrix
        'consensus' - data is (n_cells, n_cells) consistency matrix (similarity matrix)
    """
    if input_type == 'consensus':

        X = data  

       

    elif input_type == 'embed':
        X = data
    else:
        raise ValueError("input_type can only be 'embed' or 'consensus'")

    kmeans = KMeans(n_clusters=n_cluster, n_init=10, random_state=random_state)
    labels_pred = kmeans.fit_predict(X)
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



    
