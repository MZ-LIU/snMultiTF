import scipy.spatial.distance as distance
from scipy.spatial.distance import squareform
import time
import numpy as np
import torch


GPU_AVAILABLE = torch.cuda.is_available()

class UnionFind:
    def __init__(self, n):
        self.parent = np.arange(n)

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]] 
            x = self.parent[x]
        return x

    def union(self, x, y):
        rootX = self.find(x)
        rootY = self.find(y)
        if rootX != rootY:
            self.parent[rootY] = rootX
            return True
        return False

def cspcluster_pro(result_dist0):
    """Kruskal builds a minimum spanning tree (MST)"""
    N = int(np.sqrt(2 * len(result_dist0))) + 1

    idx_upper = np.triu_indices(N, k=1)
    u, v = idx_upper

    if GPU_AVAILABLE:
        dist_t = torch.from_numpy(result_dist0.astype(np.float32)).cuda()
        sorted_indices = torch.argsort(dist_t).cpu().numpy()
        del dist_t
    else:
        sorted_indices = np.argsort(result_dist0)

    u_sorted = u[sorted_indices]
    v_sorted = v[sorted_indices]
    dist_sorted = result_dist0[sorted_indices]

    uf = UnionFind(N)
    result_dij = []

    
    for i in range(len(dist_sorted)):
        node_u, node_v, d = u_sorted[i], v_sorted[i], dist_sorted[i]
        if uf.union(node_u, node_v):
            result_dij.append([node_u, node_v, d])
            if len(result_dij) == N - 1:
                break

    return np.array(result_dij)

def culcalate_level_pro(n_data, data_dij_cul):
    """Compute level (GPU parallel optimization)"""
    if not GPU_AVAILABLE:
     
        result_n_n_2 = np.zeros(n_data, dtype=np.float32)
        np.add.at(result_n_n_2, data_dij_cul[:, 0].astype(int), 1)
        np.add.at(result_n_n_2, data_dij_cul[:, 1].astype(int), 1)

        result_n_n_2_copy = result_n_n_2.copy()
        count_temp_7_31 = len(data_dij_cul)
        result_8_2_cul = np.zeros(count_temp_7_31, dtype=np.float32)
        temp_jibie = 0

        u_idx = data_dij_cul[:, 0].astype(int)
        v_idx = data_dij_cul[:, 1].astype(int)

        while count_temp_7_31 != 0:
            temp_jibie += 1
            mask = (result_8_2_cul == 0) & ((result_n_n_2_copy[u_idx] == 1) | (result_n_n_2_copy[v_idx] == 1))
            result_8_2_cul[mask] = temp_jibie
            count_temp_7_31 -= np.sum(mask)

            update_mask = (result_8_2_cul == temp_jibie)
            np.subtract.at(result_n_n_2_copy, u_idx[update_mask], 1)
            np.subtract.at(result_n_n_2_copy, v_idx[update_mask], 1)
        return result_8_2_cul

   
    u_idx = torch.from_numpy(data_dij_cul[:, 0].astype(np.int64)).cuda()
    v_idx = torch.from_numpy(data_dij_cul[:, 1].astype(np.int64)).cuda()

    
    degrees = torch.zeros(n_data, dtype=torch.float32, device='cuda')
    degrees.scatter_add_(0, u_idx, torch.ones_like(u_idx, dtype=torch.float32))
    degrees.scatter_add_(0, v_idx, torch.ones_like(v_idx, dtype=torch.float32))

    degrees_copy = degrees.clone()
    n_edges = len(data_dij_cul)
    result_levels = torch.zeros(n_edges, dtype=torch.float32, device='cuda')
    temp_jibie = 0

    while n_edges > 0:
        temp_jibie += 1
       
        mask = (result_levels == 0) & ((degrees_copy[u_idx] == 1) | (degrees_copy[v_idx] == 1))

        count_mask = mask.sum()
        if count_mask == 0: break # 
        result_levels[mask] = temp_jibie
        n_edges -= count_mask.item()

     
        update_u = u_idx[mask]
        update_v = v_idx[mask]
        degrees_copy.scatter_add_(0, update_u, -torch.ones_like(update_u, dtype=torch.float32))
        degrees_copy.scatter_add_(0, update_v, -torch.ones_like(update_v, dtype=torch.float32))

    res = result_levels.cpu().numpy()
    del u_idx, v_idx, degrees, degrees_copy, result_levels
    return res

def culcalate_label(len_data, data_dij_cula):
   
    uf = UnionFind(len_data)
    for edge in data_dij_cula:
        uf.union(int(edge[0]), int(edge[1]))

    
    final_labels = np.array([uf.find(i) for i in range(len_data)])
    unique_labels = np.unique(final_labels)
    label_map = {old: new for new, old in enumerate(unique_labels)}
    return np.array([label_map[l] for l in final_labels])

def cluster_demo_pro(data, xunhuancishu, input_type='coordinates'):
    xunhuancishu -= 1
    n_samples = data.shape[0]

    if input_type == 'consensus':
        distance_matrix = 1.0 - data
        result_dist0 = squareform(distance_matrix)
    else:
       
        if GPU_AVAILABLE and n_samples < 20000:
            data_t = torch.from_numpy(data.astype(np.float32)).cuda()
            dist_t = torch.pdist(data_t)
            result_dist0 = dist_t.cpu().numpy()
            del data_t, dist_t
        else:
            result_dist0 = distance.pdist(data, 'euclidean')

    data_dij = cspcluster_pro(result_dist0)

    beilv = []
    
    if GPU_AVAILABLE:
        weights_t = torch.from_numpy(data_dij[:, 2].astype(np.float32)).cuda()

    for temp_8_4 in range(xunhuancishu):
       
        result_8_2 = culcalate_level_pro(n_samples, data_dij)

        if GPU_AVAILABLE:
            levels_t = torch.from_numpy(result_8_2.astype(np.float32)).cuda()
            temp_9_5_t = levels_t * weights_t

            count_xunhuan_shanchu = torch.argmax(temp_9_5_t).item()
            max_val = temp_9_5_t[count_xunhuan_shanchu].item()
            beilv.append(max_val)

          
            data_dij = np.delete(data_dij, count_xunhuan_shanchu, axis=0)
            weights_t = torch.cat([weights_t[:count_xunhuan_shanchu], weights_t[count_xunhuan_shanchu+1:]])
            del levels_t, temp_9_5_t
        else:
            temp_9_5 = result_8_2 * data_dij[:, 2]
            count_xunhuan_shanchu = np.argmax(temp_9_5)
            beilv.append(np.max(temp_9_5))
            data_dij = np.delete(data_dij, count_xunhuan_shanchu, axis=0)

    labels = culcalate_label(n_samples, data_dij)
    print('cluster_demo_pro function execution clustering is completed (GPU acceleration has been applied)')
    return labels, beilv
