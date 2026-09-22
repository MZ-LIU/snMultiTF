import scipy.spatial.distance as distance
import time
import numpy as np

try:
    import torch
    HAS_GPU = torch.cuda.is_available()
except ImportError:
    torch = None
    HAS_GPU = False

class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
    
    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]
    
    def union(self, x, y):
        rootX = self.find(x)
        rootY = self.find(y)
        if rootX != rootY:
            self.parent[rootY] = rootX
def cspcluster_pro(result_dist0):
    # time0 = time.time()
    N = int(np.sqrt(2 * len(result_dist0))) + 1

   
    if HAS_GPU:
        device = torch.device("cuda")
        
        tri_idx = torch.triu_indices(N, N, offset=1, device=device)
        u_gpu = tri_idx[0].to(torch.int32)
        v_gpu = tri_idx[1].to(torch.int32)
        dist_gpu = torch.as_tensor(result_dist0, dtype=torch.float32, device=device)

       
        sorted_indices = torch.argsort(dist_gpu)
        u_sorted = u_gpu[sorted_indices]
        v_sorted = v_gpu[sorted_indices]
        dist_sorted = dist_gpu[sorted_indices]

        
        u_cpu = u_sorted.cpu().numpy()
        v_cpu = v_sorted.cpu().numpy()
        dist_cpu = dist_sorted.cpu().numpy()
    else:
       
        idx_upper = np.triu_indices(N, k=1)
        u, v = idx_upper
        sorted_indices = np.argsort(result_dist0)
        u_cpu = u[sorted_indices]
        v_cpu = v[sorted_indices]
        dist_cpu = result_dist0[sorted_indices]

    # Initialize union-find structure
    uf = UnionFind(N)
    result_dij = []

    # Kruskal's algorithm: add edges while avoiding cycles 
    for i in range(len(dist_cpu)):
        ui, vi, d = int(u_cpu[i]), int(v_cpu[i]), float(dist_cpu[i])
        if uf.find(ui) != uf.find(vi):
            uf.union(ui, vi)
            result_dij.append([ui, vi, d])
        if len(result_dij) == N - 1:
            break
    return np.array(result_dij)
def culcalate_level_pro(n_data, data_dij_cul):
    '''
    # Calculate level
    '''
    time0 = time.time()

    if HAS_GPU:
        device = torch.device("cuda")

        data_gpu = torch.as_tensor(data_dij_cul, dtype=torch.float32, device=device)
        u_idx = data_gpu[:, 0].to(torch.long)
        v_idx = data_gpu[:, 1].to(torch.long)

        result_n_n_2 = torch.zeros(n_data, dtype=torch.float32, device=device)
        ones = torch.ones(len(data_dij_cul), dtype=torch.float32, device=device)
        result_n_n_2.index_add_(0, u_idx, ones)
        result_n_n_2.index_add_(0, v_idx, ones)

        result_n_n_2_copy = result_n_n_2.clone()
        count_temp_7_31 = len(data_dij_cul)
        result_8_2_cul = torch.zeros(count_temp_7_31, dtype=torch.float32, device=device)
        temp_jibie = 0

        while count_temp_7_31 != 0:
            temp_jibie += 1
          
            mask = (result_8_2_cul == 0) & ((result_n_n_2_copy[u_idx] == 1) |
                                            (result_n_n_2_copy[v_idx] == 1))
            result_8_2_cul[mask] = float(temp_jibie)
            count_temp_7_31 -= int(torch.sum(mask))

            update_mask = (result_8_2_cul == float(temp_jibie))
            ones_update = torch.ones(int(torch.sum(update_mask)), dtype=torch.float32, device=device)
            result_n_n_2_copy.index_add_(0, u_idx[update_mask], -ones_update)
            result_n_n_2_copy.index_add_(0, v_idx[update_mask], -ones_update)

        result_cpu = result_8_2_cul.cpu().numpy()
    else:
       
        result_n_n_2 = np.zeros(n_data, dtype=np.float32)
        np.add.at(result_n_n_2, data_dij_cul[:, 0].astype(int), 1)
        np.add.at(result_n_n_2, data_dij_cul[:, 1].astype(int), 1)

        result_n_n_2_copy = result_n_n_2.copy()
        count_temp_7_31 = len(data_dij_cul)
        result_8_2_cul = np.zeros(count_temp_7_31, dtype=np.float32)
        temp_jibie = 0

        while count_temp_7_31 != 0:
            temp_jibie += 1
            mask = (result_8_2_cul == 0) & ((result_n_n_2_copy[data_dij_cul[:, 0].astype(int)] == 1) |
                                            (result_n_n_2_copy[data_dij_cul[:, 1].astype(int)] == 1))
            result_8_2_cul[mask] = temp_jibie
            count_temp_7_31 -= np.sum(mask)

            update_mask = (result_8_2_cul == temp_jibie)
            np.subtract.at(result_n_n_2_copy, data_dij_cul[update_mask, 0].astype(int), 1)
            np.subtract.at(result_n_n_2_copy, data_dij_cul[update_mask, 1].astype(int), 1)
        result_cpu = result_8_2_cul

    time1 = time.time()
    # print('level', time1 - time0)
    return result_cpu

def culcalate_label(len_data,data_dij_cula) :
    '''
   
    '''
    data_cluster = []
    for i in range(len(data_dij_cula)) :
        temp_9_in = [-1,-1]
        for j in range(len(data_cluster)) :
            # print(data_dij[i,0],'bnm',data_cluster[j])
            if data_dij_cula[i,0] in data_cluster[j] :
                temp_9_in[0] = j
            if data_dij_cula[i,1] in data_cluster[j] :
                temp_9_in[1] = j
        if temp_9_in[0] == -1 and temp_9_in[1] == -1:
            data_cluster.append([data_dij_cula[i,0],data_dij_cula[i,1]])
        elif temp_9_in[0] != -1 and temp_9_in[1] == -1 :
            data_cluster[temp_9_in[0]].append(data_dij_cula[i,1])
        elif temp_9_in[0] == -1 and temp_9_in[1] != -1 :
            data_cluster[temp_9_in[1]].append(data_dij_cula[i,0])
        elif temp_9_in[0] != -1 and temp_9_in[1] != -1 :
            min_temp = min(temp_9_in)
            max_temp = max(temp_9_in)
            data_cluster.append(data_cluster[min_temp][:] + data_cluster[max_temp][:])
            data_cluster.pop(max_temp)
            data_cluster.pop(min_temp)


    labels = np.zeros((len_data), dtype=np.int32)
    label_count = 0  # label value
    # label_count1 = 0 
    for i in range(len(data_cluster)) :
        for j in range(len(data_cluster[i])) :
            labels[int(data_cluster[i][j])] = label_count
        label_count += 1
    return labels

def cluster_demo_pro(data,xunhuancishu) :
    '''
    For clustering, the input is a two-dimensional array numpy pure python
    '''
    xunhuancishu -= 1

    if HAS_GPU:
        device = torch.device("cuda")
        data_torch = torch.as_tensor(data, dtype=torch.float32, device=device)
        
        dmat = torch.cdist(data_torch, data_torch, p=2)
        tri_idx = torch.triu_indices(dmat.shape[0], dmat.shape[1], offset=1, device=device)
        result_dist0 = dmat[tri_idx[0], tri_idx[1]].cpu().numpy()
    else:
        result_dist0 = distance.pdist(data, 'euclidean')

    data_dij = cspcluster_pro(result_dist0)

    draw_level = []
    draw_level1 = []
    draw_level2 = []
    beilv = [] #

    if HAS_GPU:
        device = torch.device("cuda")
        data_dij_gpu = torch.as_tensor(data_dij, dtype=torch.float32, device=device)

        for temp_8_4 in range(xunhuancishu) :
           
            result_8_2 = culcalate_level_pro(len(data), data_dij_gpu.cpu().numpy())
            result_8_2_gpu = torch.as_tensor(result_8_2, dtype=torch.float32, device=device)

            temp_9_5_gpu = result_8_2_gpu * data_dij_gpu[:, 2]

            max_idx = int(torch.argmax(temp_9_5_gpu))
            max_val = float(temp_9_5_gpu[max_idx])

            beilv.append(max_val)
            draw_level.append(max_val)
            draw_level1.append(float(result_8_2_gpu[max_idx]))
            draw_level2.append(float(data_dij_gpu[max_idx, 2]))

           
            mask = torch.ones(len(data_dij_gpu), dtype=torch.bool, device=device)
            mask[max_idx] = False
            data_dij_gpu = data_dij_gpu[mask]

        final_data_dij = data_dij_gpu.cpu().numpy()
    else:
        for temp_8_4 in range(xunhuancishu) :
            result_8_2 = culcalate_level_pro(len(data),data_dij)
            temp_9_5 = (result_8_2 )*data_dij[:,2]

            count_xunhuan_shanchu = np.argmax(temp_9_5)
            beilv.append(np.max(temp_9_5))
            draw_level.append(np.max(temp_9_5))
            draw_level1.append(result_8_2[count_xunhuan_shanchu])
            draw_level2.append(data_dij[count_xunhuan_shanchu,2])
            data_dij = np.delete(data_dij,count_xunhuan_shanchu, axis=0)
        final_data_dij = data_dij

    labels = culcalate_label(len(data), final_data_dij)
    print('clusterr_demo_pro finished')
    return labels, beilv
