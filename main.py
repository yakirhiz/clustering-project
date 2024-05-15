import kmeans_pp as km
import SpecClust as sc
import numpy as np
from sklearn.datasets import make_blobs

import plot_module as plotmod
import text_module as textmod


# This file is the glue behind the Spectral Clustering Project


# main - implements Spectral Clustering
# input:    K = the number of clusters
#           N = the number of points
#           rand (bool) = if true draw K and N from (MAX_CAPACITY/2, MAX_CAPACITY)
# output:   an ndarray Nxd of N points of dimension d
def main(K, N, rand):
    max_cap_k = 20
    max_cap_n = 1000
    MAX_ITER = 300
    
    print(f"Max capacity: k={max_cap_k}, n={max_cap_n}.")

    if rand:
        N = np.random.randint(max_cap_n // 2, max_cap_n)
        K = np.random.randint(max_cap_k // 2, min(max_cap_k, N)) # K cannot be larger than N
        
    K_used_in_blobs = K

    # create points to be clustered
    points, sklearn_clusters = generate_random_points(N, K)

    d = points.shape[1] # the dimension of the points

    if rand:
        # cluster with Spectral Clustering
        SC_clusters, sc_clusterIndecies = sc.cluster(points, N, points.shape[1])  # if there is no K, we use the heuristic
        K = len(SC_clusters)
    else:
        # cluster with Spectral Clustering
        SC_clusters, sc_clusterIndecies = sc.cluster(points, N, points.shape[1], K)
        
    # cluster with K-Means
    KM_cents_list = km.k_means_pp(K, N, d, MAX_ITER, points)
    KM_cents = np.reshape(np.array(KM_cents_list), (K, d)) # reshape
    # km_clusterIndecies is an ndarray s.t if points[i] is closet to KM_cents[j] then km_clusterIndecies[i] == j
    km_clusterIndecies = np.apply_along_axis(lambda row: sc.findClosestCentIndex(row, KM_cents), 1, points)  
    KM_clusters = [points[i==km_clusterIndecies,:] for i in range(K)] #fancy indexing, puts ndarrays of points in a python list

    # Visual Output
    plotmod.plot_results(K, K_used_in_blobs, N, d, SC_clusters, KM_clusters, points, sklearn_clusters, sc_clusterIndecies, km_clusterIndecies)
    
    # Textual Output
    textmod.output_results(K, N, d, points, sklearn_clusters, sc_clusterIndecies, km_clusterIndecies)

# generate_random_points - creates the random points to be sent to the clustering algorithms.
#
# input:    N = number of points
#           K = the number of blobs
# output:   an ndarray (Nxd) of N points of dimension d
def generate_random_points(N, K):
    dim = np.random.randint(2, 4)  # The dimension is chosen randomly between 2 and 3.
    return make_blobs(N, dim, centers=K)