import mySpecClust as sc
import kmeans_pp as km
import numpy as np

# This function implements Spectral Clustering
# input:    points - ndarray of n points of dimension d
#           N - number of points
#           d - dimension of each point
#           K - number of clusters
# output:   a list "clusters" of size K s.t. clusters[i] = an ndarray of points in the cluster
#           clusterIndecies is an ndarray s.t if points[i] is closet to SC_cents[j] then clusterIndecies[i] == j
def cluster(points, N, d, K=None):

    # execute the first 5 steps of the algorithm
    if (K==None):
        T_list = sc.specClust_c(points.tolist(), N, d, 0) # matrix T from step 5 if the algo (as a list)
        K = len(T_list[0])  # k from eigengap heuristic
    else:
        T_list = sc.specClust_c(points.tolist(), N, d, K)

    T_reshaped = np.array(T_list)

    # step 6 in Spectral Clustering
    MAX_ITER = 300 # as described in the bottom of page 5
    SC_cents_list, _ = km.k_means_pp(K, N, K, MAX_ITER, T_reshaped) #cents for clustering of matrix T (points in R^k)

    SC_cents = np.array(SC_cents_list)

    # step 7 in Spectral Clustering
    # for each point in T_reshaped, find the index of the closest centroid.
    clusterIndecies = np.argmin(np.linalg.norm(T_reshaped[:, None, :] - SC_cents[None, :, :], axis=2), axis=1)

    # clusters[i] is an ndarray of all rows in points that correspond to the indexing of T_reshaped 
    clusters = [points[i==clusterIndecies] for i in range(K)]

    return clusters, clusterIndecies
