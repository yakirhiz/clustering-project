# This module create the text files output for the project
import os
import numpy as np


# output_results - create the text files output (create data.txt & clusters.txt)
# input:    K - the K used in the clustering
#           N - the number of points
#           dim - the dimension of the points
#           list_of_points - list of the points
#           sklearn_clusters - the default clustering from make_blobs
#           sc_clusterIndecies - ndarray s.t if list_of_points[i] is in cluster matrices1[j] then sc_clusterIndecies[i] == j
#           km_clusterIndecies - ndarray s.t if list_of_points[i] is in cluster matrices2[j] then km_clusterIndecies[i] == j
# output:   There is no output
def output_results(K, N, dim, list_of_points, sklearn_clusters, sc_clusterIndecies, km_clusterIndecies):
    # points_of_spectral - list of lists, where each list is represent a cluster, the list will contain the indecies
    #                      of the points (the original indecies) that in each cluster, according to Spectral algorithm
    # points_of_kmeans - list of lists, where each list is represent a cluster, the list will contain the indecies
    #                      of the points (the original indecies) that in each cluster, according to K-Means algorithm
    points_of_spectral = []
    points_of_kmeans = []
    
    for i in range(K):
        points_of_spectral.append((np.where(sc_clusterIndecies == i)[0].tolist()))
        points_of_kmeans.append((np.where(km_clusterIndecies == i)[0].tolist()))

    # data.txt Part
    data_output = ""

    for i in range(N):
        for j in range(dim):
            data_output += str(list_of_points[i][j]) + ","
        if i < N - 1:
            data_output += str(sklearn_clusters[i]) + "\n"
        else:
            data_output += str(sklearn_clusters[i])

    with open(os.path.join("output", "data.txt"), "w") as data_output_file:
        data_output_file.write(data_output)

    # clusters.txt Part
    clusters_output = str(K) + "\n"  # First line contains number of clusters

    # First Algorithm
    for i in range(K):
        length = len(points_of_spectral[i])
        for j in range(length - 1):
            clusters_output += str(points_of_spectral[i][j]) + ","
        if i < K - 1:
            clusters_output += str(points_of_spectral[i][length - 1]) + "\n"
        else:
            clusters_output += str(points_of_spectral[i][length - 1])

    clusters_output += "\n"  # Newline between the results of the different algorithms

    # Second Algorithm
    for i in range(K):
        length = len(points_of_kmeans[i])
        for j in range(length - 1):
            clusters_output += str(points_of_kmeans[i][j]) + ","
        if i < K - 1:
            clusters_output += str(points_of_kmeans[i][length - 1]) + "\n"
        else:
            clusters_output += str(points_of_kmeans[i][length - 1])

    with open(os.path.join("output", "clusters.txt"), "w") as clusters_output_file:
        clusters_output_file.write(clusters_output)
