# This module create the plot for the project

import numpy as np
import matplotlib.pyplot as plt


# plot_results - create the plot (create clusters.pdf)
# input:    K - the K used in the clustering
#           K_used_in_blobs - the K used while generating the points
#           N - the number of points
#           dim - the dimension of the points
#           matrices1 - list of clusters, where each cluster is numpy 2darray (list of points in the cluster), according to Spectral algorithm
#           matrices2 - list of clusters, where each cluster is numpy 2darray (list of points in the cluster), according to K-Means algorithm
#           list_of_points - list of the points
#           sklearn_clusters - the default clustering from make_blobs
#           sc_clusterIndecies - ndarray s.t if list_of_points[i] is in cluster matrices1[j] then sc_clusterIndecies[i] == j
#           km_clusterIndecies - ndarray s.t if list_of_points[i] is in cluster matrices2[j] then km_clusterIndecies[i] == j
# output:   There is no output
def plot_results(K, K_used_in_blobs, N, dim, sc_Clusters, km_Clusters, list_of_points, sklearn_clusters, sc_clusterIndecies, km_clusterIndecies):
    fig = plt.figure(figsize=(10, 7))

    fig.suptitle("Normalized Spectral Clustering VS. K-Means")  # The title of the figure

    bottom_text = f"Data was generated from values: n = {N}, k = {K_used_in_blobs}\n"
    bottom_text += f"The k that was used for both algorithms was {K}\n"

    # Calculating the jaccard measures
    jms, jmk = JaccardMeasure(sc_clusterIndecies, km_clusterIndecies, sklearn_clusters, N)

    bottom_text += f"The Jaccard measure for Spectral Clustering: {jms}\n"
    bottom_text += f"The Jaccard measure for K-Means: {jmk}"

    # Add bottom_text to the figure
    plt.figtext(0.5, 0.01, bottom_text, ha='center', fontsize=16, va='top', fontstyle='oblique')

    # Output the graphs
    if dim == 2:
        # First subplot
        subplot1 = fig.add_subplot(121)

        for cluster in sc_Clusters:
            subplot1.scatter(cluster[:, 0], cluster[:, 1], cmap='paired')

        subplot1.set_title("Normalized Spectral Clustering")
        subplot1.set_xlabel('X')
        subplot1.set_ylabel('Y')

        # Second subplot
        subplot2 = fig.add_subplot(122)

        for cluster in km_Clusters:
            subplot2.scatter(cluster[:, 0], cluster[:, 1], cmap='paired')

        subplot2.set_title("K-Means")
        subplot2.set_xlabel('X')
        subplot2.set_ylabel('Y')
    else:
        # First subplot
        subplot1 = fig.add_subplot(121, projection='3d')

        for cluster in sc_Clusters:
            subplot1.scatter(cluster[:, 0], cluster[:, 1], cluster[:, 2], cmap='paired')

        subplot1.set_title("Normalized Spectral Clustering")
        subplot1.set_xlabel('X')
        subplot1.set_ylabel('Y')

        # Second subplot
        subplot2 = fig.add_subplot(122, projection='3d')

        for cluster in km_Clusters:
            subplot2.scatter(cluster[:, 0], cluster[:, 1], cluster[:, 2], cmap='paired')

        subplot2.set_title("K-Means")
        subplot2.set_xlabel('X')
        subplot2.set_ylabel('Y')

    plt.tight_layout()
    plt.savefig("clusters.pdf", bbox_inches='tight', dpi=100)


# JaccardMeasure - claculating the jaccard measures
# input:    sc_clusterIndecies - 
#           km_clusterIndecies - 
#           sklearn_clusters - the default clustering from make_blobs
#           N - number of points
# output:   (jms, jmk) - tuple of the jaccard measures
def JaccardMeasure(sc_clusterIndecies, km_clusterIndecies, sklearn_clusters, N):
    num_of_pairs_spectral_and_sklearn = 0
    num_of_pairs_spectral_or_sklearn = 0

    num_of_pairs_kmeans_and_sklearn = 0
    num_of_pairs_kmeans_or_sklearn = 0

    in_spectral_or_sklearn = False
    in_kmeans_or_sklearn = False

    for i in range(N):
        for j in range(i + 1, N):
            if sklearn_clusters[i] == sklearn_clusters[j]:
                num_of_pairs_spectral_or_sklearn += 1
                num_of_pairs_kmeans_or_sklearn += 1
                in_spectral_or_sklearn = True
                in_kmeans_or_sklearn = True

            if sc_clusterIndecies[i] == sc_clusterIndecies[j]:
                if in_spectral_or_sklearn == True:
                    num_of_pairs_spectral_and_sklearn += 1
                if in_spectral_or_sklearn == False:
                    num_of_pairs_spectral_or_sklearn += 1

            if km_clusterIndecies[i] == km_clusterIndecies[j]:
                if in_kmeans_or_sklearn == True:
                    num_of_pairs_kmeans_and_sklearn += 1
                if in_kmeans_or_sklearn == False:
                    num_of_pairs_kmeans_or_sklearn += 1

            in_spectral_or_sklearn = False
            in_kmeans_or_sklearn = False

    jms = num_of_pairs_spectral_and_sklearn / num_of_pairs_spectral_or_sklearn
    jmk = num_of_pairs_kmeans_and_sklearn / num_of_pairs_kmeans_or_sklearn
    
    return (jms, jmk)