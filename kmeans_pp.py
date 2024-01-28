import argparse
import pandas as pd
import numpy as np
import mykmeanssp as km


# This function reads the cmd line.
# returns: tuple (K, N, d, MAX_ITER) where:
# K – the number of clusters required.
# N – the number of observations in the file
# d – the dimension of each observation and initial centroids
# MAX_ITER – the maximum number of iterations of the K-means algorithm
def readargs():
    parser = argparse.ArgumentParser()
    K = parser.add_argument("K", type=int, help="K – Number of clusters required.")
    N = parser.add_argument("N", type=int, help="N – Number of observations in the file.")
    d = parser.add_argument("d", type=int, help="d – Dimension of each observation and centroid.")
    MAX_ITER = parser.add_argument("MAX_ITER", type=int,
                                   help="MAX_ITER – Maximum number of iterations of the K-means algorithm.")
    filename = parser.add_argument("filename", help="filename - Input file with the observations.")
    args = parser.parse_args()
    return args.K, args.N, args.d, args.MAX_ITER, args.filename


# Input:
#   K – Number of clusters required
#   N – Number of observations in the file
#   d – Dimension of each observation and initial centroids
#   MAX_ITER – Maximum number of iterations of the K-means algorithm
#   obs - Observations form the input file
# Returns:
#   a python 1d array of centroids
def k_means_pp(K, N, d, MAX_ITER, obs):
    np.random.seed(0)
    indices = np.empty(K, dtype="int32")
    insert_in_index = 0

    # My Algo
    cents = np.empty(shape=(K, d), dtype="float64")

    # The first centroid is chosen randomly
    index_of_first = np.random.choice(N, 1)[0]
    cents[0] = obs[index_of_first]
    indices[insert_in_index] = index_of_first
    insert_in_index += 1

    distances = np.full(N, np.inf, dtype="float64")

    for j in range(1, K):
        distances = np.minimum(distances, np.sum(np.power(cents[j - 1] - obs, 2), axis=1))
        probs = distances / np.sum(distances)  # Calculate probability of each observation
        index = np.random.choice(N, 1, p=probs)[0]  # Choose index randomly
        cents[j] = obs[index]
        indices[insert_in_index] = index  # Update the array of indices
        insert_in_index += 1

    indices.sort()

    # Call the API
    return km.kmeanspp_c(K, N, d, MAX_ITER, obs.tolist(), cents.tolist()) # returns a python 1d array of centroids



# debug
def main():
    K, N, d, MAX_ITER, filename = readargs()
    df = pd.read_csv(filename, header=None)
    obs = df.to_numpy(dtype="float64")

    if K < 1 or N < 1 or d < 1 or MAX_ITER < 1:
        print("Error in arguments")
    elif N < K:
        print("Error in arguments")
    elif obs.shape != (N, d):
        print("Error in arguments")

    k_means_pp(K, N, d, MAX_ITER, obs)
