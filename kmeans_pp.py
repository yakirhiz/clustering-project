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
    parser.add_argument("K", type=int, help="Number of clusters required.")
    parser.add_argument("N", type=int, help="Number of observations in the file.")
    parser.add_argument("d", type=int, help="Dimension of each observation and centroid.")
    parser.add_argument("MAX_ITER", type=int, help="Maximum number of iterations of the K-means algorithm.")
    parser.add_argument("filename", help="Input file with the observations.")
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

    # My Algo
    cents = np.empty(shape=(K, d), dtype="float64")

    index_of_first = np.random.randint(N)  # The first centroid is chosen randomly
    cents[0] = obs[index_of_first]

    distances = np.full(N, np.inf, dtype="float64")

    for j in range(1, K):
        distances = np.minimum(distances, np.sum(np.power(cents[j - 1] - obs, 2), axis=1))
        probs = distances / np.sum(distances)  # Calculate probability of each observation
        index = np.random.choice(N, 1, p=probs)[0]  # Choose index randomly
        cents[j] = obs[index]

    # Call the API
    return km.kmeanspp_c(K, N, d, MAX_ITER, obs.tolist(), cents.tolist()) # returns a python 1d array of centroids


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

    centroids_list = k_means_pp(K, N, d, MAX_ITER, obs)
    print(np.array(centroids_list).reshape(K, d))


if __name__ == '__main__':
    main()
