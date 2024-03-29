#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "km_header.h"

double *calcNewCentroids(int K, int N, int d, double *cents, double *observations) {
    int i, j, closestCent, *obsInCluster;
    double *newCents, *observation;

    newCents = (double*)calloc((K*d), sizeof(double));
    assert(newCents != NULL);
    obsInCluster = (int*)calloc(K, sizeof(int));
    assert(obsInCluster != NULL);

    for (i=0; i<N; i++){
        observation = &observations[i*d];
        closestCent = findClosestCent(d, K, cents, observation);
        for (j=0; j<d; j++) {
            newCents[closestCent*d + j] += observations[i*d + j];
        }
        obsInCluster[closestCent] += 1;
    }
/*   calculate centroids by dividing by |S_j|, while making sure we dont devide by zero*/
    for (i=0; i<K; i++){
        if (obsInCluster[i] == 0){
            continue;
        }
        for (j=0; j<d; j++) {
            newCents[i*d + j] /= obsInCluster[i];
        }
    }

    free(obsInCluster);
    return newCents;
}

/* returns the index of the closest centroid*/
int findClosestCent(int d, int K, double *cents, double *obs){
    int i;
    int index = 0;
    double dist, newDist;

    dist = squared_euclidean_distance(&cents[0], obs, d);
    for (i=1; i<K; i++) {
        newDist = squared_euclidean_distance(&cents[i*d], obs, d);
        if (newDist < dist){
            dist = newDist;
            index = i;
        }
    }
    return index;
}

/* input:   2 vectors of the same size
 * returns: the euclidean distance between vec1 and vec2 */
double squared_euclidean_distance(double vec1[], double vec2[], int size){
    int i;
    double sum = 0;

    for (i=0; i<size; i++){
        sum += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);
    }
    return sum;
}

int centsChanged(int K, int d, double *cents, double *newCents){
    int diff = 0;
    int i, j;

    for (i=0; i<K; i++){
        for (j=0; j<d; j++){
            diff = cents[i*d + j] - newCents[i*d + j];
            if (diff != 0){ 
                return 1;
            }
        }
    }
    return 0;
}

/* Prints the centroids (K x d) */
void print_centroids(double *cents, int K, int d){
    int i, j;

    for (i=0; i<K; i++){
        for (j=0; j<d; j++){
            if (j == d-1){
                printf("%f\n", cents[i*d + j]);
            } else {
                printf("%f,", cents[i*d + j]);
            }
        }
    }
}

/*
 * K – the number of clusters required.
 * N – the number of observations in the file.
 * d – the dimension of each observation and initial centroids.
 * MAX_ITER – the maximum number of iterations of the K-means algorithm.
 * observations - array containing all observations.
 * cents - array containing initial centroids.
 */
double *kmeans(int K, int N, int d, int MAX_ITER, double *observations, double *cents) {
    double *newCents;
    int iter = 0;

    while (iter < MAX_ITER){
        iter += 1;
        newCents = calcNewCentroids(K, N, d, cents, observations);
        if (!centsChanged(K, d, cents, newCents)){
            free(newCents);
            break;
        }
        free(cents);
        cents = newCents;
    }

    return cents;
}


// Reads observations from a file
double *readFile(int N, int d, char *path) {
    double *observations;
    int i, j;

    observations = (double*) malloc((d*N) * sizeof(double));
    assert(observations != NULL);

    FILE *f = fopen(path, "r");

    for (i = 0; i < N; i++) {
        for (j = 0; j < d; j++) {
            fscanf(f ,"%lf,", &observations[i*d + j]);
        }
    }

    fclose(f);
    return observations;
}

/* input:   observations = matrix of size Nxd
 * returns: cents = matrix of size Kxd
 */
double *initializeCentroids(int K, int d, double *observations) {
    int i, j;
    double *cents;

    cents = (double*) malloc((K*d) * sizeof(double));
    assert(cents != NULL);
    /* initialize cents as the first k observations*/
    for (i = 0; i < K; i++) {
        for (j = 0; j < d; j++) {
            cents[i*d + j] = observations[i*d + j];
        }
    }

    return cents;
}

// Reads observations from stdin
double *readStdin(int N, int d) {
    double *observations;
    int i, j;

    observations = (double*) malloc((d*N) * sizeof(double));
    assert(observations != NULL);

    for (i = 0; i < N; i++) {
        for (j = 0; j < d; j++) {
            scanf("%lf,", &observations[j + i*d]);
        }
    }

    return observations;
}

int main(int argc, char **argv) {
    int K, N, d, MAX_ITER;
    double *observations, *cents;
    assert(argc == 5);
    /* Command line arguments:
     * K – the number of clusters required.
     * N – the number of observations in the file
     * d – the dimension of each observation and initial centroids
     * MAX_ITER – the maximum number of iterations of the K-means algorithm
     */

    K = atoi(argv[1]);
    N = atoi(argv[2]);
    d = atoi(argv[3]);
    MAX_ITER = atoi(argv[4]);

    observations = readStdin(N, d);
    cents = initializeCentroids(K, d, observations);
    cents = kmeans(K, N, d, MAX_ITER, observations, cents);

    print_centroids(cents, K, d);
    free(cents);
    free(observations);
    return 0;
}