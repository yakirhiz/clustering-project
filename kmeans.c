#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "km_header.h"

double *calcNewCentroids(int K, int N, int d, double *cents, double *N_obs) {
    int i, j, closestCent;
    double *newCents, *obs;
    int *obsInCluster;

    newCents = (double*)calloc((K*d), sizeof(double));
    assert(newCents != NULL);
    obsInCluster = (int*)calloc(K, sizeof(int));
    assert(obsInCluster != NULL);
/*  for each obs in N_observations:
 *  -find the closest centroid (closestCent)
 *  -add obs to newCents[closestCent]
 *  -incerement obsInCluster[closestCent] by 1 */
    for (i=0; i<N; i++){
        obs = &N_obs[i*d];
        closestCent = findClosestCent(d, K, cents, obs);
        for (j=0; j<d; j++) {
            newCents[closestCent*d + j] += N_obs[i*d + j];
        }
        obsInCluster[closestCent] += 1;
    }
/*   calculate centroids by dividing by |S_j|, while making sure we dont devide by zero*/
    for (i=0; i<K; i++){
        if (obsInCluster[i] == 0){
            continue;
        }
        for (j=0; j<d; j++) {
            newCents[i*d + j] = newCents[i*d + j]/obsInCluster[i];
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

    dist = euclidDist(&cents[0], obs, d);
    for (i=1; i<K; i++) {
        newDist = euclidDist(&cents[i*d], obs, d);
        if (newDist < dist){
            dist = newDist;
            index = i;
        }
    }
    return index;
}

/* input:   2 vectors of the same size
 * returns: the euclidean distance between vec1 and vec2 */
double euclidDist(double vec1[], double vec2[], int size){
    int i;
    double sum = 0;

    for (i=0; i<size; i++){
        sum += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);
    }
    return sum;
}

int centsChanged(int K, int d, double *cents, double *newCents){
    int diff = 0;
    int i = 0;
    int j = 0;

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

/* Prints a Matrix of dimentions kxd*/
void printCentroids(int K, int d, double *cents){
    int i;
    int j;

    for (i=0; i<K; i++){
        for (j=0; j<d; j++){
            /* we reached a new row*/
            if (j == d-1){
                printf("%f\n", cents[i*d + j]);
            } else {
                printf("%f,", cents[i*d + j]);
            }
        }
    }
}

double *kmeans(int K, int N, int d, int MAX_ITER, double *N_obs, double *cents) {
    double *newCents;
    int curr_iter = 0;
    /* 
     * K – the number of clusters required.
     * N – the number of observations in the file
     * d – the dimension of each observation and initial centroids
     * MAX_ITER – the maximum number of iterations of the K-means algorithm
     */

    while (curr_iter < MAX_ITER){
        curr_iter += 1;
/*      calculate new centroids*/
        newCents = calcNewCentroids(K, N, d, cents, N_obs);
/*      check if cluster centroids change from previous iteration*/
        if (!centsChanged(K, d, cents, newCents)){
            free(newCents);
            break;
        }
        free(cents);
        cents = newCents;
    }

    free(N_obs);
    return cents;
}


// reads observations from a file (same format as previous hw)
// used for debugging
double *readFile(int d, int N, char *path) {
    double *N_observations;
    int i, j;

/*  N_observations = matrix of size N*d*/
    N_observations = (double*) malloc((d*N) * sizeof(double));
    MALLOC_CHECK(N_observations)

    if (N_observations  == NULL){
        puts("\nProblem in reading Observations file\n");
    }

    FILE *f = fopen(path, "r");

    /* Fill N_observations according to the input file */
    for (i = 0; i < N; i++) {
        for (j = 0; j < d; j++) {
            fscanf(f ,"%lf,", &N_observations[i*d + j]);
        }
    }
    fclose(f);
    return N_observations;
}

/* Prints a Matrix of dimentions nxd*/
void printMat(const double *mat, int n, int d){
    int i;
    int j;

    for (i=0; i<n; i++){
        for (j=0; j<d; j++){
            /* we reached a new row*/
            if (j == d-1){
                printf("%lf\n", mat[i*d + j]);
            } else {
                printf("%lf,", mat[i*d + j]);
            }
        }
    }
}

/* input:   N_observations = matrix of size Nxd
 * returns: cents = matrix of size Kxd
 */
double *initalizeCentroids(int K, int d, double *N_observations) {
    int i, j;
    double *cents;

    cents = (double*) malloc((K*d) * sizeof(double));
    assert(cents != NULL);
    /* initialize cents as the first k observations*/
    for (i = 0; i < K; i++) {
        for (j = 0; j < d; j++) {
            cents[i*d + j] = N_observations[i*d + j];
        }
    }

    return cents;
}

/* input:   none, expects 4 cmd arguments and an input stream
 * returns: 2d matrix (implemented as 1d) where:
 * matrix[i] is pointer to the i-th row i.e. the i-th observation.
 * matrix[i][j] is the j-th coord of the i-th row.
 */
double *readStdin(int N, int d) {
    double *N_observations;
    int i, j;

/*  N_observations = matrix of size N*d*/
    N_observations = (double*) malloc((d*N) * sizeof(double));
    assert(N_observations != NULL);

    /* Fill N_observations according to the input file */
    for (i = 0; i < N; i++) {
        for (j = 0; j < d; j++) {
            scanf("%lf,", &N_observations[j + i*d]);
        }
    }

    return N_observations;
}

int main(int argc, char **argv) {
    int K, N, d, MAX_ITER;
    double *newCents, *N_observations, *cents;
    int curr_iter = 0;
    assert(argc == 5);
    /* Command line arguments:
     * K – the number of clusters required.
     * N – the number of observations in the file
     * d – the dimension of each observation and initial centroids
     * MAX_ITER – the maximum number of iterations of the K-means algorithm
     */

    /* Does not check number of arguments or type */
    K = atoi(argv[1]);
    N = atoi(argv[2]);
    d = atoi(argv[3]);
    MAX_ITER = atoi(argv[4]);

    N_observations = readStdin(N, d);
    cents = initalizeCentroids(K, d, N_observations);
    while (curr_iter < MAX_ITER){
        curr_iter += 1;
/*      calculate new centroids*/
        newCents = calcNewCentroids(K, N, d, cents, N_observations);
/*      check if cluster centroids change from previous iteration*/
        if (!centsChanged(K, d, cents, newCents)){
            break;
        }
        free(cents);
        cents = newCents;
        printf("\n");
    }

    printCentroids(K, d, cents);
    free(N_observations);
    return 0;
}