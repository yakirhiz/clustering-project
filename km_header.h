#ifndef KMEANS_H
#define KMEANS_H

#define EPSILON 0.0001
#define MALLOC_CHECK(p) \
    do { \
        if ((p) == NULL) { \
            printf("Memory allocation failed\n"); \
            exit(1); \
        } \
    } while (0)

double *kmeans(int K, int N, int d, int MAX_ITER, double *N_observations, double *cents);
void print_centroids(double *cents, int K, int d);
double squared_euclidean_distance(double vec1[], double vec2[], int size);
int findClosestCent(int d, int K, double *cents, double *obs);
double *calcNewCentroids(int K, int N, int d, double *cents, double *N_observation);
int centsChanged(int K, int d, double *cents, double *newCents);
double *readFile(int d, int N, char *path);
double *readStdin(int N, int d);
double *initializeCentroids(int K, int d, double *N_observations);

#endif /* KMEANS_H */
