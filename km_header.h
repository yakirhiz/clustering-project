#ifndef KMEANS_H
#define KMEANS_H

#define EPSILON 0.0001
#define MALLOC_CHECK(x) if(x==NULL){puts("\nMalloc Fail\n"); exit(0);}

double *kmeans(int K, int N, int d, int MAX_ITER, double *N_observations, double *cents);
void printCentroids(int K, int d, double *cents);
double squared_euclidean_distance(double vec1[], double vec2[], int size);
int findClosestCent(int d, int K, double *cents, double *obs);
double *calcNewCentroids(int K, int N, int d, double *cents, double *N_observation);
int centsChanged(int K, int d, double *cents, double *newCents);
double *readFile(int d, int N, char *path);
void printMat(const double *mat, int n, int d);
double *readStdin(int N, int d);
double *initalizeCentroids(int K, int d, double *N_observations);

#endif /* KMEANS_H */
