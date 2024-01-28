#ifndef SOFTPROJ_FINAL_TOOLS_H
#define SOFTPROJ_FINAL_TOOLS_H

#define EPSILON 0.0001
#define MALLOC_CHECK(x) if(x==NULL){puts("\nMalloc Fail\n"); exit(0);}

//in Algos.c
double *formAdjMat(double *points, int d, int n);
double *formLapMat(double *adj_mat, double *diagdeg_mat, int n);
double **modGrahamSchmidt(double *A, int n);
double **QRIteration(const double *A, int n);
int findK(const double *diag_mat, int n);
void normalizeRows(double *mat, int k, int n);
double *formDiagDegreeMat(double *adj_mat, int n, int pow_minus_half);
double **specClust(double *points, int n, int d, int k_from_user);
int tup_comparefunc (const void * a, const void * b);

//in Tools.c
void printHello(void); //for debugging
double *readFile(int d, int N, char *path); // for debugging
void printMat(const double *mat, int d, int n);
double euclidDist(const double vec1[], const double vec2[], int size);
double sumRow(double *row, int d);
double normColumn(const double *mat, int rows, int col);
double dotCols(const double *A, int colA, const double *B, int colB, int n);
double *multMat(double *mat1,int n1, int d1, double *mat2, int n2, int d2);
int comparefunc (const void * a, const void * b);

#endif //SOFTPROJ_FINAL_TOOLS_H
