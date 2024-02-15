#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sc_header.h"


/* This module holds the main part of the Spectural Clustering Algorithm, as described
 * in the homework document, sections 3-5*/

// forms adjacency matrix from a set of points.
// input:   a matrix "points" whose rows are points in R^d space
// output:  adjacency matrix for the points
double *formAdjMat(double *points, int d, int n) {
    int i, j;
    double *p1, *p2, weight;

    double *adj_mat = (double *) calloc((n*n), sizeof(double));
    MALLOC_CHECK(adj_mat)

    //for each point, calc the distance to the rest of the points.
    // we only fill out top half of the mat as it is symmetric (undirected graph)
    for (i=0; i<n; i++){
        adj_mat[i*n + i] = 0;
        p1 = &points[i*d];
        for (j=i+1; j<n; j++){
            p2 = &points[j*d];
            weight = exp((-1.0) * (squared_euclidean_distance(p1, p2, d) / 2));
            adj_mat[i*n + j] = weight;
            adj_mat[j*n + i] = weight;
        }
    }

    return adj_mat;
}

// input:   Adjacentcy matrix, pow_minus_half boolean
// output:  "Diagonal Degree Matrix" diagdeg_mat[i][i] = degree of vertex i in graph represented by adj_mat
//          if pow_minus_half==True, diagdeg_mat[i][i] = 1/sqrt(degree)
double *formDiagDegreeMat(double *adj_mat, int n, int pow_minus_half){
    int i;

    double *diagdeg_mat = (double *) calloc((n*n), sizeof(double));
    MALLOC_CHECK(diagdeg_mat)

    if (pow_minus_half) {
        for (i = 0; i < n; i++) {
            double sumOfRow_i = sumRow(&adj_mat[i * n], n);
            if (sumOfRow_i < 0){
                puts("an error has occured in diagonal degree matrix");
                exit(0);
            } else if (sumOfRow_i == 0){
                diagdeg_mat[i * n + i] = (double) 0;
            } else {
                diagdeg_mat[i * n + i] = 1 / sqrt(sumOfRow_i);
            }
        }
    } else { //mainly for debugging
        for (i=0; i<n; i++){
            diagdeg_mat[i*n + i] = sumRow(&adj_mat[i*n], n);
        }
    }

    return diagdeg_mat;

}

// input:   "Weighted Adjacency Matrix nxn"-adj_mat
//          "Diagonal Degree Matrix nxn"-diagdeg_mat- D^(.5)
// output:  I-DWD where I is the identity matrix
double *formLapMat(double *adj_mat, double *diagdeg_mat, int n) {
    int i, j;

    double *lap_mat = (double *)calloc(n*n, sizeof(double));
    MALLOC_CHECK(lap_mat)

    double *DW = multMat(diagdeg_mat, n, n, adj_mat, n, n);
    double *DWD = multMat(DW, n, n, diagdeg_mat, n, n);

    // calculating lap_mat = I-DWD
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            if (i==j){
                lap_mat[i*n +j] = 1 - DWD[i*n +j];
            } else {
                lap_mat[i*n +j] = (-1)*DWD[i*n +j];
            }
        }
    }

    free(DW);
    free(DWD);
    return lap_mat;
}

// "The Modified Gram-Schmidt Algorithm"
// input:   a matrix "A" nxn
// output:  a tuple, "res", of two points. res[0] = Q, res[1] = R
//          s.t A = QR & Q is orthogonal & R is upper triangular
double **modGrahamSchmidt(double *A, int n){
    int i, j, k, l, n_sq;
    double R_ii, Q_ki;

    double *Q = (double *)calloc(n*n, sizeof(double));
    MALLOC_CHECK(Q)
    double *R = (double *)calloc(n*n, sizeof(double));
    MALLOC_CHECK(R)

    double *U = (double *)malloc(n*n*sizeof(double)); // U = A
    MALLOC_CHECK(U)
    n_sq = n*n;
    for (i=0; i<n_sq; i++){ // copying U=A
        U[i] = A[i];
    }

    for (i=0; i<n; i++){

        R_ii = normColumn(U, n, i); // calculates the norm of the ith column
        R[i*n + i] = R_ii; //R_ii = ||U_i||^2



        for (k=0; k<n; k++){ // Q_i = U_i/R_ii
            if (R_ii != 0){ // make sure we dont devide by zero
                Q_ki = (1/R_ii)*U[k*n + i];
            } else {
                Q_ki = 0;
            }
            Q[k*n + i] = Q_ki;
        }

        for (j=i+1; j<n; j++){
            R[i*n + j] = dotCols(Q, i, U, j, n); // R_ij = Q_i^t * U_j (dot product of columns)

            for (l=0; l<n; l++){ //U_j = U_j - R_ij*Q_i
                double U_lj = (U[l*n + j] - (R[i*n + j] * Q[l*n + i]));
                U[l*n + j] = U_lj;
            }
        }
    }


    free(U);

    double **res = (double **)calloc(2, sizeof(double *));
    MALLOC_CHECK(res)
    res[0] = Q;
    res[1] = R;
    
    return res;
}

// This function runs the "QR Iteration Algorithm"
// input:   Laplacian matrix ( or any real, symmetric, full rank matrix)
// output:  res[0] = diagonal mat "A-bar" whose diagonal elements approach lap_mat's eigenvalues
//          res[1] = orthogonal mat "Q-bar" whose columns approach lap_mat's eigenvectors
double **QRIteration(const double *A, int n){
    int i, j, converged;
    double **QR_tup;

    //A_bar = A
    double *A_bar =  (double *) calloc(n * n, sizeof(double));
    MALLOC_CHECK(A_bar)
    for (i = 0; i < n * n; i++){
        A_bar[i] = A[i];
    }

    //Q_bar = I
    double *Q_bar =  (double *) calloc(n * n, sizeof(double));
    MALLOC_CHECK(Q_bar)
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++) {
            if (i == j){
                Q_bar[i*n + j] = (double) 1;
            } else {
                Q_bar[i*n + j] = (double) 0;
            }
        }
    }

    double **res = (double **)calloc(2, sizeof(double *));
    MALLOC_CHECK(res)
    res[0] = A_bar;
    res[1] = Q_bar;


    for (i=0; i<n; i++){
        QR_tup = modGrahamSchmidt(res[0], n); //Obtain Q, R for A from the Modified Gram-Schmidt algorithm

        A_bar = multMat(QR_tup[1],n,n,QR_tup[0],n,n); // A_bar = R*Q
        free(res[0]);
        res[0] = A_bar;

        double *QbarQ = multMat(res[1],n,n,QR_tup[0],n,n); // Q_bar*Q

        // done with Q and R from modGrahamSchmidt
        free(QR_tup[0]);
        free(QR_tup[1]);
        free(QR_tup);

        //check if ||Q_bar|-|(Q_bar*Q)|| < epsilon (convergence), if so return.
        converged = 1; //innocent until proven guilty
        for(j=0; j<n*n; j++){
            if (fabs(fabs(res[1][i])-fabs(QbarQ[i]))>EPSILON){ // fabs returns the absolute val of a double
                converged=0;
                break;
            }
        }

        if (converged) {
            free(QbarQ);
            return res;
        } else { // another iteration is required
            free(res[1]);
            res[1] = QbarQ;
        }

    }

    return res;
}


// input: a diagonal matrix nxn "diag_mat" (eigenvalues on the diagonal) 
// output: (int) k according to the eigengap heuristic
int findK(const double *diag_mat, int n) {
    int i, k;

    double *eigenvals = (double *) calloc(n, sizeof(double)); // 1D array of the main diagonal
    MALLOC_CHECK(eigenvals)
    for (i = 0; i < n; ++i) { // filling the vector
        eigenvals[i] = diag_mat[i * n + i];
    }

    qsort(eigenvals, n, sizeof(double), comparefunc); // sort the eigenvalues in ascending order

    double *deltas = (double *) calloc((n-1), sizeof(double));
    MALLOC_CHECK(deltas)

    for (i=0; i<(n-1); ++i) {
        deltas[i] = fabs(eigenvals[i] - eigenvals[i+1]);
    }

    k=0;
    double max_delta_i = deltas[0];
    for (i=1; i<(n/2); ++i) {
        if (deltas[i]>max_delta_i){
            k = i;
        }
    }

    free(eigenvals);
    free(deltas);

    return (k+1);
}

// input:   nxk matrix
// output:  NULL - mutates all rows to be of unit length
void normalizeRows(double *mat, int n, int k) {
    int i, j;

    for (i=0; i<n; i++){
        double rowNorm = 0;
        for (j=0; j<k; j++){ // obtaining the norm of the row
            rowNorm += (mat[i*k + j]*mat[i*k + j]);
        }

        rowNorm = sqrt(rowNorm); // l2 norm

        if (rowNorm != 0){ // rowNorm==0 iff every coordinate is 0
            for (j=0; j<k; j++){ // deviding each cordinate by the norm
                mat[i*k + j] = mat[i*k + j]/rowNorm;
            }
        }
    }

}

// input:    k = the number of columns in the result
//          values = a nxn matrix where the main diagonal is eiganvalues
//          vectors =  a nxn matrix where column i is an eigenvector corresponding to eigenvalue values[i][i]
// output:   nxk matrix where the k columns are the "first" (in sorted order acc. to values) eigenvectors
double *obtainFirstKEigen(int k, double *values, double *vectors, int n){
    int i=0, j=0;

    // build an array of tuples (tup) s.t. eigenvals_tup[i] = tup s.t. tup[0] = values[i][i], tup[1] = i
    double **eigenvals_tup = (double **) calloc(n, sizeof(double*));
    MALLOC_CHECK(eigenvals_tup)

    for (i=0; i<n; ++i) { // filling the vector
        double *tup_i = (double *) calloc(2, sizeof(double));
        MALLOC_CHECK(tup_i)
        tup_i[0] = values[i*n + i];
        tup_i[1] = (double) i;
        eigenvals_tup[i] = tup_i;
    }

    qsort(eigenvals_tup, n, sizeof(double *), tup_comparefunc); // sort the eigenvalues in ascending order

    double *res =  (double *) calloc(k*n, sizeof(double));
    MALLOC_CHECK(res)

    // go through the sorted arr of tups in ascending order. for each tup i, column i of res is the same as vectors
    for (i=0; i<k; i++){
        int col = (int) eigenvals_tup[i][1];
        for (j=0; j<n; j++){
            res[j*k + i] = vectors[j*n + col];
        }
        free(eigenvals_tup[i]);
    }

    for (i=k; i<n; i++){
        free(eigenvals_tup[i]);
    }

    free(eigenvals_tup);

    return res;
}

// coparison func for the fucntion obtainFirstKEigen.
// sorts an array of tuples according to the first value of each tuple.
int tup_comparefunc (const void * a, const void * b){
    double *A_pointer = *((double **) a);
    double *B_pointer = *((double **) b);


    double res = A_pointer[0]-B_pointer[0];
    if (fabs(res) < EPSILON){
        return 0;
    } else if (res<0){
        return -1;
    } else {
        return 1;
    }
}




// this function implements the Spectral Clustering Algorithm
// input:   a set of points as a 1D matrix (nxd) and a k from the user (=0 if k is not provided)
// output:  a tuple Tk_tup s.t. Tk_tup[0] = "T" (nxk matrix from step 5 of algorithm 3)
//          Tk_tup[1] = k from step 3 of the algrithm
double **specClust(double *points, int n, int d, int k_from_user){
    int k;

    double *adj_mat = formAdjMat(points, d, n); //form adj matrix 

    double *diagdeg_mat = formDiagDegreeMat(adj_mat, n, 1); 

    double *lap_mat = formLapMat(adj_mat, diagdeg_mat, n); //form laplacian matrix 

    double **AQ_tuple = QRIteration(lap_mat, n); // AQ_tuple[0] = A_hat (diagonal)  ,AQ_tuple[1] = Q_hat (orthogonal)

    if (k_from_user == 0){
       k = findK(AQ_tuple[0], n); 
    } else {
        k = k_from_user;
    }

    // obtain the first
    double *U = obtainFirstKEigen(k, AQ_tuple[0], AQ_tuple[1], n);
    
    normalizeRows(U, n, k);

    double **Tk_tup = (double **)calloc(2, sizeof(double *));
    Tk_tup[0] = U;

    // TK_tup[1] = k 
    double *k_pointer = malloc(sizeof(double));
    MALLOC_CHECK(k_pointer)
    k_pointer[0] = (double) k;
    Tk_tup[1] = k_pointer;


    free(adj_mat);
    free(diagdeg_mat);
    free(lap_mat);
    free(AQ_tuple[0]);
    free(AQ_tuple[1]);
    free(AQ_tuple);

    return Tk_tup;
}