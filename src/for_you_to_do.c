#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use
    /*add your code here */
    return 8;

}

/**
 *
 * this function computes LU factorization
 * for a square matrix
 *
 * syntax
 *
 *  input :
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally
 *
 **/


int mydgetrf(double *A, int *ipiv, int n)
{
    int i, t, j, k;

    for (i=0; i<n-1; i++) {
        int maxind = i;
        int max = abs(A[i*n + i]);
        // Find greatest value in column
        for(t=i+1; t<n; t++) {
            if (abs(A[t*n + i]) > max) {
                maxind = t;
                max = abs(A[t*n + i]);
            }
        }
        if (max == 0)
            return -1; // LU Factorizations failed, coefficient matrix is singular
        if (maxind != i) {
            // Save pivoting information
            int temps = ipiv[i];
            ipiv[i] = ipiv[maxind];
            ipiv[maxind] = temps;
            // Swap rows
            for(j = 0; j<n; j++) {
                double tempv = A[i*n + j];
                A[i*n + j] = A[maxind*n + j];
                A[maxind*n + j] = tempv;
            }
        }
        // Factorization
        for(j=i+1; j<n; j++) {
            A[j*n + i] = A[j*n + i] / A[i*n + i];
            for (k=i+1; k<n; k++) {
                A[j*n + k] = A[j*n + k] - A[j*n + i] * A[i*n + k];
            }
        }
    }
    return 0;
}

/**
 *
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 *
 * syntax
 *
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 *
 *      A     n by n     , square matrix
 *
 *      B     1 by n     , vector
 *
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 *
 *      n                , length of vector / size of matrix
 *
 *  output :
 *      none
 *
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    int i, j;
    double sum;
    double *y = malloc(n * sizeof(double));

    if (UPLO == 'L') { // Forward?
        y[0] = B[ipiv[0]];
        for(i=1; i<n; i++) {
            sum = 0;
            for(j=0; j<i; j++) {
                sum += y[j] * A[i*n + j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
    } else if (UPLO == 'U') { // Backward?
        y[n-1] = B[n-1] / A[(n-1)*n + (n-1)];
        for(i=n-2; i>=0; i--) {
            sum = 0;
            for(j=i+1; j<n; j++) {
                sum += y[j] * A[i*n + j];
            }
            y[i] = (B[i] - sum) / A[i*n + i];
        }
    }

    for(i=0; i<n; i++) {
        B[i] = y[i];
    }

    free(y);

    return;
}

/**
 *
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 *
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
    return;
}

/**
 *
 * this function computes LU decomposition
 * for a square matrix using block gepp introduced in course
 * lecture .
 *
 * just implement the block algorithm you learned in class.
 *
 * syntax
 *
 *  input :
 *
 *
 *      A     n by n     , square matrix
 *
 *
 *
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 *
 *      n                , length of vector / size of matrix
 *
 *      b                , block size
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally
 *
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b)
{
    return 0;
}

