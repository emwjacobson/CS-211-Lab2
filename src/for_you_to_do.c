#include "../include/for_you_to_do.h"

int get_block_size(){
    //return the block size you'd like to use
    // Should be a multiple of 3, as the DGEMM uses blockes of 3x3
    return 3*10;

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
        double max = fabs(A[i*n + i]);
        // Find greatest value in column
        for(t=i+1; t<n; t++) {
            if (fabs(A[t*n + i]) > max) {
                maxind = t;
                max = fabs(A[t*n + i]);
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
            for(j=0; j<n; j++) {
                double tempv = A[i*n + j];
                A[i*n + j] = A[maxind*n + j];
                A[maxind*n + j] = tempv;
            }
        }
        // Factorization
        for(j=i+1; j<n; j++) {
            A[j*n + i] = A[j*n + i] / A[i*n + i];
            for (k=i+1; k<n; k++) {
                A[j*n + k] = A[j*n + k] - (A[j*n + i] * A[i*n + k]);
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

    if (UPLO == 'L') { // Forward
        y[0] = B[ipiv[0]];
        for(i=1; i<n; i++) {
            sum = 0;
            for(j=0; j<i; j++) {
                sum += y[j] * A[i*n + j];
            }
            y[i] = B[ipiv[i]] - sum;
        }
    } else if (UPLO == 'U') { // Backward
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
 * Block size must be multiple of 3
 **/
void mydgemm(double *A, double *B, double *C, int n, int in_i, int in_j, int in_k, int b)
{
    int i, j, k;
    int i1, j1, k1;

    for(i=0; i<in_i; i+=b) {
        for(j=0; j<in_j; j+=b) {
            for(k=0; k<in_k; k+=b) {
                // Cache reuse block (multiple of 3x3)
                for (i1 = i; i1 < i+b; i1+=3) {
                    for (j1 = j; j1 < j+b; j1+=3) {
                        int ij = i1*in_i+j1; int ijn = ij+in_i; int ijnn = ijn+in_i;

                        register double c00 = C[ij]; register double c01 = C[ij+1]; register double c02 = C[ij+2];
                        register double c10 = C[ijn]; register double c11 = C[ijn+1]; register double c12 = C[ijn+2];
                        register double c20 = C[ijnn]; register double c21 = C[ijnn+1]; register double c22 = C[ijnn+2];
                        for (k1 = k; k1 < k+b; k1+=3) {
                            register int ik = i1*n+k1; register int ikn = ik+n; register int iknn = ikn+n;
                            register int kj = k1*n+j1; register int kjn = kj+n; register int kjnn = kjn+n;

                            register double a00 = A[ik]; register double a10 = A[ikn]; register double a20 = A[iknn];
                            register double b00 = B[kj]; register double b01 = B[kj+1]; register double b02 = B[kj+2];

                            c00 += a00*b00;
                            c01 += a00*b01;
                            c02 += a00*b02;
                            c10 += a10*b00;
                            c11 += a10*b01;
                            c12 += a10*b02;
                            c20 += a20*b00;
                            c21 += a20*b01;
                            c22 += a20*b02;

                            a00 = A[ik+1]; a10 = A[ikn+1]; a20 = A[iknn+1];
                            b00 = B[kjn]; b01 = B[kjn+1]; b02 = B[kjn+2];

                            c00 += a00*b00;
                            c01 += a00*b01;
                            c02 += a00*b02;
                            c10 += a10*b00;
                            c11 += a10*b01;
                            c12 += a10*b02;
                            c20 += a20*b00;
                            c21 += a20*b01;
                            c22 += a20*b02;

                            a00 = A[ik+2]; a10 = A[ikn+2]; a20 = A[iknn+2];
                            b00 = B[kjnn]; b01 = B[kjnn+1]; b02 = B[kjnn+2];

                            c00 += a00*b00;
                            c01 += a00*b01;
                            c02 += a00*b02;
                            c10 += a10*b00;
                            c11 += a10*b01;
                            c12 += a10*b02;
                            c20 += a20*b00;
                            c21 += a20*b01;
                            c22 += a20*b02;
                        }
                        C[ij] = c00; C[ij+1] = c01; C[ij+2] = c02;
                        C[ijn] = c10; C[ijn+1] = c11; C[ijn+2] = c12;
                        C[ijnn] = c20; C[ijnn+1] = c21; C[ijnn+2] = c22;
                    }
                }
            }
        }
    }
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
    int ib, i, t, j, k;
    for(ib=0; ib<n-1; ib+=b) {
        int end = ib+b;

        // (1) BLAS2 version to calculate P*L*U

        // Process each column in block
        for(i=ib; i<end; i++) {
            // First, do pivoting
            int maxind = i;
            double max = fabs(A[i*n + i]);
            // Find greatest value in column
            for(t=i+1; t<n; t++) {
                if (fabs(A[t*n + i]) > max) {
                    maxind = t;
                    max = fabs(A[t*n + i]);
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
                for(j=0; j<n; j++) {
                    double tempv = A[i*n + j];
                    A[i*n + j] = A[maxind*n + j];
                    A[maxind*n + j] = tempv;
                }
            }

            // Second, factorize remaining columns
            for(j=i+1; j<n; j++) {
                register double a = A[j*n + i] = A[j*n + i] / A[i*n + i];
                for (k=i+1; k<end; k++) {
                    A[j*n + k] = A[j*n + k] - (a * A[i*n + k]);
                }
            }
        }

        // (2) Update the rows
        int count=1;
        for(i=ib+1; i<end; i++) { // Row (+1 because we shouldnt have to do anything for the first row)
            for(k=0; k<count; k++) {
                register double temp = A[i*n + (k+ib)];
                for(j=end; j<n; j++) { // Col
                    A[i*n + j] -= (temp * A[(k+ib)*n + j]);
                }
            }
            count++;
        }

        // (3) Update rest of matrix (green area)
        // mydgemm(cont double *A, const double *B, double *C, int n, int i, int j, int k, int b);
        // Need to multiply A[end+1:n, ib:end] by A[ib:end, end+1:n]
        int ii = n-end;
        int kk = end-ib;

        double *C = (double *)calloc(n * n, sizeof(double));
        mydgemm(&A[end*n + ib], &A[ib*n + end], C, n, ii, ii, kk, b);

        for(i=ib+b; i<n; i++) { // Row
            for(j=ib+b; j<n; j++) { // Col
                A[i*n + j] -= C[(i-end)*ii + (j-end)];
            }
        }

        free(C);
    }
    return 0;
}

