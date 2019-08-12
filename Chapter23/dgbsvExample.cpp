#include <cstdlib>
#include <iostream>
using namespace std;

#define F77NAME(x) x##_
extern "C" {
// LAPACK routine for solving systems of linear equations
void F77NAME(dgbsv)(const int& n,
    const int& nrhs,
    const double* A,
    const int& lda,
    int* ipiv,
    double* B,
    const int& ldb,
    int& info);
// BLAS routine for performing matrix-vector multiplication
void F77NAME(dgemv)(const char& trans,
    const int& m,
    const int& n,
    const double& alpha,
    const double* A,
    const int& lda,
    const double* x,
    const int& incx,
    const double& beta,
    double* y,
    const int& incy);
}

int main()
{

    int n = 50;                 // Problem size
    int kl = 1;                 // Lower diagonal bandwidth
    int ku = 2;                 // Upper diagonal bandwidth
    int nrhs = 1;               // Number of RHS vectors
    int ldab = 1 + 2 * kl + ku; // Number of rows in compressed matrix
    int ldb = n;                // Size of RHS vector
    int info;
    double* A = new double[ldab * n]; // Banded matrix + ’fill’ row
    int* piv = new int[n];            // Pivot data
    double* b = new double[n];        // RHS vector / output vector

    // ... Populate matrix A in banded format and vector b here

    // Solve banded matrix system Ax=b.
    // On input, b contains RHS vector. On output, b contains solution.
    F77NAME(dgbsv)(n, kl, ku, nrhs, A, ldab, piv, b, ldb, info);
    if(info) {
	cout << "ERROR: An error occurred in DGBSV: " << info << endl;
    }
}