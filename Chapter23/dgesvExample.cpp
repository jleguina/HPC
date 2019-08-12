#include <cstdlib>
#include <iostream>
using namespace std;

#define F77NAME(x) x##_
extern "C" {
// LAPACK routine for solving systems of linear equations
void F77NAME(dgesv)(const int& n,
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
    const int n = 10; // Define system size
    const int nrhs = 1;
    int info = 0;
    double* A = new double[n * n]; // Allocate matrix and vectors
    double* x = new double[n];
    double* y = new double[n];
    int* ipiv = new int[n]; // Vector for pivots

    // Populate A and x with random numbers
    srand(time(0));
    for(int i = 0; i < n * n; ++i)
	A[i] = (double)(rand()) / RAND_MAX;
    for(int i = 0; i < n; ++i)
	x[i] = (double)(rand()) / RAND_MAX;

    // Compute product y
    F77NAME(dgemv)('N', n, n, 1.0, A, n, x, 1, 0.0, y, 1);

    // Solve system in-place (y is RHS on input and solution on output)
    F77NAME(dgesv)(n, nrhs, A, n, ipiv, y, n, info);

    // Print out result and compare with original input
    for(int i = 0; i < n; ++i){
	cout << x[i] << ", " << y[i] << endl;
	}
	
	
	
}