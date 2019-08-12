/*
 * High-performance Computing
 *
 * Conjugate Gradient algorithm.
 */
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <math.h>
using namespace std;
#include "cblas.h"
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

// Print a vector
void displayVector(const int n, const double* x)
{
    for(int i = 0; i < n; ++i) {
	if(i == 0) {
	    cout << "[";
	}
	if(i >= 0 && i < (n - 1)) {
	    cout << x[i] << endl;
	}
	if(i == (n - 1)) {
	    cout << x[i] << "]" << endl;
	}
    }
}

// Print an r * c array
void displayArray(const int r, const int c, const double* A)
{
    for(int i = 0; i < r; ++i) {
	for(int j = 0; j < c; ++j) {

	    if(j == 0) {
		cout << "[";
	    }
	    cout << A[i * r + j];

	    if(j >= 0 && j < (c - 1)) {
		cout << "  ";
	    }

	    if(j == (c - 1)) {
		cout << "]" << endl;
	    }
	}
    }
}

// Calculate forcing function
void ForcingFunction(const int n, const double lambda, const double Dx, double* f)
{
    double ffff;
    double x;
    const long double pi = 3.141592653589793238462643383279502884L;
    double step = Dx;

    for(int i = 0; i < n; ++i) {
	x = i * step;
	// ffff = -(lambda + (pi * pi)) * sin(pi * x);
	ffff = (-lambda + (pi * pi)) * cos(pi * x);
	f[i] = ffff;
    }
}

// Solve the linear system
void SolveConjugateGradient(int n, double* A, double* b, double* x)
{
    double* r = new double[n];
    double* p = new double[n];
    double* t = new double[n]; // temp
    int k;
    double alpha;
    double beta;
    double eps;
    double tol = 0.00001;

    cblas_dcopy(n, b, 1, r, 1);                                             // r_0 = b (i.e. b)
    cblas_dsymv(CblasRowMajor, CblasUpper, n, -1.0, A, n, x, 1, 1.0, r, 1); // r_0 = b - A x_0
    cblas_dcopy(n, r, 1, p, 1);                                             // p_0 = r_0
    k = 0;
    do {
	cblas_dsymv(CblasRowMajor, CblasUpper, n, 1.0, A, n, p, 1, 0.0, t, 1); // t = A p_k
	alpha = cblas_ddot(n, t, 1, p, 1);                                     // alpha = p_k^T A p_k
	alpha = cblas_ddot(n, r, 1, r, 1) / alpha;                             // compute alpha_k
	beta = cblas_ddot(n, r, 1, r, 1);                                      // r_k^T r_k

	cblas_daxpy(n, alpha, p, 1, x, 1);  // x_{k+1} = x_k + alpha_k p_k
	cblas_daxpy(n, -alpha, t, 1, r, 1); // r_{k+1} = r_k - alpha_k A p_k

	eps = cblas_dnrm2(n, r, 1);
	cout << "Iteration " << k << ": eps=" << eps << endl;
	if(eps < tol * tol) {
	    break;
	}
	beta = cblas_ddot(n, r, 1, r, 1) / beta;

	cblas_dcopy(n, r, 1, t, 1);
	cblas_daxpy(n, beta, p, 1, t, 1);
	cblas_dcopy(n, t, 1, p, 1);

	k++;
    } while(k < 5000); // Set a maximum number of iterations

    delete[] r;
    delete[] p;
    delete[] t;
}

int main()
{

    const int n = 5;
    const double a = 0.;
    // const double b = 2.;
    const double b = 3.;
    const double lambda = 1.;
    double Dx = (b - a) / (n - 1);
    double alpha = (-2 / (Dx * Dx)) - lambda;
    double beta = 1 / (Dx * Dx);

    double* M = new double[(n - 2) * (n - 2)];
    double* y = new double[(n - 2)];
    double* f = new double[n];
    double* u = new double[n];

    for(int i = 0; i < (n - 2); ++i) {     // index rows
	for(int j = 0; j < (n - 2); ++j) { // index cols
	    if(i == j) {
		M[i * (n - 2) + j] = alpha;
	    }
	    if(i == j + 1) {
		M[i * (n - 2) + j] = beta;
	    }
	    if(i == j - 1) {
		M[i * (n - 2) + j] = beta;
	    }
	}
    }

    cout << alpha << endl;
    cout << beta << endl;
    cout << Dx << endl;
    displayArray(n - 2, n - 2, M);

    ForcingFunction(n, lambda, Dx, f);
    cout << "Vector F is:" << endl;
    displayVector(n, f);

    int k = 1;
    for(int i = 0; i < (n - 2); ++i) { // index rows
	if(i == 0) {
	    y[i] = f[k] - (beta * 0); // Since u[0] = 1
	}
	if(i == (n - 3)) {
	    y[i] = f[k] - (beta * (0)); // Since u[3] = -1
	}
	if(i != 0 && i != (n - 3)) {
	    y[i] = f[k];
	}
	k++;
    }

    cout << "Vector Y is:" << endl;
    displayVector((n - 2), y);

    // Solve
    // SolveConjugateGradient((n - 2), M, y, u);

    const int nrhs = 1;
    int info = 0;
    int* ipiv = new int[n];
	
   //F77NAME(dgemv)('N', n, n, 1.0, M, n, u, 1, 0.0, y, 1);
   //F77NAME(dgesv)(n, nrhs, M, n, ipiv, y, n, info);

    cout << "Vector U is:" << endl;
    displayVector((n - 2), u);

    delete[] M;
    delete[] y;
    delete[] f;
    delete[] u;

    return 0;
}
