/*
 * High-performance Computing
 *
 * Using BLAS Level 1 routines
 */
#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;

// Directly calling Fortran API equally valid
#include "cblas.h"

// Generate Symmetric matrix
void ransymMat(int n, double* A)
{
    for(int i = 0; i < n; i++) {      // Rows & Columns of C
	for(int j = 0; j <= i; j++) { // Rows of C
	    A[i * n + j] = (2 * (double)rand() / RAND_MAX) - 1;
	    A[j * n + i] = A[i * n + j];
	}
    }
}

void ranVec(int n, double* x)
{
    for(int i = 0; i < n; i++) {
	x[i] = (2 * (double)rand() / RAND_MAX) - 1;
    }
}

// Call general matrix-matrix multiplication BLAS routine
void multiplyBLAS(const int n, double* A, double* B, double* C)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, B, n, 0.0, C, n);
}

// Call general traspose matrix-matrix multiplication BLAS routine
void transMultiplyBLAS(const int n, double* A, double* B, double* C)
{
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1.0, A, n, B, n, 0.0, C, n);
}

// Call general matrix-vector multiplication routine
void matvecMultiply(const int n, double* A, double* x, double* y)
{
    for(int i = 0; i < n; ++i) {
	for(int j = 0; j < n; ++j) {
	    y[i] += A[i * n + j] * x[j];
	}
    }
}

// Display array function
void displayArray(const int r, const int c, const double* A)
{
    for(int i = 0; i < r; ++i) {
	for(int j = 0; j < c; ++j) {

	    if(j == 0) {
		cout << "[";
	    }
	    cout << A[i * c + j];

	    if(j >= 0 && j < (c - 1)) {
		cout << "  ";
	    }

	    if(j == (c - 1)) {
		cout << "]" << endl;
	    }
	}
    }
}

// Display vector function
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

// Transpose of a matrix
void transpose(const int r, const int c, double* A, double* At)
{

    // cout << "The matrix is:" << endl;
    // displayArray(r, c, A);

    for(int i = 0; i < r; ++i)
	for(int j = 0; j < c; ++j) {
	    At[j * r + i] = A[i * r + j];
	}
    // cout << "The transpose of the matrix is:" << endl;
    // displayArray(r, c, At);
}

// Conjugate gradient algorithm using BLAS routines
void gradientAlgorithm(const int n, double* A, double* b, double* x, double* x2)
{
    double* r = new double[n];
    double* r2 = new double[n];
    double* rr = new double[n];
    double* Ap = new double[n];
    double* p = new double[n];
    double* pt = new double[n];
    double* ptA = new double[n];
    double* ptAp = new double[n];
    double* rt = new double[n];
    double* rtr = new double[n];
    double* r2t = new double[n];
    double* r2tr2 = new double[n];
    double* p2 = new double[n];
    double* x_0 = new double[n];
    double l, alpha, beta;
    double eps = 0.0001;
    double mod = 10000;

    ranVec(n, x_0);

    // displayVector(n, x_0);

    matvecMultiply(n, A, x_0, r);

    for(int i = 0; i < n; ++i) {
	r[i] = b[i] - r[i];
    }

    // cout << "Vector r is:" << endl;
    // displayVector(n, r);

    p = r;
    int k = 0;

    while(mod > eps) {

	transpose(n, 1, p, pt);

	// cout << "Vector r is:" << endl;
	// displayVector(n, r);

	multiplyBLAS(n, pt, A, ptA); // These two are wrong, need to find vec-mat multiplication
	multiplyBLAS(n, ptA, p, ptAp);

	transpose(n, 1, r, rt);
	multiplyBLAS(n, rt, r, rtr); // Wrong, need to find vec-mat multiplication

	alpha = rtr[1] / ptAp[1];

	cout << "Alpha is:" << endl;
	cout << alpha << endl;

	multiplyBLAS(n, A, p, Ap);

	for(int i = 0; i < n; ++i) {
	    x2[i] = x[i] + alpha * p[i];
	    r2[i] = r[i] - alpha * Ap[i];
	}

	for(int i = 0; i < n; ++i) {
	    l += r2[i] * r2[i];
	}

	mod = sqrt(l);

	transpose(n, 1, r2, r2t);
	multiplyBLAS(n, r2t, r2, r2tr2);

	beta = r2tr2[1] / rtr[1];

	for(int i = 0; i < n; ++i) {
	    p2[i] = r2[i] + beta * p[k];
	}

	p = p2;
	r = r2;
	x = x2;

	cout << "K is:" << k++ << endl;
    }
}

int main()
{
    // srand(time(0));

    // Declare and define number of elemnents in matrix 1
    const int n = 3;

    // Declare arrays
    double* A = new double[n * n];
    // double* B = new double[n * n];
    // double* C = new double[n * n];
    double* M = new double[n * n];
    double* x = new double[n];
    double* x2 = new double[n];
    double* b = new double[n];

    // Generate and display array M
    ransymMat(n, M);
    // cout << "Matrix M is:" << endl;
    // displayArray(n, n, M);

    // Generate and display vector x
    ranVec(n, x);
    // cout << "Vector x is:" << endl;
    // displayVector(n, x);

    // Find positive-definite matrix
    transMultiplyBLAS(n, M, M, A);
    // cout << "Matrix A is:" << endl;
    // displayArray(n, n, A);

    // Matrix-vector product
    matvecMultiply(n, A, x, b);
    // cout << "Vector b is:" << endl;
    // displayVector(n, b);

    gradientAlgorithm(n, A, b, x, x2);

    cout << "Vector x2 is:" << endl;
    displayVector(n, x2);
}
