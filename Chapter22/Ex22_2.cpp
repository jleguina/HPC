/*
 * High-performance Computing
 *
 * Using BLAS Level 1 routines
 */
#include <cmath>
#include <iostream>

using namespace std;

// Directly calling Fortran API equally valid
#include "cblas.h"

// Multiply the two matrices A and B into the matrix C
void multiply(int n, double* A, double* B, double* C)
{
    for(int i = 0; i < n; i++) {     // Columns of C
	for(int j = 0; j < n; j++) { // Rows of C
	    C[i * n + j] = 0.0;
	    for(int k = 0; k < n; k++) { // Along row of A / Column of B
		C[i * n + j] += A[i * n + k] * B[k * n + j];
	    }
	}
    }
}

// Print a matrix
void printMatrix(int n, double* A)
{
    for(int i = 0; i < n; i++) {
	for(int j = 0; j < n; j++) {
	    cout << A[i * n + j] << "  ";
	}
	cout << endl;
    }
}

// Call general matrix-matrix multiplication BLAS routine
void multiplyBLAS(int n, double* A, double* B, double* C)
{
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, B, n, 0.0, C, n);
}

void displayArray(const int n, const double* A)
{
    for(int i = 0; i < n; ++i) {
	for(int j = 0; j < n; ++j) {

	    if(j == 0) {
		cout << "[";
	    }
	    cout << A[i * n + j];

	    if(j >= 0 && j < (n - 1)) {
		cout << "  ";
	    }

	    if(j == (n - 1)) {
		cout << "]" << endl;
	    }
	}
    }
}

int main()
{
    // Declare and define number of elemnents in matrix 1
    // const int n = 3;
    // cout << "Please input number n, for an n*n array (n>0): " << endl;
    // cin >> const int n;

    const int n = 2;

    // Declare arrays
    double* A = new double[n * n];
    double* B = new double[n * n];
    double* C1 = new double[n * n];
    double* C2 = new double[n * n];

    // double A[n][n];
    // double B[n][n];
    // double C2[n][n];
    // double C1[n][n];

    // Generate and display array A
    for(int i = 0; i < n; ++i) {
	for(int j = 0; j < n; ++j) {
	    A[i * n + j] = 10 * (double)rand() / RAND_MAX;
	}
    }

    // Generate array B
    for(int i = 0; i < n; ++i) {
	for(int j = 0; j < n; ++j) {
	    B[i *n + j] = 10 * (double)rand() / RAND_MAX;
	}
    }

    cout << "Array A is: " << endl;
    displayArray(n, A);
    //printMatrix(n, A);
    cout << "Array B is: " << endl;
    displayArray(n, B);
    //printMatrix(n, B);

    multiply(n, A, B, C2);
    cout << "Array C2 is: " << endl;
    displayArray(n, C2);
    
	//printMatrix(n, C2);
    cout << "Array C1 is: " << endl;
    
	// Multiply A and B into D using BLAS
    multiplyBLAS(n, A, B, C1);
    displayArray(n, C1);
    //printMatrix(n, C1);
}
