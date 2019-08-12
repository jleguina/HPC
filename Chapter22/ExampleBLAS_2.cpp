#include <cstdlib>
#include <iostream>
#include <cstdlib>
using namespace std;
#define F77NAME(x) x##_
extern "C" {
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

const int N = 3;               // Number of columns
const int M = 2;               // Number of rows
double* A = new double[N * M]; // Matrix allocation
double* x = new double[N];     // Input/output vectors
double* y = new double[M];

// Print a matrix function
void printMatrix(double* A[][M])
{
    for(int i = 0; i < M; i++) {
	for(int j = 0; j < N; j++) {
	    cout << A[i][j] << "  ";
	}
	cout << endl;
    }
}

int main()
{

    for(int i = 0; i < N; ++i) { // Populate matrix and vector
	for(int j = 0; j < M; ++j) {
	    A[i * M + j] = double(rand()) / RAND_MAX;
	}
	x[i] = double(rand()) / RAND_MAX;
    }

    // Perform matrix-vector product.
    F77NAME(dgemv)('N', M, N, 1.0, A, M, x, 1, 0.0, y, 1);


for(int i = 0; i < N; ++i) { // Display matrix
	for(int j = 0; j < M; ++j) {
	     cout << "Element " << j+1 << ", " << i+1 << " of matrix A = x*y is: " << A[i * j] << endl;
	}
 }

for(int i = 0; i < N; ++i) { // Display vector x
	     cout << "Element " << i+1 << " of vextor x is: " << x[i] << endl;
}

for(int i = 0; i < M; ++i) { // Display vector y
	     cout << "Element " << i+1 << " of vextor y is: " << y[i] << endl;
}
 }