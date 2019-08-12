#include <cstdlib>
#include <iostream>
#include <cstdlib>
using namespace std;

#include <cblas.h>

int main()
{
    const int N = 8; // Columns
    const int M = 6; // Rows
    double* A = new double[N * M];
    double* x = new double[N];
    double* y = new double[M];

    for(int i = 0; i < N; ++i) {
	for(int j = 0; j < M; ++j) {
	    A[i * M + j] = double(rand()) / RAND_MAX;
	}
	x[i] = double(rand()) / RAND_MAX;
    }

    cblas_dgemv(CblasColMajor, CblasNoTrans, M, N, 1.0, A, M, x, 1, 0.0, y, 1);


for(int i = 0; i < N; ++i) { // Display matrix
	for(int j = 0; j < M; ++j) {
	     cout << "Element " << j+1 << ", " << i+1 << " of matrix A = x*y is: " << A[i * j] << endl;
	}
 }

for(int i = 0; i < N; ++i) { // Display vector x
	     cout << "Element " << i+1 << " of vector x is: " << x[i] << endl;
}

for(int i = 0; i < M; ++i) { // Display vector y
	     cout << "Element " << i+1 << " of vector y is: " << y[i] << endl;
}


}