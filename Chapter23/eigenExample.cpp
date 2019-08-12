#include <cstdlib>
#include <iostream>
using namespace std;

#define F77NAME(x) x##_
extern "C" {
void F77NAME(dsyev)(const char& v,
    const char& ul,
    const int& n,
    double* a,
    const int& lda,
    double* w,
    double* work,
    const int& lwork,
    int* info);
}
int main()
{
    const int n = 10, lda = n, ldv = n;
    int info = 0, lwork = -1;
    double wkopt;
    double* A = new double[n * n];
    double* w = new double[n];
    double* work;

    srand(time(0));
    for(int i = 0; i < n; ++i) // Populate random symmetric matrix
	for(int j = i; j < n; ++j)
	    A[i * n + j] = A[j * n + i] = (double)(rand()) / RAND_MAX;

    // Query for optimal workspace size
    F77NAME(dsyev)('V', 'U', n, A, lda, w, &wkopt, lwork, &info);

    // Allocate workspace
    lwork = (int)wkopt;
    work = new double[lwork];

    // Compute eigenvalues - A is replaced with eigenvectors
    F77NAME(dsyev)('V', 'U', n, A, lda, w, work, lwork, &info);

    // Print eigenvalues.
    for(int i = 0; i < n; ++i) {
	cout << "EV: " << i << ": " << w[i] << endl;
    }

    delete[] A, w, work;
}