#include <cstdlib>
#include <iostream>
using namespace std;
#define F77NAME(x) x##_
extern "C" {
double F77NAME(ddot)(const int& n, const double* x, const int& incx, const double* y, const int& incy);
}

int main()
{
    const int N = 64;          // Size of vectors
    double* x = new double[N]; // Allocate vectors
    double* y = new double[N];
    double z = 0.0;

    srand(time(0));                         // Initialise random numbers
    for(int i = 0; i < N; ++i) {            // Populate vectors
	x[i] = (double)(rand()) / RAND_MAX; // with numbers in [0,1]
	y[i] = (double)(rand()) / RAND_MAX;
    }

    z = F77NAME(ddot)(N, x, 1, y, 1); // Compute dot-product using BLAS

    cout << "Dot product result: " << z << endl;

    delete[] x; // Deallocate memory!
    delete[] y;
}