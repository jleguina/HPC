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

double L2_mod(const int n, const double* A)
{
    double d = 0;
    double d2 = 0;
    for(int i = 0; i < n; ++i) {
	d2 += pow(A[i], 2);
    }
    d = sqrt(d2);
    return d;
}

double sum_vector(const int n, const double* A)
{
    double sum = 0;
    for(int i = 0; i < n; ++i) {
	sum += A[i];
    }
    return sum;
}

double max_vector(const int n, const double* A)
{
    double max = 0;
    for(int i = 0; i < (n); ++i) {
	if(A[i + 1] > max) {
	    max = A[i + 1];
	}

	if(A[0] > max) {
	    max = A[0];
	}
    }

    return max;
}

int main()
{
    // Declare and define number of elemnents in vector
    int n;
    cout << "Please input number of elements in array: " << endl;
    cin >> n;

    // Declare array
    double* A = new double[n];

    // Generate and display array
    cout << "Array A is: " << endl;
    for(int i = 0; i < n; ++i) {
	A[i] = (double)rand() / RAND_MAX;
	cout << A[i] << endl;
    }

    // Use BLAS routines
    // Use BLAS to generate L2-norm
    double norm = cblas_dnrm2(n, A, 1);

    // Since x[i] >= 0, we can use DASUM to compute the sum
    double sum = cblas_dasum(n, A, 1);

    // Use the BLAS routine IDAMAX to find the index of entry with maximum value
    int imax = cblas_idamax(n, A, 1);
    double max = A[imax];

    // Use custom funcions
    double max_cus = max_vector(n, A);
    double sum_cus = sum_vector(n, A);
    double L2_cus = L2_mod(n, A);

    // Compare BLAS and manual values
    cout << "Maximum value. BLAS: " << max << ". Custom: " << max_cus << endl;
    cout << "Summation value. BLAS: " << sum << ". Custom: " << sum_cus << endl;
    cout << "L2 norm value. BLAS: " << norm << ". Custom: " << L2_cus << endl;
}