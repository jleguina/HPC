// C++ program to find adjoint and inverse of a matrix
#include <bits/stdc++.h>
using namespace std;

// Function to get cofactor of A[n][n] in temp[][]. n is current
// dimension of A[][]
void getCofactor(const int n, double* A, double* temp)
{
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for(int row = 0; row < n; row++) {
	for(int col = 0; col < n; col++) {
	    //  Copying into temporary matrix only those element
	    //  which are not in given row and column
	    if(row != n && col != n) {
		j++;
		temp[i * n + j] = A[row * col];

		// Row is filled, so increase row index and
		// reset col index
		if(j == n - 1) {
		    j = 0;
		    i++;
		}
	    }
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

/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
double determinant(double* A, const int n)
{
    int D = 0; // Initialize result

    //  Base case : if matrix contains single element
    if(n == 1)
	return A[0 * 0];

    double* temp[n * n]; // To store cofactors

    int sign = 1; // To store sign multiplier

    // Iterate for each element of first row
    for(int f = 0; f < n; f++) {
	// Getting Cofactor of A[0][f]
	getCofactor(A, temp, n);
	D += sign * A[0 * f] * determinant(temp, n - 1);

	// terms are to be added with alternate sign
	sign = -sign;
    }

    return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(const int n, double* A, double* adj)
{
    if(n == 1) {
	adj[0 * 0] = 1;
	return;
    }

    // temp is used to store cofactors of A[][]
    int sign = 1, temp[n * n];

    for(int i = 0; i < n; i++) {
	for(int j = 0; j < n; j++) {
	    // Get cofactor of A[i][j]
	    getCofactor(A, temp, i, j, n);

	    // sign of adj[j][i] positive if sum of row
	    // and column indexes is even.
	    sign = ((i + j) % 2 == 0) ? 1 : -1;

	    // Interchanging rows and columns to get the
	    // transpose of the cofactor matrix
	    adj[n * j + i] = (sign) * (determinant(temp, n - 1));
	}
    }
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool inverse(int A[n * n], float inverse[n * n])
{
    // Find determinant of A[][]
    int det = determinant(A, N);
    if(det == 0) {
	cout << "Singular matrix, can't find its inverse";
	return false;
    }

    // Find adjoint
    int adj[n * n];
    adjoint(A, adj);

    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for(int i = 0; i < n; i++)
	for(int j = 0; j < n; j++)
	    inverse[i * n + j] = adj[i * n + j] / float(det);

    return true;
}

// Generic function to display the matrix.  We use it to display
// both adjoin and inverse. adjoin is integer matrix and inverse
// is a float.
template <class T> void display(T A[n * n])
{
    for(int i = 0; i < n; i++) {
	for(int j = 0; j < n; j++)
	    cout << A[n * i + j] << " ";
	cout << endl;
    }
}

// Driver program
int main()
{
    const int n = 4;
    double* A = new double[n];
    double* adj = new double[n]; // To store adjoint of A[][]
    double* inv = new double[n]; // To store inverse of A[][]

    for(int i = 1; i < n : i++) {
	for(int j = 1; j < n : j++) {
	    A[n * i + j] = rand();
	}
    }

    cout << "Input matrix is :\n";
    displayArray(n, n, A)

            cout
        << "\nThe Adjoint is :\n";
    adjoint(A, adj);
    display(adj);

    cout << "\nThe Inverse is :\n";
    if(inverse(A, inv))
	display(inv);

    return 0;
}