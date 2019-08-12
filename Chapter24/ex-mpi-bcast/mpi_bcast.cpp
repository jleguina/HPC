/*
 * High-performance Computing
 *
 * Broadcast a number to all processes.
 */
#include <iostream>
#include <sstream>
using namespace std;

#include <mpi.h>

int main(int argc, char* argv[]) {

    int n    = 0;
    int rank = 0;

    // Initialise MPI.
    int err = MPI_Init(&argc, &argv);
    if (err != MPI_SUCCESS) {
        cout << "Failed to initialise MPI" << endl;
        return -1;
    }

    // Get the rank and comm size on each process.
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // We ask the user for a number, but only on rank 0.
    if (rank == 0) {
        string input;
        cout << "Enter an integer: ";
        cin >> input;
        stringstream S(input);
        S >> n;
    }

    // Broadcast the number from rank 0 to all other processes.
    // All processes must call this function to avoid deadlock!
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Print out rank and the number on all processes.
    cout << "Rank: " << rank << ", n=" << n << endl;

    // Finalise MPI.
    MPI_Finalize();

    return 0;
}
