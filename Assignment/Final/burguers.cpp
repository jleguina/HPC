#include <chrono>
#include <iostream>

#include "Burgers.h"
#include "BurgersSingle.h"
#include "Model.h"
#include "mpi.h"

int main(int argc, char* argv[])
{

    Model m(argc, argv);

    Burgers b(m);

    BurgersSingle c(m);

    // Remove velocity file
    // remove("./velField.txt");

    // Call code to initialise the problem here

    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;

    hrc::time_point start = hrc::now();

    int Px = GetPx();
    int Py = GetPy();
    int rank;

    // Initialise MPI.
    int err = MPI_Init(&argc, &argv);
    if(err != MPI_SUCCESS) {
	std::cout << "Failed to initialise MPI" << std::endl;
	return -1;
    }

    if(Px == 1 && Py == 1) {

	double E = c.Solve();

	std::cout << E << std::endl;

    } else {

	// Find rank
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double E = b.Solve();

	if(rank == 0) {
	    std::cout << "E = " << E << std::endl;
	}

	MPI_Finalize();
    }

    hrc::time_point end = hrc::now();
    std::cout << "Time taken to finish process with rank = " << rank
              << " is: " << std::chrono::duration_cast<ms>(end - start).count() << " ms" << std::endl;

    return 0;
}