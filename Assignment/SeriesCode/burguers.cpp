#include <chrono>
#include <iostream>

#include "Burgers.h"
#include "Model.h"

int main(int argc, char* argv[])
{

    Model m(argc, argv);

    Burgers b(m);


    // Call code to initialise the problem here

    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;

    hrc::time_point start = hrc::now();

    double E = b.Solve();

    std::cout << E << std::endl;

    hrc::time_point end = hrc::now();
	
	std::cout << "Time taken = " << std::chrono::duration_cast<ms>(end-start).count() << std::endl;

    return 0;
}