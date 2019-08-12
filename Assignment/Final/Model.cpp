/*
 * High-performance Computing
 *
 * Model class implementaiton.
 */

#include <iostream> /* cerr*/
#include <stdlib.h> /* atof */

#include "Model.h"
#include "mpi.h"

Model::Model(int argc, char* argv[])
{

    ParseParameters(argc, argv);
    // PrintParameters();
}

Model::Model()
{
}

Model::~Model()
{
    // nothing to do
}

void Model::ParseParameters(int argc, char* argv[])
{
    x0 = atof(argv[1]);
    y0 = atof(argv[2]);
    Lx = atof(argv[3]);
    Ly = atof(argv[4]);
    T = atof(argv[5]);
    Nx = atof(argv[6]);
    Ny = atof(argv[7]);
    Nt = atof(argv[8]);

    ax = atof(argv[9]);
    ay = atof(argv[10]);
    b = atof(argv[11]);
    c = atof(argv[12]);

    if(argc == 15) {
	Px = atof(argv[13]);
	Py = atof(argv[14]);
    }

    dx = Lx / (Nx - 1);
    dy = Ly / (Ny - 1);
    dt = T / (Nt - 1);

    if(argc != 14 || argc != 12 || !(int)Nx || !(int)Ny || !(int)Nt || dx < 0 || dy < 0 || dt < 0) {
	validity = false;
    } else {
	validity = true;
    }
}

void Model::PrintParameters()
{
    std::cout << "ax = " << ax << std::endl;
    std::cout << "ay = " << ay << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "x0 = " << x0 << std::endl;
    std::cout << "y0 = " << y0 << std::endl;
    std::cout << "Lx = " << Lx << std::endl;
    std::cout << "Ly = " << Ly << std::endl;
    std::cout << "T = " << T << std::endl;
    std::cout << "Nx = " << Nx << std::endl;
    std::cout << "Ny = " << Ny << std::endl;
    std::cout << "Nt = " << Nt << std::endl;
    std::cout << "dx = " << dx << std::endl;
    std::cout << "dy = " << dy << std::endl;
    std::cout << "dt = " << dt << std::endl;
}
