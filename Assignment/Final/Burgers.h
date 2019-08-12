/*
 * High-performance Computing
 *
 * Burgers class header file.
 */

#ifndef CLASS_BURGERS
#define CLASS_BURGERS

#include "Burgers.h"
#include "Model.h"
#include "mpi.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h> /* sqrt */

class Burgers
{

    friend class Model;

public:
    Burgers(Model m);
    ~Burgers();

    // Initialise functions
    void displayArray(double* U, int Nx, int Ny);

    void TimeIntegrate(int Nx_,
        int Ny_,
        int Lx_,
        int Ly_,
        double T_,
        double dx_,
        double dy_,
        double dt_,
        double ax_,
        double ay_,
        double b_,
        double c_,
        int Px_,
        int Py_);

    void writeFile(int Nx_, int Ny_, int Px_, int Py_);

    double energy(int Nx_, int Ny_, double dx_, double dy_, int Px_, int Py_);

    double Solve();

private:
    double Ly_;
    double Lx_;
    double T_;
    int Nx_;
    int Ny_;
    int Nt_;
    double dx_;
    double dy_;
    double dt_;

    double ax_;
    double ay_;
    double b_;
    double c_;

    int Px_;
    int Py_;

    double E;

    double* U;
    double* V;

    Model ModM;
};

#endif