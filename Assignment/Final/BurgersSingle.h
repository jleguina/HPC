/*
 * High-performance Computing
 *
 * Burgers class header file.
 */

#ifndef CLASS_BURGERSSINGLE
#define CLASS_BURGERSSINGLE

#include "BurgersSingle.h"
#include "Model.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h> /* sqrt */
//#include "mpi.h"

class BurgersSingle
{

    friend class Model;

public:
    BurgersSingle(Model m);
    ~BurgersSingle();

    // Initialise functions
    void InitialVelocity(double* U, double* V, int Nx_, int Ny_, double Lx_, double Ly_, double dx_, double dy_);

    void displayArray(double* U, int Nx, int Ny);

    void TimeIntegrate(double* U,
        double* V,
        double* U_temp,
        double* V_temp,
        int Nx_,
        int Ny_,
        double T_,
        double dx_,
        double dy_,
        double dt_,
        double ax_,
        double ay_,
        double b_,
        double c_);

    void writeFile(double* U, double* V, int Nx_, int Ny_);

    double energy(double* U, double* V, int Nx_, int Ny_, double dx_, double dy_);

    double Solve();

private:
    // Model m;

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

    double E;

    double* U = new double[Ny_ * Nx_];
    double* V = new double[Ny_ * Nx_];

    double* U_temp = new double[Ny_ * Nx_];
    double* V_temp = new double[Ny_ * Nx_];

    Model ModMS;
};

#endif