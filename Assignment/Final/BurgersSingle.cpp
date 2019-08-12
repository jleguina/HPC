/*
 * High-performance Computing
 *
 *
 * BurgersSingle class implementaiton.
 */

#include "BurgersSingle.h"
#include "Model.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h> /* sqrt */
//#include "mpi.h"

BurgersSingle::BurgersSingle(Model m)
    : ModMS(m)
{
}

BurgersSingle::~BurgersSingle()
{
    // nothing to do
}

double BurgersSingle::Solve()
{
    int Nx_ = ModMS.GetNx();
    int Ny_ = ModMS.GetNy();
    // int Nt_ = m.GetNt(); Not needed

    double Lx_ = ModMS.GetLx();
    double Ly_ = ModMS.GetLy();
    double c_ = ModMS.GetC();
    double b_ = ModMS.GetB();
    double ax_ = ModMS.GetAx();
    double ay_ = ModMS.GetAy();
    double T_ = ModMS.GetT();
    double dx_ = ModMS.GetDx();
    double dy_ = ModMS.GetDy();
    double dt_ = ModMS.GetDt();

    double E;

    double* U = new double[Ny_ * Nx_];
    double* V = new double[Ny_ * Nx_];

    double* U_temp = new double[Ny_ * Nx_];
    double* V_temp = new double[Ny_ * Nx_];

    InitialVelocity(U, V, Nx_, Ny_, Lx_, Ly_, dx_, dy_);

    TimeIntegrate(U, V, U_temp, V_temp, Nx_, Ny_, T_, dx_, dy_, dt_, ax_, ay_, b_, c_);

    writeFile(U, V, Nx_, Ny_);

    E = energy(U, V, Nx_, Ny_, dx_, dy_);

    return E;
}

void BurgersSingle::InitialVelocity(double* U,
    double* V,
    int Nx_,
    int Ny_,
    double Lx_,
    double Ly_,
    double dx_,
    double dy_)
{
    double r;
    double x;
    double y;

    for(int j = 0; j < Ny_; j++) {
	for(int i = 0; i < Nx_; i++) {
	    x = -(Lx_ / 2) + i * dx_;
	    y = (Ly_ / 2) - j * dy_;
	    r = sqrt(x * x + y * y);
	    if(r <= 1) {
		U[Nx_ * j + i] = 2 * pow((1 - r), 4) * (4 * r + 1);
	    }
	    if(r > 1.0) {
		U[Nx_ * j + i] = 0.0;
	    }
	}
    }
    memcpy(&V, &U, sizeof(V));
}

void BurgersSingle::displayArray(double* U, int Nx_, int Ny_)
{

    for(int j = 0; j < Ny_; ++j) {
	for(int i = 0; i < Nx_; ++i) {
	    std::cout << U[Nx_ * j + i] << ",";
	}
	std::cout << std::endl;
    }
}

void BurgersSingle::TimeIntegrate(double* U,
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
    double c_)
{
    // Initial U = Initial V
    memcpy(&V, &U, sizeof(V));

    // Copy velocity vectors to temporal variables
    memcpy(&U_temp, &U, sizeof(U_temp));
    memcpy(&V_temp, &V, sizeof(V_temp));

    // Start Time Loop
    for(double t = 0; t < T_; t = t + dt_) {

	for(int j = 1; j < Ny_ - 1; j++) {
	    for(int i = 1; i < Nx_ - 1; i++) {

		U_temp[Nx_ * j + i] = U[Nx_ * j + i] +
		    dt_ * (c_ * ((U[Nx_ * j + i + 1] - 2 * U[Nx_ * j + i] + U[Nx_ * j + i - 1]) * pow(dx_, -2) +
		                    (U[Nx_ * (j + 1) + i] - 2 * U[Nx_ * j + i] + U[Nx_ * (j - 1) + i]) * pow(dy_, -2)) -
		              (ax_ + b_ * U[Nx_ * j + i]) * (U[Nx_ * j + i] - U[Nx_ * j + i - 1]) * pow(dx_, -1) -
		              (ay_ + b_ * V[Nx_ * j + i]) * (U[Nx_ * j + i] - U[Nx_ * (j - 1) + i]) * pow(dy_, -1));

		V_temp[Nx_ * j + i] = V[Nx_ * j + i] +
		    dt_ * (c_ * ((V[Nx_ * j + i + 1] - 2 * V[Nx_ * j + i] + V[Nx_ * j + i - 1]) * pow(dx_, -2) +
		                    (V[Nx_ * (j + 1) + i] - 2 * V[Nx_ * j + i] + V[Nx_ * (j - 1) + i]) * pow(dy_, -2)) -
		              (ax_ + b_ * U[Nx_ * j + i]) * (V[Nx_ * j + i] - V[Nx_ * j + i - 1]) * pow(dx_, -1) -
		              (ay_ + b_ * V[Nx_ * j + i]) * (V[Nx_ * j + i] - V[Nx_ * (j - 1) + i]) * pow(dy_, -1));
	    }
	}
	// Trnasfer temporal variables to U & V
	memcpy(&U, &U_temp, sizeof(U));
	memcpy(&V, &V_temp, sizeof(V));
    }

    // Display
    displayArray(U, Nx_, Ny_);
    displayArray(V, Nx_, Ny_);
}

void BurgersSingle::writeFile(double* U, double* V, int Nx_, int Ny_)
{

    using namespace std;

    // Open file and check it
    ofstream vOut("velField.txt", ios::out | ios::trunc);
    if(vOut.good()) {
	for(int j = 0; j < Ny_; j++) {
	    for(int i = 0; i < Nx_; i++) {

		// Print
		vOut << "(" << U[Nx_ * j + i] << ", " << V[Nx_ * j + i] << ")"
		     << "  ";
	    }
	    vOut << endl;
	}
    }
    // Close file
    vOut.close();
    return;
}

double BurgersSingle::energy(double* U, double* V, int Nx_, int Ny_, double dx_, double dy_)
{
    double energy1 = 0;
    double* E = new double[Ny_ * Nx_];

    // Calculate modulus
    for(int j = 0; j < Ny_; j++) {
	for(int i = 0; i < Nx_; i++) {

	    E[Nx_ * j + i] = pow(U[Nx_ * j + i], 2) + pow(V[Nx_ * j + i], 2);
	}
    }

    // Integrate modulus over the range by means of a trapezoidal integration scheme
    for(int j = 0; j < (Ny_ - 1); j++) {
	for(int i = 0; i < (Nx_ - 1); i++) {

	    energy1 = +0.25 * dx_ * dy_ *
	        (E[Nx_ * j + i] + E[Nx_ * j + i + 1] + E[Nx_ * (j + 1) + i] + E[Nx_ * (j + 1) + i + 1]);
	}
    }
    energy1 = 0.5 * energy1;
    return energy1;
}