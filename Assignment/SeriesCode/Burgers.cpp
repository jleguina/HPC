/*
 * High-performance Computing
 *
 *
 * Burger class implementaiton.
 */

#include "Burgers.h"
#include "Model.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h> /* sqrt */
//#include "mpi.h"

Burgers::Burgers(Model m)
		:ModM (m) 
{
    // Define constant values / operations here

    // Why not use friendship??
       
    // Must delete all variables after
}

Burgers::~Burgers()
{
    // nothing to do
}

double Burgers::Solve(){
	int Nx_ = ModM.GetNx();
    int Ny_ = ModM.GetNy();
    // int Nt_ = m.GetNt(); Not needed

    double Lx_ = ModM.GetLx();
    double Ly_ = ModM.GetLy();
    double c_ = ModM.GetC();
    double b_ = ModM.GetB();
    double ax_ = ModM.GetAx();
    double ay_ = ModM.GetAy();
    double T_ = ModM.GetT();
    double dx_ = ModM.GetDx();
    double dy_ = ModM.GetDy();
    double dt_ = ModM.GetDt();

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

void Burgers::InitialVelocity(double* U, double* V, int Nx_, int Ny_, double Lx_, double Ly_, double dx_, double dy_)
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

    //displayArray(U, Nx_, Ny_);
    //displayArray(V, Nx_, Ny_);
}

void Burgers::displayArray(double* U, int Nx_, int Ny_)
{

    for(int j = 0; j < Ny_; ++j) {
	for(int i = 0; i < Nx_; ++i) {
	    std::cout << U[Nx_ * j + i] << ",";
	}
	std::cout << std::endl;
    }
}

void Burgers::TimeIntegrate(double* U,
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
	memcpy(&V, &U, sizeof(V));
	
    //displayArray(U, Nx_, Ny_);
    //displayArray(V, Nx_, Ny_);
	
    memcpy(&U_temp, &U, sizeof(U_temp));
    memcpy(&V_temp, &V, sizeof(V_temp));

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
	memcpy(&U, &U_temp, sizeof(U));
	memcpy(&V, &V_temp, sizeof(V));
    }

     displayArray(U, Nx_, Ny_);
     displayArray(V, Nx_, Ny_);
}

void Burgers::writeFile(double* U, double* V, int Nx_, int Ny_)
{

    using namespace std;

    ofstream vOut("velField.txt", ios::out | ios::trunc);
    if(vOut.good()) {
	for(int j = 0; j < Ny_ ; j++) {
	    for(int i = 0; i < Nx_; i++) {

		vOut << "(" << U[Nx_ * j + i] << ", " << V[Nx_ * j + i] << ")"
		     << "  ";
	    }
	    vOut << endl;
	}
    }
    vOut.close();
    return;
}

double Burgers::energy(double* U, double* V, int Nx_, int Ny_, double dx_, double dy_)
{
    double energy1 = 0;
    double* E = new double[Ny_ * Nx_];

    for(int j = 0; j < Ny_; j++) {
	for(int i = 0; i < Nx_; i++) {

	    E[Nx_ * j + i] = pow(U[Nx_ * j + i], 2) + pow(V[Nx_ * j + i], 2);
	}
    }

    for(int j = 0; j < (Ny_ - 1); j++) {
	for(int i = 0; i < (Nx_ - 1); i++) {

	    energy1 =+ 0.25 * dx_ * dy_ *
	        (E[Nx_ * j + i] + E[Nx_ * j + i + 1] + E[Nx_ * (j + 1) + i] + E[Nx_ * (j + 1) + i + 1]);
	}
    }
    energy1 = 0.5 * energy1;
    return energy1;
}