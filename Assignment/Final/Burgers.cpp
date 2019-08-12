/*
 * High-performance Computing
 *
 *
 * Burger class implementaiton.
 */

#include "Burgers.h"
#include "Model.h"
#include "mpi.h"
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h> /* sqrt */

Burgers::Burgers(Model m)
    : ModM(m)
{
}

Burgers::~Burgers()
{
}

double Burgers::Solve()
{
    int Nx_ = ModM.GetNx();
    int Ny_ = ModM.GetNy();
    int Lx_ = ModM.GetLx();
    int Ly_ = ModM.GetLx();

    double c_ = ModM.GetC();
    double b_ = ModM.GetB();
    double ax_ = ModM.GetAx();
    double ay_ = ModM.GetAy();
    double T_ = ModM.GetT();
    double dx_ = ModM.GetDx();
    double dy_ = ModM.GetDy();
    double dt_ = ModM.GetDt();

    int Px_ = ModM.GetPx();
    int Py_ = ModM.GetPy();

    double E;

    int rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    TimeIntegrate(Nx_, Ny_, Lx_, Ly_, T_, dx_, dy_, dt_, ax_, ay_, b_, c_, Px_, Py_);

    E = energy(Nx_, Ny_, dx_, dy_, Px_, Py_);

    writeFile(Nx_, Ny_, Px_, Py_);

    if(rank == 0) {
	return E;
    }
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

void Burgers::TimeIntegrate(int Nx_,
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
    int Py_)
{

    // Declare variables
    int rank;
    int ini_x, ini_y, fin_x, fin_y, size_y, size_x; // Sizes and bounds
    int i, j, k, l;                                 // Counters
    double r, x, y;                                 // Radius and coordinates
    int copy_place;
    double t;

    memcpy(&V[0], &U[0], (size_x * size_y) * sizeof(double));

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Submatrix row position
    int row = floor(rank / Px_);
    // Submatrix column position
    int col = rank - Px_ * row;

    // Define bounds for each submatrix in region I - first to penultimate rows & columns
    if(row < (Py_ - 1) && col < (Px_ - 1)) {

	ini_x = int(ceil((1.0 * Nx_) / (1.0 * Px_)) * col);
	fin_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * (col + 1));
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * (row + 1));
    }
    // Define bounds for each submatrix in region II - first to penultimate rows & last column
    if(row < (Py_ - 1) && col == (Px_ - 1)) {

	ini_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * col);
	fin_x = int(Nx_);
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * (row + 1));
    }
    // Define bounds for each submatrix in region III - last row & first to penultimate columns
    if(row == (Py_ - 1) && col < (Px_ - 1)) {
	ini_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * col);
	fin_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * (col + 1));
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(Ny_);
    }
    // Define bounds for each submatrix in region IV - last row & last column
    if(row == (Py_ - 1) && col == (Px_ - 1)) {

	ini_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * col);
	fin_x = int(Nx_);
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(Ny_);
    }

    // Determine size of augmented submatrix
    if(row == 0) {
	size_y = fin_y - ini_y + 1;
    }
    if(col == 0) {
	size_x = fin_x - ini_x + 1;
    }
    if(col > 0 && col < Px_ - 1) {
	size_x = fin_x - ini_x + 2;
    }
    if(row > 0 && row < Py_ - 1) {
	size_y = fin_y - ini_y + 2;
    }
    if(row == Py_ - 1) {
	size_y = fin_y - ini_y + 1;
    }
    if(col == Px_ - 1) {
	size_x = fin_x - ini_x + 1;
    }

    // Generate augmented submatrix
    U = new double[size_x * size_y];
    V = new double[size_x * size_y];

    // Define augmented temporal variables
    double* U_temp = new double[size_x * size_y];
    double* V_temp = new double[size_x * size_y];

    k = 0;
    l = 0;
    if(row == 0) {
	j = ini_y;
    }
    if(col == 0) {
	i = ini_x;
    }
    if(col > 0 && col < Px_ - 1) {
	i = ini_x - 1;
    }
    if(row > 0 && row < Py_ - 1) {
	j = ini_y - 1;
    }
    if(row == Py_ - 1) {
	j = ini_y - 1;
    }
    if(col == Px_ - 1) {
	i = ini_x - 1;
    }

    // Define initial velocity field
    for(int k1 = 0; k1 < size_y; k1++) {
	for(int k2 = 0; k2 < size_x; k2++) {
	    x = -(Lx_ / 2) + i * dx_;
	    y = (Ly_ / 2) - j * dy_;
	    r = sqrt(x * x + y * y);
	    if(r <= 1) {
		U[size_x * k1 + k2] = 2 * pow((1 - r), 4) * (4 * r + 1);
	    }
	    if(r > 1.0) {
		U[size_x * k1 + k2] = 0.0;
	    }
	    i++;
	}

	if(col == 0) {
	    i = ini_x;
	}
	if(col > 0 && col < Px_ - 1) {
	    i = ini_x - 1;
	}
	if(col == Px_ - 1) {
	    i = ini_x - 1;
	}

	j++;
    }

    // Same initial velocity in U and V
    memcpy(&V[0], &U[0], (size_x * size_y) * sizeof(double));

    t = 0;

    // Start time loop
    while(t < T_) {

	memcpy(&U_temp[0], &U[0], (size_x * size_y) * sizeof(double));
	memcpy(&V_temp[0], &V[0], (size_x * size_y) * sizeof(double));

	// Integrate from first to last positions of the augmented submatrix
	// Note: this is valid since first and last row/column of augmented matrices are either 0 (B.C.) or part of the
	// augmentation
	for(j = 1; j < size_y - 1; j++) {
	    for(i = 1; i < size_x - 1; i++) {

		U_temp[size_x * j + i] = U[size_x * j + i] +
		    dt_ *
		        (c_ * ((U[size_x * j + i + 1] - 2 * U[size_x * j + i] + U[size_x * j + i - 1]) * pow(dx_, -2) +
		                  (U[size_x * (j + 1) + i] - 2 * U[size_x * j + i] + U[size_x * (j - 1) + i]) *
		                      pow(dy_, -2)) -
		            (ax_ + b_ * U[size_x * j + i]) * (U[size_x * j + i] - U[size_x * j + i - 1]) *
		                pow(dx_, -1) -
		            (ay_ + b_ * V[size_x * j + i]) * (U[size_x * j + i] - U[size_x * (j - 1) + i]) *
		                pow(dy_, -1));

		V_temp[size_x * j + i] = V[size_x * j + i] +
		    dt_ *
		        (c_ * ((V[size_x * j + i + 1] - 2 * V[size_x * j + i] + V[size_x * j + i - 1]) * pow(dx_, -2) +
		                  (V[size_x * (j + 1) + i] - 2 * V[size_x * j + i] + V[size_x * (j - 1) + i]) *
		                      pow(dy_, -2)) -
		            (ax_ + b_ * U[size_x * j + i]) * (V[size_x * j + i] - V[size_x * j + i - 1]) *
		                pow(dx_, -1) -
		            (ay_ + b_ * V[size_x * j + i]) * (V[size_x * j + i] - V[size_x * (j - 1) + i]) *
		                pow(dy_, -1));
	    }
	}

	// Make real variables same as temporal variables
	memcpy(&U[0], &U_temp[0], (size_x * size_y) * sizeof(double));
	memcpy(&V[0], &V_temp[0], (size_x * size_y) * sizeof(double));

	// Reestablish size as size of reduced submatrix for clarity
	size_x = fin_x - ini_x;
	size_y = fin_y - ini_y;

	// Initialise matrices that will hold expanded rows and columns
	double* U_send_x = new double[size_x]; // Row
	double* U_send_y = new double[size_y]; // Column

	double* V_send_x = new double[size_x]; // Row
	double* V_send_y = new double[size_y]; // Column

	if(row == 0 && col == 0) {

	    // Restart counter or index
	    k = 0;

	    // Go from first to last value in x - for elements in penultimate row only
	    j = size_y - 1;
	    for(i = 0; i < size_x; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 1) * j + i];
		V_send_x[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k++;
	    }

	    // Restart counter or index
	    k = 0;

	    // Go from first to penultimate value in y - for elements in penultimate column only
	    i = size_x - 1;
	    for(j = 0; j < size_y; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 1) * j + i];
		V_send_y[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix below
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);

	    // Send vectors to submatrix to the right
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

	    // Receive from matrix below and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Receive from matrix to the right and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in last column of U & V
	    i = size_x;
	    for(j = 0; j < size_y; j++) {
		U[(size_x + 1) * j + i] = U_send_y[k];
		V[(size_x + 1) * j + i] = V_send_y[k];
		k++;
	    }

	    // Determine where to copy U_send_x & V_send_x
	    int copy_place = (size_x + 1) * (size_y + 1 - 1);

	    // Copy U_send & V_send at the end of expanded matrix
	    memcpy(&U[copy_place], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[copy_place], &V_send_x[0], size_x * sizeof(double));

	    // Reset reduced matrix size variables as augmented matrix size variables
	    size_x = size_x + 1;
	    size_y = size_y + 1;
	}

	if(row == 0 && col > 0 && col < Px_ - 1) {

	    // Restart counter or index
	    k = 0;

	    // Go from 2nd to penultimate (avoid corners of augmented matrix)...
	    // values in x - for elements in penultimate row only
	    j = size_y - 1;
	    for(i = 1; i < size_x + 1; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 2) * j + i];
		V_send_x[k] = V[(size_x + 2) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix below
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);

	    k = 0; // Restart counter or index

	    // Go from first to penultimate value in y - for elements in second column only
	    i = 1;
	    for(j = 0; j < size_y; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 2) * j + i];
		V_send_y[k] = V[(size_x + 2) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix to the left
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

	    k = 0; // Restart counter or index

	    // Go from first to penultimate value in y - for elements in penultimate column only
	    i = size_x;
	    for(j = 0; j < size_y; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 2) * j + i];
		V_send_y[k] = V[(size_x + 2) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix to the right
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

	    // Receive from matrix below and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Receive from matrix to the left and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in first column of U & V
	    i = 0;
	    for(j = 0; j < size_y; j++) {
		U[(size_x + 2) * j + i] = U_send_y[k];
		V[(size_x + 2) * j + i] = V_send_y[k];
		k++;
	    }

	    // Determine where to copy U_send_x & V_send_x
	    copy_place = (size_x + 2) * (size_y + 1 - 1) + 1; // Note: +1 to avoid lower left corner

	    // Copy U_send & V_send at the end of expanded matrix
	    memcpy(&U[copy_place], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[copy_place], &V_send_x[0], size_x * sizeof(double));

	    // displayArray(U_send_x,size_x,1);

	    // displayArray(U,size_x+2,size_y+1);

	    // Receive from matrix to the right and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in last column of U & V
	    i = size_x + 1;
	    for(j = 0; j < size_y; j++) {
		U[(size_x + 2) * j + i] = U_send_y[k];
		V[(size_x + 2) * j + i] = V_send_y[k];
		k++;
	    }

	    // Reset reduced matrix size variables as augmented matrix size variables
	    size_x = size_x + 2;
	    size_y = size_y + 1;
	}

	if(row == 0 && col == Px_ - 1) {

	    // Restart counter or index
	    k = 0;

	    // Go from second to last value in x - for elements in penultimate row only
	    j = size_y - 1;
	    for(i = 1; i < size_x + 1; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 1) * j + i];
		V_send_x[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Restart counter or index
	    k = 0;

	    // Go from first to penultimate value in y - for elements in second column only
	    i = 1;
	    for(j = 0; j < size_y; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 1) * j + i];
		V_send_y[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix below
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);

	    // Send vectors to submatrix to the left
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

	    // Receive from matrix below and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Receive from matrix to the left and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in first column of U & V
	    i = 0;
	    for(j = 0; j < size_y; j++) {
		U[(size_x + 1) * j + i] = U_send_y[k];
		V[(size_x + 1) * j + i] = V_send_y[k];
		k++;
	    }

	    // Determine where to copy U_send_x & V_send_x
	    copy_place = (size_x + 1) * (size_y + 1 - 1) + 1; // Note: +1 to avoid corner

	    // Copy U_send & V_send at the end of expanded matrix
	    memcpy(&U[copy_place], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[copy_place], &V_send_x[0], size_x * sizeof(double));

	    // Reset reduced matrix size variables as augmented matrix size variables
	    size_x = size_x + 1;
	    size_y = size_y + 1;
	}

	if(row > 0 && row < Py_ - 1 && col == 0) {

	    // Restart counter or index
	    k = 0;

	    // Go from second to penultimate value in y - for elements in penultimate column only
	    i = size_x - 1;
	    for(j = 1; j < size_y + 1; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 1) * j + i];
		V_send_y[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Restart counter or index
	    k = 0;

	    // Go from first to penultimate value in x - for elements in second row only
	    j = 1;
	    for(i = 0; i < size_x; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 1) * j + i];
		V_send_x[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix above
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);

	    // Send vectors to submatrix to the right
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

	    // Restart counter or index
	    k = 0;

	    // Go from first to penultimate value in x - for elements in penultimate row only
	    j = size_y;
	    for(i = 0; i < size_x; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 1) * j + i];
		V_send_x[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix below
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);

	    // Receive from matrix above and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Receive from matrix to the right and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Determine where to copy U_send_x & V_send_x
	    int copy_place = 0;

	    // Copy U_send & V_send at the end of expanded matrix
	    memcpy(&U[copy_place], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[copy_place], &V_send_x[0], size_x * sizeof(double));

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in last column of U & V starting at row 1
	    i = size_x;
	    for(j = 1; j < size_y + 1; j++) {
		U[(size_x + 1) * j + i] = U_send_y[k];
		V[(size_x + 1) * j + i] = V_send_y[k];
		k++;
	    }

	    // Receive from matrix below and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Determine where to copy U_send_x & V_send_x
	    copy_place = (size_x + 1) * (size_y + 2 - 1);

	    // Copy U_send & V_send at the end of expanded matrix
	    memcpy(&U[copy_place], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[copy_place], &V_send_x[0], size_x * sizeof(double));

	    // Reset reduced matrix size variables as augmented matrix size variables
	    size_x = size_x + 1;
	    size_y = size_y + 2;
	}

	if(row > 0 && row < Py_ - 1 && col > 0 && col < Px_ - 1) {

	    // Restart counter or index
	    k = 0;

	    // Go from second to penultimate value in y - for elements in penultimate column only
	    i = size_x;
	    for(j = 1; j < size_y + 1; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 2) * j + i];
		V_send_y[k] = V[(size_x + 2) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Restart counter or index
	    k = 0;

	    // Go from second to penultimate value in x - for elements in penultimate row only
	    j = size_y;
	    for(i = 1; i < size_x + 1; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 2) * j + i];
		V_send_x[k] = V[(size_x + 2) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix to the right
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

	    // Send vectors to submatrix below
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);

	    // Restart counter or index
	    k = 0;

	    // Go from second to penultimate value in x - for elements in second row only
	    j = 1;
	    for(i = 1; i < size_x + 1; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 2) * j + i];
		V_send_x[k] = V[(size_x + 2) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    k = 0; // Restart counter or index

	    // Go from second to penultimate value in y - for elements in second column only
	    i = 1;
	    for(j = 1; j < size_x + 1; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 2) * j + i];
		V_send_y[k] = V[(size_x + 2) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix to the left
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

	    // Send vectors to submatrix above
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);

	    // displayArray(U_send_x,size_x,1);

	    // Receive from matrix above and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // displayArray(U_send_x,size_x,1);

	    // Receive from matrix to the left and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Determine where to copy U_send_x & V_send_x
	    copy_place = 1; // Skip corner

	    // Copy U_send_x & V_send_x at beginning +1 of expanded matrix
	    memcpy(&U[copy_place], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[copy_place], &V_send_x[0], size_x * sizeof(double));

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in first column of U & V starting at row 1
	    i = 0;
	    for(j = 1; j < size_y + 1; j++) {
		U[(size_x + 2) * j + i] = U_send_y[k];
		V[(size_x + 2) * j + i] = V_send_y[k];
		k++;
	    }

	    // Receive from matrix below and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Receive from matrix to the right and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Determine where to copy U_send_x & V_send_x
	    copy_place = (size_x + 2) * (size_y + 2 - 1) + 1; // Note: +1 to avoid lower left corner

	    // Copy U_send & V_send at the end of expanded matrix
	    memcpy(&U[copy_place], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[copy_place], &V_send_x[0], size_x * sizeof(double));

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in last column of U & V starting at row 1
	    i = size_x + 1;
	    for(j = 1; j < size_y + 1; j++) {
		U[(size_x + 2) * j + i] = U_send_y[k];
		V[(size_x + 2) * j + i] = V_send_y[k];
		k++;
	    }

	    // Reset reduced matrix size variables as augmented matrix size variables
	    size_x = size_x + 2;
	    size_y = size_y + 2;
	}

	if(row > 0 && row < Py_ - 1 && col == Px_ - 1) {

	    // Restart counter or index
	    k = 0;

	    // Go from second to last value in x - for elements in second row only
	    j = 1;
	    for(i = 1; i < size_x + 1; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 1) * j + i];
		V_send_x[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix above
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);

	    // Restart counter or index
	    k = 0;

	    // Go from second to last value in x - for elements in penultimate row only
	    j = size_y;
	    for(i = 1; i < size_x + 1; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 1) * j + i];
		V_send_x[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix below
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD);

	    // Restart counter or index
	    k = 0;

	    // Go from second to penultimate value in y - for elements in first column only
	    i = 0;
	    for(j = 1; j < size_y + 1; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 1) * j + i];
		V_send_y[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix to the left
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

	    // Receive from matrix above and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Receive from matrix to the left and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Determine where to copy U_send_x & V_send_x
	    copy_place = 1; // Avoid upper left corner

	    // Copy U_send & V_send at the beginning (+1) of expanded matrix
	    memcpy(&U[copy_place], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[copy_place], &V_send_x[0], size_x * sizeof(double));

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in first column of U & V starting at row 1
	    i = 0;
	    for(j = 1; j < size_y + 1; j++) {
		U[(size_x + 1) * j + i] = U_send_y[k];
		V[(size_x + 1) * j + i] = V_send_y[k];
		k++;
	    }

	    // Receive from matrix below and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank + Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Determine where to copy U_send_x & V_send_x
	    copy_place = (size_x + 1) * (size_y + 2 - 1) + 1; // Note: +1 to avoid lower left corner

	    // Copy U_send & V_send at the end of expanded matrix
	    memcpy(&U[copy_place], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[copy_place], &V_send_x[0], size_x * sizeof(double));

	    // Reset reduced matrix size variables as augmented matrix size variables
	    size_x = size_x + 1;
	    size_y = size_y + 2;
	}

	if(row == Py_ - 1 && col == 0) {

	    // Restart counter or index
	    k = 0;

	    // Go from first to penultimate value in x - for elements in second row only
	    j = 1;
	    for(i = 0; i < size_x; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 1) * j + i];
		V_send_x[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Restart counter or index
	    k = 0;

	    // Go from second to last value in y - for elements in penultimate column only
	    i = size_x - 1;
	    for(j = 1; j < size_y + 1; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 1) * j + i];
		V_send_y[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix to the right
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

	    // Send vectors to submatrix above
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);

	    // Receive from matrix above and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Receive from matrix to the right and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Copy U_send_x & V_send_x at the beginning of expanded matrix
	    memcpy(&U[0], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[0], &V_send_x[0], size_x * sizeof(double));

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in last column of U & V avoinding first row
	    i = size_x;
	    for(j = 1; j < size_y + 1; j++) {
		U[(size_x + 1) * j + i] = U_send_y[k];
		V[(size_x + 1) * j + i] = V_send_y[k];
		k++;
	    }

	    // Reset reduced matrix size variables as augmented matrix size variables
	    size_x = size_x + 1;
	    size_y = size_y + 1;
	}

	if(row == Py_ - 1 && col > 0 && col < Px_ - 1) {

	    // Restart counter or index
	    k = 0;

	    // Go from second to penultimate value in x - for elements in second row only
	    j = 1;
	    for(i = 1; i < size_x + 1; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 2) * j + i];
		V_send_x[k] = V[(size_x + 2) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix above
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);

	    // Restart counter or index
	    k = 0;

	    // Go from second to last value in y - for elements in second column only
	    i = 1;
	    for(j = 1; j < size_y + 1; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 2) * j + i];
		V_send_y[k] = V[(size_x + 2) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix to the left
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

	    // Restart counter or index
	    k = 0;

	    // Go from second to last value in y - for elements in penultimate column only
	    i = size_x;
	    for(j = 1; j < size_y + 1; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 2) * j + i];
		V_send_y[k] = V[(size_x + 2) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix to the left
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);

	    // Receive from matrix above and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Receive from matrix to the left and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Determine where to copy U_send_x & V_send_x
	    copy_place = 1; // Note: avoid top left corner

	    // Copy U_send & V_send at the end of expanded matrix
	    memcpy(&U[copy_place], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[copy_place], &V_send_x[0], size_x * sizeof(double));

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in first column of U & V skipping first row
	    i = 0;
	    for(j = 1; j < size_y + 1; j++) {
		U[(size_x + 2) * j + i] = U_send_y[k];
		V[(size_x + 2) * j + i] = V_send_y[k];
		k++;
	    }

	    // Receive from matrix to the right and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in last column of U & V skipping first row
	    i = size_x + 1;
	    for(j = 1; j < size_y + 1; j++) {
		U[(size_x + 2) * j + i] = U_send_y[k];
		V[(size_x + 2) * j + i] = V_send_y[k];
		k++;
	    }

	    // Reset reduced matrix size variables as augmented matrix size variables
	    size_x = size_x + 2;
	    size_y = size_y + 1;
	}

	if(row == Py_ - 1 && col == Px_ - 1) {

	    // Restart counter or index
	    k = 0;

	    // Go from second to last value in x - for elements in second row only
	    j = 1;
	    for(i = 1; i < size_x + 1; i++) {

		// Store values to be sent
		U_send_x[k] = U[(size_x + 1) * j + i];
		V_send_x[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Restart counter or index
	    k = 0;

	    // Go from second to last value in y - for elements in second column only
	    i = 1;
	    for(j = 1; j < size_y + 1; j++) {

		// Store values to be sent
		U_send_y[k] = U[(size_x + 1) * j + i];
		V_send_y[k] = V[(size_x + 1) * j + i];

		// Counter to define element index in array
		k = k + 1;
	    }

	    // Send vectors to submatrix to the left
	    MPI_Send(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

	    // Send vectors to submatrix above
	    MPI_Send(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);
	    MPI_Send(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD);

	    // Receive from matrix above and overwrite into U_send_x & V_send_x
	    MPI_Recv(U_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_x, size_x, MPI_DOUBLE, rank - Px_, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Receive from matrix to the left and overwrite into U_send_y & V_send_y
	    MPI_Recv(U_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    MPI_Recv(V_send_y, size_y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	    // Copy U_send_x & V_send_x at the beginning of expanded matrix avoiding first position
	    memcpy(&U[1], &U_send_x[0], size_x * sizeof(double));
	    memcpy(&V[1], &V_send_x[0], size_x * sizeof(double));

	    // Reset counter
	    k = 0;

	    // Insert U_send_y & V_send_y in first column of U & V avoinding first row
	    i = 0;
	    for(j = 1; j < size_y + 1; j++) {
		U[(size_x + 1) * j + i] = U_send_y[k];
		V[(size_x + 1) * j + i] = V_send_y[k];
		k++;
	    }
	}

	t = t + dt_;
    } // End time loop

} // End function

void Burgers::writeFile(int Nx_, int Ny_, int Px_, int Py_)
{

    using namespace std;

    int rank;
    int ini_x, fin_x, ini_y, fin_y, size_x, size_y, start_x, start_y, end_x, end_y;

    double* Receiver_U;
    double* Receiver_V;
    double* U_sorted;
    double* V_sorted;
    int* sizes;
    int* sizes_x;
    int* sizes_y;
    int* offset;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Submatrix row position
    int row = floor(rank / Px_);
    // Submatrix column position
    int col = rank - Px_ * row;

    // Define bounds for each submatrix in region I - first to penultimate rows & columns
    if(row < (Py_ - 1) && col < (Px_ - 1)) {

	ini_x = int(ceil((1.0 * Nx_) / (1.0 * Px_)) * col);
	fin_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * (col + 1));
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * (row + 1));
    }
    // Define bounds for each submatrix in region II - first to penultimate rows & last column
    if(row < (Py_ - 1) && col == (Px_ - 1)) {

	ini_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * col);
	fin_x = int(Nx_);
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * (row + 1));
    }
    // Define bounds for each submatrix in region III - last row & first to penultimate columns
    if(row == (Py_ - 1) && col < (Px_ - 1)) {
	ini_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * col);
	fin_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * (col + 1));
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(Ny_);
    }
    // Define bounds for each submatrix in region IV - last row & last column
    if(row == (Py_ - 1) && col == (Px_ - 1)) {

	ini_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * col);
	fin_x = int(Nx_);
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(Ny_);
    }

    int size_x_red = fin_x - ini_x;
    int size_y_red = fin_y - ini_y;
    int size_tot_red = size_x_red * size_y_red;

    // Determine size of augmented submatrix
    if(row == 0) {
	size_y = fin_y - ini_y + 1;
	start_y = 0;
	end_y = 1;
    }
    if(col == 0) {
	size_x = fin_x - ini_x + 1;
	start_x = 0;
	end_x = 1;
    }
    if(col > 0 && col < Px_ - 1) {
	size_x = fin_x - ini_x + 2;
	start_x = 1;
	end_x = 1;
    }
    if(row > 0 && row < Py_ - 1) {
	size_y = fin_y - ini_y + 2;
	start_y = 1;
	end_y = 1;
    }
    if(row == Py_ - 1) {
	size_y = fin_y - ini_y + 1;
	start_y = 1;
	end_y = 0;
    }
    if(col == Px_ - 1) {
	size_x = fin_x - ini_x + 1;
	start_x = 1;
	end_x = 0;
    }

    double* U_red = new double[size_tot_red];
    double* V_red = new double[size_tot_red];

    int k = 0;
    int l = 0;
    int j, i;

    // Reduce augmented matrices
    for(int j = start_y; j < size_y - end_y; j++) {
	for(int i = start_x; i < size_x - end_x; i++) {
	    U_red[size_x_red * l + k] = U[size_x * j + i];
	    V_red[size_x_red * l + k] = V[size_x * j + i];
	    // cout << U_red[size_x_red * l + k] << endl;
	    k++;
	}
	k = 0;
	l++;
    }

    if(rank == 0) {
	Receiver_U = new double[Nx_ * Ny_];
	Receiver_V = new double[Nx_ * Ny_];
	sizes = new int[Px_ * Py_];
	sizes_x = new int[Px_ * Py_];
	sizes_y = new int[Px_ * Py_];
	offset = new int[Px_ * Py_];
	U_sorted = new double[Nx_ * Ny_];
	V_sorted = new double[Nx_ * Ny_];
    }

    MPI_Gather(&size_tot_red, 1, MPI_INT, sizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&size_x_red, 1, MPI_INT, sizes_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&size_y_red, 1, MPI_INT, sizes_y, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank == 0) {
	for(int j = 0; j < Px_ * Py_; j++) {
	    offset[j] = 0;
	    for(int i = 0; i < j; i++) {
		offset[j] += sizes[i];
	    }
	}
    }

    MPI_Gatherv(U_red, size_tot_red, MPI_DOUBLE, Receiver_U, sizes, offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(V_red, size_tot_red, MPI_DOUBLE, Receiver_V, sizes, offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0) {

	j = 0;
	l = 0;
	cout << "sdfggegrgrgre" << endl;
	for(k = 0; k < Px_ * Py_; k++) {
	    for(i = 0; i < sizes_x[k]; i++) {
		U_sorted[Nx_ * k + l] = Receiver_U[j + i];
		V_sorted[Nx_ * k + l] = Receiver_V[j + i];
		l++;
	    }
	    j = +sizes[k];
	    l = 0;
	}

	ofstream vOut("velField.txt", ios::out | ios::trunc);
	if(vOut.good()) {
	    for(int j = 0; j < Ny_; j++) {
		for(int i = 0; i < Nx_; i++) {
		    vOut << "Element"
		         << "(" << j << ", " << i << ")"
		         << " for (U, V) is "
		         << "(" << U_sorted[Nx_ * j + i] << ", " << V_sorted[Nx_ * j + i] << ")" << endl;
		}
	    }
	}
	vOut.close();
	return;
    }
}

double Burgers::energy(int Nx_, int Ny_, double dx_, double dy_, int Px_, int Py_)
{
    int rank;
    int ini_x, fin_x, ini_y, fin_y, size_x, size_y, start_x, start_y, end_x, end_y;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Submatrix row position
    int row = floor(rank / Px_);
    // Submatrix column position
    int col = rank - Px_ * row;

    // Define bounds for each submatrix in region I - first to penultimate rows & columns
    if(row < (Py_ - 1) && col < (Px_ - 1)) {

	ini_x = int(ceil((1.0 * Nx_) / (1.0 * Px_)) * col);
	fin_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * (col + 1));
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * (row + 1));
    }
    // Define bounds for each submatrix in region II - first to penultimate rows & last column
    if(row < (Py_ - 1) && col == (Px_ - 1)) {

	ini_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * col);
	fin_x = int(Nx_);
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * (row + 1));
    }
    // Define bounds for each submatrix in region III - last row & first to penultimate columns
    if(row == (Py_ - 1) && col < (Px_ - 1)) {
	ini_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * col);
	fin_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * (col + 1));
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(Ny_);
    }
    // Define bounds for each submatrix in region IV - last row & last column
    if(row == (Py_ - 1) && col == (Px_ - 1)) {

	ini_x = int(ceil((Nx_ * 1.0) / (Px_ * 1.0)) * col);
	fin_x = int(Nx_);
	ini_y = int(ceil((Ny_ * 1.0) / (Py_ * 1.0)) * row);
	fin_y = int(Ny_);
    }

    int size_x_red = fin_x - ini_x;
    int size_y_red = fin_y - ini_y;

    // Determine size of augmented submatrix
    if(row == 0) {
	size_y = fin_y - ini_y + 1;
	start_y = 0;
	end_y = 1;
    }
    if(col == 0) {
	size_x = fin_x - ini_x + 1;
	start_x = 0;
	end_x = 1;
    }
    if(col > 0 && col < Px_ - 1) {
	size_x = fin_x - ini_x + 2;
	start_x = 1;
	end_x = 1;
    }
    if(row > 0 && row < Py_ - 1) {
	size_y = fin_y - ini_y + 2;
	start_y = 1;
	end_y = 1;
    }
    if(row == Py_ - 1) {
	size_y = fin_y - ini_y + 1;
	start_y = 1;
	end_y = 0;
    }
    if(col == Px_ - 1) {
	size_x = fin_x - ini_x + 1;
	start_x = 1;
	end_x = 0;
    }

    int k;
    int l;

    double energy1 = 0;
    double* E = new double[(size_x_red) * (size_y_red)];

    // Calculate from start to end positions of each augmented submatrix
    for(int j = start_y, k = 0; j < (size_y - end_y); j++, k++) {
	for(int i = start_x, l = 0; i < (size_x - end_x); i++, l++) {
	    E[(size_x_red)*k + l] = pow(U[size_x * j + i], 2) + pow(V[size_x * j + i], 2);
	}
    }

    for(int j = 0; j < size_x_red; j++) {
	for(int i = 0; i < size_y_red; i++) {

	    energy1 = +0.25 * dx_ * dy_ * (E[size_x_red * j + i] + E[size_x_red * j + i + 1] +
	                                      E[size_x_red * (j + 1) + i] + E[size_x_red * (j + 1) + i + 1]);
	}
    }

    energy1 = 0.5 * energy1;

    double energy_final;

    MPI_Reduce(&energy1, &energy_final, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0) {
	return energy_final;
    }
}