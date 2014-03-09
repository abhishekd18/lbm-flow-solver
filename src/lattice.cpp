
#include "lattice.h"

/*!
 *	Create and initialize the lattice points with structure to hold flow variables
 */
void lbmLattice::createLattice(input* inputData){

	// Number of interior lattice points x-dir
	int Nx = inputData->getNx();

	// Number of interior lattice points y-dir
	int Ny = inputData->getNy();

	// Create the lattice points including ghost points
	latpoint = new latticePoint[(Nx+1)*(Ny+1)]();

	// Initialize the lattice variables
	
	// Rho
	double Rho = inputData->getDen();
	
	// Velocity components
	double Ux = inputData->getUx();
	double Uy = inputData->getUy();
	double U_sqr;

	// Initialize flow variables
	for(int j=1;j<Ny;j++){
		for(int i=1;i<Nx;i++){

			latpoint[(Nx+1)*j + i].rho = Rho;

			latpoint[(Nx+1)*j + i].P = 1;

			latpoint[(Nx+1)*j + i].u[0] = Ux;
			latpoint[(Nx+1)*j + i].u[1] = Uy;

			U_sqr = (Ux*Ux + Uy*Uy)/(2*cs_sqr);	
			
			// Initialize the distribution f
			latpoint[(Nx+1)*j + i].f[0] = Rho*w_a*( 1 - U_sqr );

			latpoint[(Nx+1)*j + i].f[1] = Rho*w_b*( 1 + (Ux/cs_sqr) - U_sqr + (Ux*Ux/(2*cs_sqr*cs_sqr)));
			latpoint[(Nx+1)*j + i].f[2] = Rho*w_b*( 1 + (Uy/cs_sqr) - U_sqr + (Uy*Uy/(2*cs_sqr*cs_sqr)));
			latpoint[(Nx+1)*j + i].f[3] = Rho*w_b*( 1 - (Ux/cs_sqr) - U_sqr + (Ux*Ux/(2*cs_sqr*cs_sqr)));
			latpoint[(Nx+1)*j + i].f[4] = Rho*w_b*( 1 - (Uy/cs_sqr) - U_sqr + (Uy*Uy/(2*cs_sqr*cs_sqr)));

			latpoint[(Nx+1)*j + i].f[5] = Rho*w_c*( 1 + ((Ux+Uy)/cs_sqr) - U_sqr + ((Ux+Uy)*(Ux+Uy)/(2*cs_sqr*cs_sqr)));
			latpoint[(Nx+1)*j + i].f[6] = Rho*w_c*( 1 + ((-Ux+Uy)/cs_sqr) - U_sqr + ((-Ux+Uy)*(-Ux+Uy)/(2*cs_sqr*cs_sqr)));
			latpoint[(Nx+1)*j + i].f[7] = Rho*w_c*( 1 - ((Ux+Uy)/cs_sqr) - U_sqr + ((Ux+Uy)*(Ux+Uy)/(2*cs_sqr*cs_sqr)));
			latpoint[(Nx+1)*j + i].f[8] = Rho*w_c*( 1 + ((Ux-Uy)/cs_sqr) - U_sqr + ((Ux-Uy)*(Ux-Uy)/(2*cs_sqr*cs_sqr)));
		}
	}

return;
}
