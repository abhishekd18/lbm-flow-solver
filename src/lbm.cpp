/*!
 *	This flow solver is based on Lattice Boltzmann Methods.
*/

#include "input.h"
#include "lattice.h"
#include "solver.h"

using namespace std;

int main(int argc, char **argv){

	input* 		inputData = new input;
	lbmLattice*     lattice   = new lbmLattice;
	lbmSolver*	solver	  = new lbmSolver;

	/// Read input file
	inputData->readInputFile();

	/// Create lattice structure and initialize with flow variables
	lattice->createLattice(inputData);	

	/// Solve for flow variables
	solver->solveFlow(lattice, inputData);

	delete inputData;
	delete lattice;
	delete solver;

return 0;
}
