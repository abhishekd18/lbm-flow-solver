#ifndef SOLVER_H_
#define SOLVER_H_

#include "input.h"
#include "lattice.h"

/*!
 * \brief This class defines the flow solver and member functions
 */

class lbmSolver
{
    private:
	lbmLattice* lattice;	// Local pointer for lattice
	input*	    inputData;	// Local pointer for input data
	
	/// Private Member functions
	double ComputeMacro(const int Nx, const int Ny);
	void ComputeEqDistribution(const int Nx, const int Ny);
	void Collision(const int Nx, const int Ny);
	void SourceTerm(const int Nx, const int Ny);
	void BoundaryCondition(const int Nx, const int Ny);
	void Streaming(const int Nx, const int Ny);
	void ExtractVel(const int Nx, const int Ny);

    protected:

    public:
        /// Default constructor
        lbmSolver(){};

        /// Destructor
        ~lbmSolver(){};

        /// Member functions
	void solveFlow(lbmLattice*, input*);
};

#endif


