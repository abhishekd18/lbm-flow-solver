/*!
 *	Class definitions for lattice points 
*/

#ifndef LATTICE_H_
#define LATTICE_H_

#include "input.h"

/*!
 *	Class for lattice point level structure
 */
class latticePoint
{
	private:

	protected:

	public:
	double* f_eq;
	double* f;
	double* f_buf;
	double  rho;
	double* u;
	int     P;
	double  x;
	double  y;

	/// Default constructor
	latticePoint(){
		rho  = 0.0;
		f_eq = new double[9]();
		f    = new double[9]();
		f_buf= new double[9]();
		u    = new double[2]();	
	};

	/// Default destructor
	~latticePoint(){
		delete [] f_eq;
		delete [] f;
		delete [] f_buf;
		delete [] u;
	};			

};

/*!
 *	Class for creating Lattice structure
 */ 
class lbmLattice
{
	private:

	protected:

	public:
	latticePoint*	latpoint;

	/// Default constructor
	lbmLattice(){};

	/// Default destructor
	~lbmLattice(){delete [] latpoint;};			
	
	/// public interface
	void createLattice(input*);
};


#endif
