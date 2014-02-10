#ifndef INPUT_H_
#define INPUT_H_

#include "constants.h"

/*! 
 * \brief This class is used for reading and storig the data from input file 
 */
class input
{
    private:
        /// PRIVATE VARIABLES
        string  title;      // title of the document
	string	wdir;	    // working directory
	double  lx;	    // length of domain in x-dir
	double  ly;	    // length of domain in x-dir
	int	Nx;	    // number of interior lattice points in x-dir
	int	Ny;	    // number of interior lattice points in y-dir
	double  den;        // density of fluid
	double  ux;	    // x-component of velocity
	double  uy;	    // y-component of velocity
        double  omega;      // viscosity parameter
        int     nIter;      // number of iterations
        
    protected:

    public:
        /// DEFAULT CONSTRUCTOR ///
        input();

        /// DESTRUCTOR
        ~input(){};

        /// GETTERS ///  
        string          getTitle()      {return title;};
        string          getWdir()       {return wdir;};
        double          getLx()         {return lx;};
        double          getLy()         {return ly;};
	int		getNx()		{return Nx;};
	int		getNy()		{return Ny;};
        double          getDen()        {return den;};
        double          getUx()         {return ux;};
        double          getUy()         {return uy;};
        double          getOmega()      {return omega;};
        int             getNIter()      {return nIter;};

        /// PUBLIC INTERFACE METHOD
        void readInputFile();
};

#endif
