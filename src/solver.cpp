
#include "solver.h"

/*!
 * \brief Solver implementation 
 */
void lbmSolver::solveFlow(lbmLattice* argLattice, input* argInput){

	double u_max, u_max_prev = 0.0, rel_residual;

	/// Assign to arguments to local pointers
	lattice = argLattice;
	inputData = argInput;
	
	/// Number of interior lattice points x-dir
	int Nx = inputData->getNx();

	/// Number of interior lattice points y-dir
	int Ny = inputData->getNy();	

	/// Generate porous structure based on given porosity function
	PorosityFunction(Nx, Ny);

	/// Loop through number of iterations
	for(int i=0;i<inputData->getNIter();i++){
		
		/// Compute Macroscopic variables
		u_max = ComputeMacro(Nx, Ny);

		/// Compute Equilibrium Distribution
		ComputeEqDistribution(Nx, Ny);

		/// Collision Operator
		Collision(Nx, Ny);

		/// Source Term
		SourceTerm(Nx, Ny);

		/// Apply Boundary Condition
		BoundaryCondition(Nx, Ny);

		/// Apply Porous Boundary condition
		ApplyPorousBC(Nx, Ny);

		/// Streaming operator
		Streaming(Nx, Ny);

		rel_residual = fabs((u_max_prev - u_max)/u_max);
		cout<<">Iteration No. : "<<i<<"\tRelative residual : "<<rel_residual<<"\tMax Velocity : "<<setprecision(17)<<u_max<<endl;

		/// Convergence criteria
		if(rel_residual< 1e-16){
			cout<<"\n> Simulation reached convergence! \n"<<endl;
			break;
		}
		
		u_max_prev = u_max;
	}

	/// Extract data
	ExtractVel(Nx, Ny);

return;
}

/*!
 * \brief Computation of macroscopic flow variables such as density and velocity 
 */
double lbmSolver::ComputeMacro(const int Nx, const int Ny){
	
	double* f = new double[9]();
	double Rho, u_max=0.0;

	// Loop through all interior lattice points
	for(int j=1;j<Ny+1;j++){
		for(int i=1;i<Nx+1;i++){

		   if(lattice->latpoint[(Nx+2)*j + i].P == 1){
			std::memcpy(f, lattice->latpoint[(Nx+2)*j + i].f, 9*sizeof(double));

			// Calculate density from the distribution
			Rho = 0.0;
			for(int k=0;k<9;k++)	Rho += f[k];

			lattice->latpoint[(Nx+2)*j + i].rho = Rho;

			// Calculate velocity components from distribution
			lattice->latpoint[(Nx+2)*j + i].u[0] = (f[1] - f[3] + f[5] - f[6] - f[7] + f[8])/Rho;
			lattice->latpoint[(Nx+2)*j + i].u[1] = (f[2] - f[4] + f[5] + f[6] - f[7] - f[8])/Rho;

			if(lattice->latpoint[(Nx+2)*j + i].u[0]>u_max)
				u_max = lattice->latpoint[(Nx+2)*j + i].u[0];
		   }else{
			lattice->latpoint[(Nx+2)*j + i].rho = 0.0;
			lattice->latpoint[(Nx+2)*j + i].u[0] = 0.0;
			lattice->latpoint[(Nx+2)*j + i].u[1] = 0.0;
		   }
		}
	}

	delete [] f;
return u_max;
}

/*!
 * \brief Computation of Equilibrium distribution 
 */
void lbmSolver::ComputeEqDistribution(const int Nx, const int Ny){
	
	double Rho, Ux, Uy, U_sqr;
	
	// Loop through all interior lattice points
	for(int j=1;j<Ny+1;j++){
		for(int i=1;i<Nx+1;i++){

		   if(lattice->latpoint[(Nx+2)*j + i].P == 1){

			Rho = lattice->latpoint[(Nx+2)*j + i].rho;
			Ux = lattice->latpoint[(Nx+2)*j + i].u[0];
			Uy = lattice->latpoint[(Nx+2)*j + i].u[1];
			U_sqr = (Ux*Ux + Uy*Uy)/(2*cs_sqr);	
			
			lattice->latpoint[(Nx+2)*j + i].f_eq[0] = Rho*w_a*( 1.0 - U_sqr );

			lattice->latpoint[(Nx+2)*j + i].f_eq[1] = Rho*w_b*( 1.0 + (Ux/cs_sqr) - U_sqr + (Ux*Ux/(2*cs_sqr*cs_sqr)));
			lattice->latpoint[(Nx+2)*j + i].f_eq[2] = Rho*w_b*( 1.0 + (Uy/cs_sqr) - U_sqr + (Uy*Uy/(2*cs_sqr*cs_sqr)));
			lattice->latpoint[(Nx+2)*j + i].f_eq[3] = Rho*w_b*( 1.0 - (Ux/cs_sqr) - U_sqr + (Ux*Ux/(2*cs_sqr*cs_sqr)));
			lattice->latpoint[(Nx+2)*j + i].f_eq[4] = Rho*w_b*( 1.0 - (Uy/cs_sqr) - U_sqr + (Uy*Uy/(2*cs_sqr*cs_sqr)));

			lattice->latpoint[(Nx+2)*j + i].f_eq[5] = Rho*w_c*( 1.0 + ((Ux+Uy)/cs_sqr) - U_sqr + ((Ux+Uy)*(Ux+Uy)/(2*cs_sqr*cs_sqr)));
			lattice->latpoint[(Nx+2)*j + i].f_eq[6] = Rho*w_c*( 1.0 + ((-Ux+Uy)/cs_sqr) - U_sqr + ((-Ux+Uy)*(-Ux+Uy)/(2*cs_sqr*cs_sqr)));
			lattice->latpoint[(Nx+2)*j + i].f_eq[7] = Rho*w_c*( 1.0 - ((Ux+Uy)/cs_sqr) - U_sqr + ((Ux+Uy)*(Ux+Uy)/(2*cs_sqr*cs_sqr)));
			lattice->latpoint[(Nx+2)*j + i].f_eq[8] = Rho*w_c*( 1.0 + ((Ux-Uy)/cs_sqr) - U_sqr + ((Ux-Uy)*(Ux-Uy)/(2*cs_sqr*cs_sqr)));
		   }
		}
	}

return;
}

/*!
 * \brief Collision operator 
 */
void lbmSolver::Collision(const int Nx, const int Ny){

	double omega = inputData->getOmega();
	double* f = new double[9]();
	double* f_eq = new double[9]();
	
	// Loop through all interior lattice points
	for(int j=1;j<Ny+1;j++){
		for(int i=1;i<Nx+1;i++){

		   if(lattice->latpoint[(Nx+2)*j + i].P == 1){

			std::memcpy(f, lattice->latpoint[(Nx+2)*j + i].f, 9*sizeof(double));
			std::memcpy(f_eq, lattice->latpoint[(Nx+2)*j + i].f_eq, 9*sizeof(double));	
			for(int k=0;k<9;k++)
				lattice->latpoint[(Nx+2)*j + i].f[k] += omega*(f_eq[k] - f[k]);
		   }

		}
	}
	
	delete [] f;
	delete [] f_eq;
return;
}

/*!
 * \brief Adding Source term 
 */
void lbmSolver::SourceTerm(const int Nx, const int Ny){

	double omega = inputData->getOmega();
	double H = inputData->getNy();
	double Rho_0 = inputData->getDen();
	double Ux = inputData->getUx();
	
	// Viscosity of fluid
	double mu = Rho_0*(1.0/omega - 0.5)*cs_sqr;
	
	// Maximum target velocity
	double u_m = Ux;

	// Pressure gradient
	double G = 8*mu*u_m/(H*H);
	
	
	// Loop through all interior lattice points
	for(int j=1;j<Ny+1;j++){
		for(int i=1;i<Nx+1;i++){

		   if(lattice->latpoint[(Nx+2)*j + i].P == 1){

			lattice->latpoint[(Nx+2)*j + i].f[1] += G/6;
			lattice->latpoint[(Nx+2)*j + i].f[5] += G/6;
			lattice->latpoint[(Nx+2)*j + i].f[8] += G/6;

			lattice->latpoint[(Nx+2)*j + i].f[3] -= G/6;
			lattice->latpoint[(Nx+2)*j + i].f[6] -= G/6;
			lattice->latpoint[(Nx+2)*j + i].f[7] -= G/6;
		   }

		}
	}
return;
}

/*!
 * \brief Applying boundary condition
 */
void lbmSolver::BoundaryCondition(const int Nx, const int Ny){
	
	// Buffer f
	double* f = new double[9]();

	// West and East neighbours
	int W, E;

	// Bounceback condition at top and bottom walls
	for(int i=1;i<Nx+1;i++){
		if(i==0) W = (i+Nx+1)%(Nx+2) - 1; else W = (i+Nx+1)%(Nx+2);
		if(i==Nx+1) E = (i+1)%(Nx+2) + 1; else E = (i+1)%(Nx+2);

		lattice->latpoint[(Nx+2)*0 + i].f[2] = lattice->latpoint[(Nx+2)*1 + i].f[4];
		lattice->latpoint[(Nx+2)*0 + i].f[5] = lattice->latpoint[(Nx+2)*1 + E].f[7];
		lattice->latpoint[(Nx+2)*0 + i].f[6] = lattice->latpoint[(Nx+2)*1 + W].f[8];

		lattice->latpoint[(Nx+2)*(Ny+1) + i].f[4] = lattice->latpoint[(Nx+2)*Ny + i].f[2];
		lattice->latpoint[(Nx+2)*(Ny+1) + i].f[7] = lattice->latpoint[(Nx+2)*Ny + W].f[5];
		lattice->latpoint[(Nx+2)*(Ny+1) + i].f[8] = lattice->latpoint[(Nx+2)*Ny + E].f[6];
	}

	// Bounceback condition at corners
	lattice->latpoint[(Nx+2)*0 + 0].f[5] = lattice->latpoint[(Nx+2)*1 + 1].f[7];
	lattice->latpoint[(Nx+2)*(Ny+1) + 0].f[8] = lattice->latpoint[(Nx+2)*Ny + 1].f[6];
	lattice->latpoint[(Nx+2)*(Ny+1) + (Nx+1)].f[7] = lattice->latpoint[(Nx+2)*Ny + Nx].f[5];
	lattice->latpoint[(Nx+2)*0 + (Nx+1)].f[6] = lattice->latpoint[(Nx+2)*1 + Nx].f[8];			

	
	// Periodic BC at left and right walls
	for(int j=1;j<Ny+1;j++){
		lattice->latpoint[(Nx+2)*j + 0].f[1] = lattice->latpoint[(Nx+2)*j + Nx].f[1];
		lattice->latpoint[(Nx+2)*j + 0].f[5] = lattice->latpoint[(Nx+2)*j + Nx].f[5];
		lattice->latpoint[(Nx+2)*j + 0].f[8] = lattice->latpoint[(Nx+2)*j + Nx].f[8];

		lattice->latpoint[(Nx+2)*j + Nx + 1].f[3] = lattice->latpoint[(Nx+2)*j + 1].f[3];
		lattice->latpoint[(Nx+2)*j + Nx + 1].f[6] = lattice->latpoint[(Nx+2)*j + 1].f[6];
		lattice->latpoint[(Nx+2)*j + Nx + 1].f[7] = lattice->latpoint[(Nx+2)*j + 1].f[7];
	}

	delete [] f;
return;
}

/*!
 * \brief Streaming operator 1
 */
void lbmSolver::Streaming(const int Nx, const int Ny){

	
	// Interior points
	for(int j=1;j<Ny+1;j++){
		for(int i=1;i<Nx+1;i++){
			lattice->latpoint[(Nx+2)*j + i].f_buf[0] = lattice->latpoint[(Nx+2)*j + i].f[0];

			lattice->latpoint[(Nx+2)*j + i + 1].f_buf[1] = lattice->latpoint[(Nx+2)*j + i].f[1];	
			lattice->latpoint[(Nx+2)*(j+1) + i].f_buf[2] = lattice->latpoint[(Nx+2)*j + i].f[2];
			lattice->latpoint[(Nx+2)*j + i - 1].f_buf[3] = lattice->latpoint[(Nx+2)*j + i].f[3];
			lattice->latpoint[(Nx+2)*(j-1) + i].f_buf[4] = lattice->latpoint[(Nx+2)*j + i].f[4];

			lattice->latpoint[(Nx+2)*(j+1) + i + 1].f_buf[5] = lattice->latpoint[(Nx+2)*j + i].f[5];
			lattice->latpoint[(Nx+2)*(j+1) + i - 1].f_buf[6] = lattice->latpoint[(Nx+2)*j + i].f[6];
			lattice->latpoint[(Nx+2)*(j-1) + i - 1].f_buf[7] = lattice->latpoint[(Nx+2)*j + i].f[7];
			lattice->latpoint[(Nx+2)*(j-1) + i + 1].f_buf[8] = lattice->latpoint[(Nx+2)*j + i].f[8];
		}
	}
			
	// Corner points
	lattice->latpoint[(Nx+2)*1 + 1].f_buf[5] = lattice->latpoint[(Nx+2)*0 + 0].f[5];		//i=0, j=0
	lattice->latpoint[(Nx+2)*Ny + 1].f_buf[8] = lattice->latpoint[(Nx+2)*(Ny+1) + 0].f[8];		//i=0, j=Ny+1
	lattice->latpoint[(Nx+2)*Ny + Nx].f_buf[7] = lattice->latpoint[(Nx+2)*(Ny+1) + (Nx+1)].f[7];	//i=Nx+1, j=Ny+1
	lattice->latpoint[(Nx+2)*1 + Nx].f_buf[6] = lattice->latpoint[(Nx+2)*0 + (Nx+1)].f[6];		//i=Nx+1, j=0

	// Bottom wall
	for(int i=1;i<Nx+1;i++){
		lattice->latpoint[(Nx+2)*1 + i].f_buf[2] = lattice->latpoint[(Nx+2)*0 + i].f[2];
		lattice->latpoint[(Nx+2)*1 + i + 1].f_buf[5] = lattice->latpoint[(Nx+2)*0 + i].f[5];
		lattice->latpoint[(Nx+2)*1 + i - 1].f_buf[6] = lattice->latpoint[(Nx+2)*0 + i].f[6];
	}

	// Top wall
	for(int i=1;i<Nx+1;i++){
		lattice->latpoint[(Nx+2)*Ny + i].f_buf[4] = lattice->latpoint[(Nx+2)*(Ny+1) + i].f[4];
		lattice->latpoint[(Nx+2)*Ny + i - 1].f_buf[7] = lattice->latpoint[(Nx+2)*(Ny+1) + i].f[7];
		lattice->latpoint[(Nx+2)*Ny + i + 1].f_buf[8] = lattice->latpoint[(Nx+2)*(Ny+1) + i].f[8];
	}

	// Left wall
	for(int j=1;j<Ny+1;j++){
		lattice->latpoint[(Nx+2)*j + 1].f_buf[1] = lattice->latpoint[(Nx+2)*j + 0].f[1];
		lattice->latpoint[(Nx+2)*(j+1) + 1].f_buf[5] = lattice->latpoint[(Nx+2)*j + 0].f[5];
		lattice->latpoint[(Nx+2)*(j-1) + 1].f_buf[8] = lattice->latpoint[(Nx+2)*j + 0].f[8];
	}
			
	// Right wall
	for(int j=1;j<Ny+1;j++){
		lattice->latpoint[(Nx+2)*j + Nx].f_buf[3] = lattice->latpoint[(Nx+2)*j + (Nx+1)].f[3];
		lattice->latpoint[(Nx+2)*(j+1) + Nx].f_buf[6] = lattice->latpoint[(Nx+2)*j + (Nx+1)].f[6];
		lattice->latpoint[(Nx+2)*(j-1) + Nx].f_buf[7] = lattice->latpoint[(Nx+2)*j + (Nx+1)].f[7];
	}

	// Copy from buffer to original
	for(int j=0;j<Ny+2;j++){
		for(int i=0;i<Nx+2;i++){
			std::memcpy(lattice->latpoint[(Nx+2)*j + i].f, lattice->latpoint[(Nx+2)*j + i].f_buf, 9*sizeof(double));
		}
	}

return;
}

/*!
 * \brief Generate the structure for solid and fluid nodes based on porosity function
 */
void lbmSolver::PorosityFunction(const int Nx, const int Ny){

	double im = Nx/2;
	double jm = Ny/2;
	double r;

	for(int i=1;i<Nx+1;i++){
		for(int j=1;j<Ny+1;j++){
        		r = ((i-im)/Nx)*((i-im)/Nx) + ((j-jm)/Ny)*((j-jm)/Ny);
       			if(r<0.05)	
				lattice->latpoint[(Nx+2)*j + i].P = 0;
		}
	}

	/// Random porosity function
/*	for(int i=1;i<Nx+1;i++){
		for(int j=1;j<Ny+1;j++){
			r = (double) rand() / (RAND_MAX);
			cout<<r<<endl;
       			if(i>im-im/2 && i<im+im/2)	
				if(r>0.9) lattice->latpoint[(Nx+2)*j + i].P = 0;
		}
	}*/
return;
}

/*!
 * \brief Applying porous boundary condition
 */
void lbmSolver::ApplyPorousBC(const int Nx, const int Ny){

for(int i=1;i<Nx+1;i++){
	for(int j=1;j<Ny+1;j++){
        
		if(lattice->latpoint[(Nx+2)*j + i].P == 0){
			if(lattice->latpoint[(Nx+2)*j + i + 1].P == 1)
       				lattice->latpoint[(Nx+2)*j + i].f[1] = lattice->latpoint[(Nx+2)*j + i + 1].f[3];
			if(lattice->latpoint[(Nx+2)*(j + 1) + i].P == 1)
       				lattice->latpoint[(Nx+2)*j + i].f[2] = lattice->latpoint[(Nx+2)*(j + 1) + i].f[4];
			if(lattice->latpoint[(Nx+2)*j + i - 1].P == 1)
       				lattice->latpoint[(Nx+2)*j + i].f[3] = lattice->latpoint[(Nx+2)*j + i - 1].f[1];
			if(lattice->latpoint[(Nx+2)*(j - 1) + i].P == 1)
       				lattice->latpoint[(Nx+2)*j + i].f[4] = lattice->latpoint[(Nx+2)*(j - 1) + i].f[2];
			if(lattice->latpoint[(Nx+2)*(j + 1) + i + 1].P == 1)
       				lattice->latpoint[(Nx+2)*j + i].f[5] = lattice->latpoint[(Nx+2)*(j + 1) + i + 1].f[7];
			if(lattice->latpoint[(Nx+2)*(j + 1) + i - 1].P == 1)
       				lattice->latpoint[(Nx+2)*j + i].f[6] = lattice->latpoint[(Nx+2)*(j + 1) + i - 1].f[8];
			if(lattice->latpoint[(Nx+2)*(j - 1) + i - 1].P == 1)
       				lattice->latpoint[(Nx+2)*j + i].f[7] = lattice->latpoint[(Nx+2)*(j - 1) + i - 1].f[5];
			if(lattice->latpoint[(Nx+2)*(j - 1) + i + 1].P == 1)
       				lattice->latpoint[(Nx+2)*j + i].f[8] = lattice->latpoint[(Nx+2)*(j - 1) + i + 1].f[6];
		}
	}
}

return;
}

/*!
 * \brief print density and velocity at lattice points
 */
void lbmSolver::ExtractVel(const int Nx, const int Ny){

	string dir = inputData->getTitle();
	string dummy;

    	mkdir(dir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	ofstream file, ux, uy, nx, ny;

	dummy = dir;
	ux.open(dummy.append("/Ux.dat").c_str(),ios::out);
	dummy = dir;
	uy.open(dummy.append("/Uy.dat").c_str(),ios::out);
	dummy = dir;
	nx.open(dummy.append("/nx.dat").c_str(),ios::out);
	dummy = dir;
	ny.open(dummy.append("/ny.dat").c_str(),ios::out);
	dummy = dir;
	file.open(dummy.append("/rho.dat").c_str(),ios::out);

	for(int j=1;j<Ny+1;j++){

		for(int i=1;i<Nx+1;i++){
			file<<lattice->latpoint[(Nx+2)*j + i].rho<<"\t";
			ux<<lattice->latpoint[(Nx+2)*j + i].u[0]<<"\t";
			uy<<lattice->latpoint[(Nx+2)*j + i].u[1]<<"\t";
			nx<<i<<"\t";
			ny<<j<<"\t";
		}
		file<<endl;
		ux<<endl;
		uy<<endl;
		nx<<endl;
		ny<<endl;
	}

	file.close();
	ux.close();
	uy.close();
	nx.close();
	ny.close();
return;
}
