Copyright @ Abhishek Deshmukh 2014 - 2099

****************************************************************************************************
2D LATTICE BOLTZMANN FLOW SOLVER
****************************************************************************************************


****************************************************************************************************
DESCRIPTION OF THE CODE
****************************************************************************************************
The code is based on Lattice Boltzmann Methods to solve Navier-Stokes equations in 2D. D2Q9 Lattice
is used. Two files, namely, Ux.dat and rho.dat, are written in the current folder as output.

****************************************************************************************************
COMPILING AND USING THE CODE (IN LINUX)
****************************************************************************************************
* Open a terminal
* Go to src folder in the terminal
* Type "make" in the terminal
* An executable named "lbm" will be created in the current folder
* Define the settings in input file "input"
* Type "./lbm" to run the program.

****************************************************************************************************
EXAMPLE INPUT FILE (settings.in)
****************************************************************************************************
# Title of the simulation
title 2D_Channel

# Working directory "../run/" or "./"
wdir ./

# Length of domain in x-direction
lx 2.0

# Length of domain in y-direction
ly 1.0

# Number of interior lattice points in x-direction
Nx 10

# Number of interior lattice points in y-direction
Ny 20

# Density of fluid
rho 1.0

# Velocity components (ux uy)
vel 0.1 0.0

# Omega (Viscosity parameter)
omega 1.25

# Number of iterations
iter 1000
