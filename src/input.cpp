
#include "input.h"

/*!
 *	Default constructor for input object
*/
input::input(){

    // Default values for the input parameters
    title = "LBM";
    lx    = 1.0;
    ly    = 1.0;
    Nx    = 10;
    Ny    = 10;
    den   = 1.0;
    ux    = 1.0;
    uy    = 0.0;
    omega = 1.0;
    nIter = 1;
}

/*!
 *	Reads Input file
*/
void input::readInputFile(){

    string lineString;
    string dummyString;
    
    ifstream inputFile;
    inputFile.open("input",ios::in);

    if (inputFile.is_open()==false){
        cout << "Unable to open input file! Aborting... " << endl;
        exit(0);
    }

    cout.precision(7);
    cout << scientific;
    cout << endl << "====== Inputs ======" << endl;
    
    while (!inputFile.eof())
    {
        // Get a line and store in lineString
        getline(inputFile, lineString, '\n');

        // If the first character of the line is not a '#'
        if (lineString.c_str()[0] != '#')
        {
            istringstream iss(lineString);
            iss >> dummyString;
            if(dummyString == "title")
                iss >> title;
            else if(dummyString == "wdir")
                iss >> wdir;
            else if(dummyString == "lx")
                iss >> lx;
            else if(dummyString == "ly")
                iss >> ly;
            else if(dummyString == "Nx")
                iss >> Nx;
            else if(dummyString == "Ny")
                iss >> Ny;
            else if(dummyString == "rho")
                iss >> den;
            else if(dummyString == "vel"){
                iss >> ux;
		iss >> uy;
            }else if(dummyString == "omega")
                iss >> omega;
            else if(dummyString == "iter")
                iss >> nIter;
            else{
                cout << endl << "Unknown keyword in the settings file : " << dummyString;
                cout << endl << "Aborting...";
                exit(0);    
            }
        }
        
    }

    // Report the settings read from the file.

    //cout << fixed;
    cout << "Title of the simulation                 : " << title << endl;
    cout << "Working Directory                       : " << wdir << endl;
    cout << "Length of physical domain               : " << lx << endl;
    cout << "Width of physical domain                : " << ly << endl;
    cout << "No. of interior lattice points in x-dir : " << Nx << endl;
    cout << "No. of interior lattice points in y-dir : " << Ny << endl;
    cout << "Density of fluid                        : " << den << endl;
    cout << "Velocity components(ux uy)              : " << ux <<"  "<< uy << endl;
    cout << "Viscosity parameter (Omega)             : " << omega << endl;
    cout << "Number of iterations                    : " << nIter  << endl;
    cout << endl << endl;

    inputFile.close();
    
    return;
}


