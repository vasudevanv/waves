/*

  Copyright 2013, Vasudevan Venkateshwaran, Garde group @ RPI
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  sim_harmonicpin.cpp

  Routines to deal with 9-3 surface. 

*/

#include <string>
#include <fstream>
#include <iostream>
#include "constants.h"
#include "sim_surface93.h"


surface93:: surface93 ()
{
    
}

surface93:: ~surface93 ()
{

}

void surface93:: init ( int argc, char **argv, std::ofstream& fp ) 
{
    if ( argc <  5 )
    {
	std::cout << "Incorrect parameters for 9-3 surface.\n"; 
	std::cout << "Format: surface93 zlocation sigma epsilon rho "
		  << "(R_0 = 0.337 nm) (rho_l = 33.33 nm**-3 ) \n";
	std::cout << "The parameters in paranthesis are optional. \n";
	exit (EXIT_FAILURE);
    }
    
    zref    = strtod(argv[1], NULL);
    sigma   = strtod(argv[2], NULL);
    epsilon = strtod(argv[3], NULL);
    rho     = strtod(argv[4], NULL);
    
    switch ( argc )
    {
    case 5:
	R_0 = 0.337; // nm
	rho_l = 33.3333; // #/nm^3
	break;
    case 6: 
	R_0 = strtod(argv[5], NULL); 
	rho_l = 33.3333; 
	break;
    case 7:
	R_0 = strtod(argv[5], NULL); 
	rho_l = strtod(argv[6], NULL); 
	break;
    default:
	std::cout << "Incorrect parameters for 9-3 surface. \n";
	std::cout << "Format: surface93 zlocation sigma epsilon rho "
		  << "(R_0 = 0.337 nm) (rho_l = 33.33 nm**-3 ) \n";
	std::cout << "The parameters in paranthesis are optional. \n";
	exit (EXIT_FAILURE);
    }
    d_skin = 0.01; // skin depth in nm
    sigmasq = sigma*sigma;
    
    std::cout << "Initializing 9-3 surface with parameters : \n";
    std::cout << "zref_surface    = " << zref    << "\n"
	      << "sigma_surface   = " << sigma   << "\n"
	      << "epsilon_surface = " << epsilon << "\n"
	      << "rho_surface     = " << rho     << "\n"
	      << "R_0             = " << R_0     << "\n"
	      << "rho_liquid      = " << rho_l   << "\n"; 

    fp << "Initializing 9-3 surface with parameters : \n";
    fp << "zref_surface    = " << zref    << "\n"
       << "sigma_surface   = " << sigma   << "\n"
       << "epsilon_surface = " << epsilon << "\n"
       << "rho_surface     = " << rho     << "\n"
       << "R_0             = " << R_0     << "\n"
       << "rho_liquid      = " << rho_l   << "\n"; 
}

void surface93:: calc_forces ( double* h, double *force, t_Grid* grid )
{
    /*
      Calculate the force due to a 9-3 surface potential
    */
    int     kx,ky,el;
    double  rx,ry,rz;
    double  dsq;
    double  rzsq, ratio, ratiosq;
    double  prefactor;
    
    prefactor = rho_l * grid->dx * grid->dy;

    for ( kx = 0; kx < grid->nx; kx++ )
    {
	for ( ky = 0; ky < grid->ny; ky++ )
	{
	    el = kx * grid->ny + ky; 
	    rz = h[el] - zref;
	    ratio = (sigma/rz);
	    ratiosq = ratio*ratio;
	    if ( rz < 0 )
	    {
		std::cout<< "rz error";
		exit(EXIT_FAILURE);
	    }

	    if ( rz > R_0 )
	    {
		force[el] += prefactor*(TWOPI)*rho*epsilon*sigmasq*
		    ratiosq*ratiosq*
		    (1 - (2.0/5.0)*(ratiosq*ratiosq*ratiosq)); 
	    }
	    else
	    {
		force[el] += prefactor*(-2.0/rho_l)* k_b * T *
		    (R_0 - rz)/(d_skin*d_skin);
		
	    }
	}
    }
}



