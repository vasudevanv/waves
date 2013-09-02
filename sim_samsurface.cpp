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
#include "sim_samsurface.h"


samsurface:: samsurface ()
{
    
}

samsurface:: ~samsurface ()
{

}

void samsurface:: init ( int argc, char **argv, std::ofstream& fp ) 
{
    if ( argc <  5 )
    {
	std::cout << "Incorrect parameters for 9-3 surface.\n"; 
	std::cout << "Format: samsurface zlocation sigma epsilon rho "
		  << "(R_0 = 0.337 nm) (rho_l = 33.33 nm**-3 ) \n";
	std::cout << "The parameters in paranthesis are optional. \n";
	exit (EXIT_FAILURE);
    }
    
    zref         = strtod(argv[1], NULL);
    sigma_head   = strtod(argv[2], NULL);
    epsilon_head = strtod(argv[3], NULL);
    rho_head     = strtod(argv[4], NULL);
    sigma_tail   = strtod(argv[5], NULL);
    epsilon_tail = strtod(argv[6], NULL);
    rho_tail     = strtod(argv[7], NULL);
    
    switch ( argc )
    {
    case 8:
	lambda = 1.0;
	xi = 0.4525; // nm
	R_0 = 0.337; // nm
	rho_l = 33.3333; // #/nm^3
	break;
    case 9:
	lambda = strtod(argv[8], NULL);
	xi = 0.4525; // nm
	R_0 = 0.337; // nm
	rho_l = 33.3333; // #/nm^3
	break;
    case 10: 
	lambda = strtod(argv[8], NULL);
	xi = strtod(argv[8], NULL);
	R_0 = 0.337; 
	rho_l = 33.3333; 
	break;
    case 11: 
	lambda = strtod(argv[8], NULL);
	xi = strtod(argv[8], NULL);
	R_0 = strtod(argv[5], NULL); 
	rho_l = 33.3333; 
	break;
    case 12:
	lambda = strtod(argv[8], NULL);
	xi = strtod(argv[8], NULL);
	R_0 = strtod(argv[5], NULL); 
	rho_l = strtod(argv[6], NULL); 
	break;
    default:
	std::cout << "Incorrect parameters for sam surface. \n";
	std::cout << "Format: samsurface zlocation sigma_head epsilon_head rho_head "
		  << "sigma_tail epsilon_tail rho_tail xi (lambda = 1.00)"
		  << "(R_0 = 0.337 nm) (rho_l = 33.33 nm**-3 ) \n";
	std::cout << "The parameters in paranthesis are optional. \n";
	exit (EXIT_FAILURE);
    }
    d_skin = 0.01; // skin depth in nm
    sigmasq_head = sigma_head*sigma_head;
    sigmacu_tail = sigma_tail*sigma_tail*sigma_tail;

    std::cout << "Initializing sam surface with parameters : \n";
    std::cout << "zref_surface    = " << zref         << "\n"
	      << "sigma_head      = " << sigma_head   << "\n"
	      << "epsilon_head    = " << epsilon_head << "\n"
	      << "rho_head        = " << rho_head     << "\n"
	      << "sigma_head      = " << sigma_tail   << "\n"
	      << "epsilon_head    = " << epsilon_tail << "\n"
	      << "rho_head        = " << rho_tail     << "\n"
	      << "lambda          = " << lambda       << "\n"
	      << "R_0             = " << R_0          << "\n"
	      << "rho_liquid      = " << rho_l        << "\n"; 

    fp << "Initializing sam surface with parameters : \n";
    fp << "zref_surface    = " << zref         << "\n"
       << "sigma_head      = " << sigma_head   << "\n"
       << "epsilon_head    = " << epsilon_head << "\n"
       << "rho_head        = " << rho_head     << "\n"
       << "sigma_head      = " << sigma_tail   << "\n"
       << "epsilon_head    = " << epsilon_tail << "\n"
       << "rho_head        = " << rho_tail     << "\n"
       << "lambda          = " << lambda       << "\n"
       << "R_0             = " << R_0          << "\n"
       << "rho_liquid      = " << rho_l        << "\n"; 
}

void samsurface:: calc_forces ( double* h, double *force, t_Grid* grid )
{
    /*
      Calculate the force due to a sam surface potential
    */
    int     kx,ky,el;
    double  rx,ry,rz;
    double  dsq;
    double  rzsq, ratiosq_head;
    double  ratiocu_tail;
    double  prefactor;
    
    prefactor = rho_l * grid->dx * grid->dy;

    for ( kx = 0; kx < grid->nx; kx++ )
    {
	for ( ky = 0; ky < grid->ny; ky++ )
	{
	    el = kx * grid->ny + ky; 
	    rz = h[el] - zref;
	    ratiosq_head = sigmasq_head/(rz*rz);
	    ratiocu_tail = sigmacu_tail/((rz+xi)*(rz+xi)*(rz+xi));
	    if ( rz < 0 )
	    {
		std::cout<< "rz error";
		exit(EXIT_FAILURE);
	    }

	    if ( rz > R_0 )
	    {
		force[el] += lambda*prefactor*(TWOPI)*rho_head*epsilon_head*sigmasq_head*
		    ratiosq_head*ratiosq_head* (1 - (2.0/5.0)*(ratiosq_head*ratiosq_head*ratiosq_head))
		    + (TWOPI/3.0)*rho_tail*epsilon_tail*sigmacu_tail*
		    (ratiocu_tail)*(1 - (2.0/15.0)*ratiocu_tail*ratiocu_tail) ; 
	    }
	    else
	    {
		force[el] += prefactor*(-2.0/rho_l)* k_b * T *
		    (R_0 - rz)/(d_skin*d_skin);
		
	    }
	}
    }
}



