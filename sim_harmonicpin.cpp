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

  Routines to deal with harmonically pinned points on the membrane. 

*/

#include <string>
#include <fstream>
#include <iostream>
#include "sim_harmonicpin.h"

harmonicpin:: harmonicpin ()
{
    
}

harmonicpin:: ~harmonicpin ()
{

}

void harmonicpin:: init ( int argc, char **argv, std::ofstream& fp ) 
{
    if ( argc <  6 )
    {
	std::cout << "Incorrect parameters for harmonic pin.\n"; 
	std::cout << "Format: harmonicpin x y z k l \n";
	exit (EXIT_FAILURE);
    }
    
    x  = strtod(argv[1], NULL);
    y  = strtod(argv[2], NULL);
    z  = strtod(argv[3], NULL);
    k  = strtod(argv[4], NULL);
    l  = strtod(argv[5], NULL);
 
    std::cout << "Initializing harmonic pin with parameters : \n";
    std::cout << "X location       = " << x << "\n"
	      << "Y location       = " << y << "\n"
	      << "Z location       = " << z << "\n"
	      << "K (Spring Const) = " << k << "\n"
	      << "l                = " << l << "\n"; 

    fp << "Initializing harmonic pin with parameters : \n";
    fp << "X location       = " << x << "\n"
       << "Y location       = " << y << "\n"
       << "Z location       = " << z << "\n"
       << "K (Spring Const) = " << k << "\n"
       << "l                = " << l << "\n"; 
}

void harmonicpin :: calc_forces ( double* h, double *force, t_Grid* grid )
{
    /*
      Calculate the force due to a harmonic pinning potential
    */
    int     kx,ky,el;
    double  rx,ry,rz;
    double  dsq;

    for ( kx = 0; kx < grid->nx; kx++ )
    {
	for ( ky = 0; ky < grid->ny; ky++ )
	{
	    el = kx * grid->ny + ky; 
	    rx = kx*grid->dx - x;
	    ry = ky*grid->dy - y;
	    if ( rx >  grid->LX/2.0 )  rx = rx - grid->LX;
	    if ( rx < -grid->LX/2.0 )  rx = rx + grid->LX;
	    if ( ry >  grid->LY/2.0 )  ry = ry - grid->LY;
	    if ( ry < -grid->LY/2.0 )  ry = ry + grid->LY;
	    dsq = rx*rx + ry*ry;
	    dsq = dsq / ( 0.25 * 0.25 * l * l );
	    force[el] += k * (h[el] - z) * exp ( -dsq ); 
	}
    }
}
