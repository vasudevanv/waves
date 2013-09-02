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

  sim_harmonicpin.h

  Header file containing the class which represents a uniform field applied 
  along the interface normal.

*/

#include <iostream>
#include <fstream>
#include "grid.h"
#include "constants.h"

#ifndef _samsurface_h
#define _samsurface_h

#include "sim_field.h"

class samsurface: public Field
{
private:
    double zref;
    double sigma_head;    // LJ sigma of head groups (nm)   
    double epsilon_head;  // LJ epsilon of head groups (kJ/mol)
    double rho_head;      // Surface density of head groups (#/nm^2)
    double sigma_tail;    // LJ sigma of tail groups (nm)   
    double epsilon_tail;  // LJ epsilon of tail groups (kJ/mol)
    double rho_tail;      // Surface density of tail groups (#/nm^2)
    double xi;            // Depth of tail groups 
    double lambda;        // Surface scaling parameter
    double R_0;
    double rho_l;
    double d_skin;
    double sigmasq_head;
    double sigmacu_tail;

public:
    samsurface();
    virtual ~samsurface();
    void init ( int argc, char **argv, std::ofstream& fp );
    void calc_forces ( double* h, double *f, t_Grid* grid);
      
};

#endif /* _samsurface_h */
