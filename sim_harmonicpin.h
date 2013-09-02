/*

  Copyright 2013, Vasudevan Venkateshwaran,  Garde group @ RPI
  
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

  Header file containing the class which represents a harmonically pinned 
  point on the membrane. 

  U_pin = (K/2.0) * (h(r)**2) * exp (-((r-R)/(0.25*l))**2)
  
  R = Pin location
  
  l = gridsize L/N.

*/


#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include "grid.h"
#include "constants.h"

#ifndef _harmonicpin_h
#define _harmonicpin_h

#include "sim_field.h"

class harmonicpin: public Field
{
private: 
    
    double x,y,z;
    double k;
    double l;
    
public:
    harmonicpin();
    virtual ~harmonicpin();
    void init ( int argc, char **argv, std::ofstream& fp );
    void calc_forces ( double* h, double *force, t_Grid* grid );      
};

#endif /* _harmonicpin_h */
