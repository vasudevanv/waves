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
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <vector>
#include "grid.h"
#include "constants.h"

#ifndef _field_h
#define _field_h

typedef struct 
{
    void *sFieldList;
    int num_fields;
} t_Field;


class Field
{

public:
    Field();
    virtual ~Field();
    virtual void init( int argc, char **argv, std::ofstream& fp );
    virtual void calc_forces( double* h, double *force, t_Grid* grid );

protected:
    double T;  // temperature in kelvin
    
    // Friend functions
    friend void init_fields ( const char* fn, t_Field* fd, double T,
			      std::ofstream& fp );
};

typedef struct
{
        int argc;
        char **argv;
        char *line;

} FieldArgs;


void init_fields( const char* fn, t_Field* fd, double T, std::ofstream& fp );

void finish_fields( t_Field* fd );

void field_potential( double* h, double* f, t_Grid* grid, t_Field* fd );

#endif /* _field_h */



