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

  constants.h
  
  Defines the constants and conversion factors needed for internal units.
  The internal units for the various dimensions are 

  LENGTH       : nm
  MASS         : amu
  TIME	       : ps
  TEMPERATURE  : K

  All other units for quantities such as energy, velocity, acceleration etc. 
  are written in these units. Have a look at .cpp for the units used in the
  input parameter file.
  
*/

#ifndef _constants_h
#define _constants_h

// Defines
#define XX 0
#define YY 1
#define ZZ 2


// Coversion factors
#define MJPM2TOAMUPPS2 0.60221415       // Conversion factor for converting
                                        // from mJ/m**2 -> amu/ps**2

// Constants
const double k_b = 0.0083144621;         // will give kbT in amu nm**2/ps**2
const double PI = 3.14159265359;
const double TWOPI = 2.0*3.14159265359;

#endif /* _constants_h */
