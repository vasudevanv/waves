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

  random.h

  Header file for routines which perform random number generation

*/


#include <vector>


#ifndef _random_h
#define _random_h


/* Routines in random.cpp */
// Generate uniform random number between 0 and 1
float ran(int iseed);

// Generate standard normal random variables using the 
// Marsaglia polar method
double marsaglia();

void marsaglia_array(std::vector<double> &x, int n);


#endif /* _random_h */
