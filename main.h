/*

  Copyright 2013, Vasudevan Venkateshwaran
  
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

  main.h

  Header file containing includes and definitions needed by all other 
  cpp files.

*/

// Headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <vector>
#include <ctime>
#include "unistd.h"
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "fftw3.h"
#include "sim_io.h"
#include "grid.h"

#ifndef _main_h
#define _main_h

/* Routines in sim_run.cpp */
// Run the simulation
void run_simulation ( inputrec* ir );

#endif /* _main_h */

