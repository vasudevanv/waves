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

  inputrec.h

  Class definition for the input record. The input record contains all
  the input variables/parameters needed for the simulation. 

*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <map>
#include <xdrfile.h>
#include <xdrfile_xtc.h>

#include "constants.h"
#include "grid.h"

#ifndef _sim_io_h
#define _sim_io_h

class inputrec
{
    
private:
    int           nsteps;       // Total number of steps
    double        gamma;        // Surface tension
    t_Grid        grid;         // The grid representing the surface
    double        T;            // Temperature 
    double        mu;           // Mass per unit area of membrane    
    double        l;            // Cut-off length for wave vectors
    double        dt;           // Timestep  
    double        tau;          // Timescale for momentum decorrelation
    double        hinit;        // Starting uniform height
    int           nxtc;         // Trajectory writing frequency
    std::string   gro;          // Output file name for initial coordinates
    std::string   traj;         // Output file name for trajectory
    std::string   log;          // Output file name for log
    std::string   grf;          // Output file name for grid information file
    std::string   fieldfile;    // File containing harmonic pin parameters
    bool          bhinit;       // Boolean to check for hinit
    bool          bfield;       // Boolean to check for external field
    bool          bgrid;        // Boolean to check if grid was given 
    bool          bdiffusion;   // Allow diffusion of membrane 

public:
    inputrec ( );
    void get_internal_units ( );
    friend int main ( int argc, char *argv[] ); 
    friend void parse_parameter_file ( const char*, inputrec* );
    friend void run_simulation ( inputrec* );
    friend void print_input_parameters ( inputrec* inp );
    friend void printlog_input_parameters ( inputrec* inp, std::ofstream& );
};


enum ParameterValue 
{
    /*
      enum defining the input parameters
    */
    evNotDefined,
    evT,
    evdt,
    evNsteps,
    evGamma,
    evTau,
    evTraj,
    evGro,
    evLog,
    evGrf,
    evMu,
    evl,
    evGrid,
    evBox,
    evHinit,
    evField,
    evNxtc,
    evDiffusion
};

static void initialize_map ( );

char *fgets2(char *line, int n, FILE *stream);

void ltrim (char *str);

void rtrim (char *str);

void trim (char *str);

std::string trim_right_copy ( const std::string& s,
			      const std::string& delimiters = " \f\n\r\t\v" );

std::string trim_left_copy ( const std::string& s, 
			     const std::string& delimiters = " \f\n\r\t\v" );

std::string trim_copy( const std::string& s,
		       const std::string& delimiters = " \f\n\r\t\v" );

void parse_parameter_file ( const char *fn, inputrec* inp );

void print_usage ( void );

void write_gro ( const char *gro, double* h, t_Grid* grid );

int writeframe(double *h, t_Grid* grid, XDRFILE* xd, int step, 
	       float xtctime, matrix box);

void print_input_parameters ( inputrec* inp );

void printlog_input_parameters ( inputrec* inp, std::ofstream& fp );


#endif
