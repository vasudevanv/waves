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

  sim_libs.h

  The main library routine definitions needed by the program

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

#ifndef _sim_libs_h
#define _sim_libs_h

/* Routines in sim_libs.cpp */
// Center of mass
double calculate_com ( double *h, int nx, int ny );

void remove_comm ( double *h, int nx, int ny, double zinit );


// FFTW FFT calculation
void duplicate_fftw_array ( fftw_complex *h, fftw_complex* hc, int nx, int ny );

void dftreal ( double *h, fftw_complex *htwiddle, int nx, int ny, int nh );

void idftreal ( fftw_complex *htwiddle, double *h, int nx, int ny, int nh ); 

void normalize_fft_array ( fftw_complex *h, int nx, int ny, int nyh );

void normalize_ifft_array ( double *h, int nx, int ny );


// Initalizing the velocity
void init_velocity ( double* , int , int , double );


// Generating noise terms for langevin damping
void generate_psitwiddle(fftw_complex *psitwiddle, int nx, int ny, int nh, 
		 double mean, double var);


// Generating the k-vector magnitudes
void generate_modksq ( double *modksq, t_Grid* grid );


// Force calculation
void calcforce(fftw_complex *htwiddle, fftw_complex *psitwiddle,
	       fftw_complex *fpintwiddle, fftw_complex *Gt, double *modksq, 
	       double gamma, double dx, double dy, int elstart, int nfftwpoints );

#endif /* _sim_libs_h */
