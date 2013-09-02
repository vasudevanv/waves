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

  sim_run.cpp

  Routine which propagates the equations of motion using a velocity-Verlet
  integrator.

*/

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

#include "constants.h"
#include "main.h"
#include "sim_io.h"
#include "sim_libs.h"
#include "grid.h"
#include "sim_field.h"

#define new_fftw(m) (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * m )   

void run_simulation ( inputrec* ir )
{
    // Program variables
    int           nx;        // Number of x-grid points
    int           ny;        // Number of y-grid points 
    int           nyh;       // Number of y-grid points in fft
    XDRFILE*      cnutraj;   // File id for output
    float         xtctime;   // Output file time
    int           step;      // Step number
    double        beta;      // 1/(k_B T)
    double        eta;       // Langevin damping parameter
    int           kx,ky;     // Counter variables;
    double        mean;      // Mean of gaussian random noise
    double        var;       // Variance of gaussian random noise
    matrix        xtcbox;    // Box needed for xtc
    double        m;         // mass of elemental area
    double        cutoffk;   // Cutoff for k vector
    double        zcominit;  // initial z center of mass
    std::ofstream logfile;
    std::ofstream gridfile;

    // Coefficients to simplify computation 
    double        vxcoeff;   
    double        fxcoeff;
    double        vvcoeff;
    double        fvcoeff;
    double        vvdenom;

    // Allocatables
    
    double*        h;                // h(x,y) defines the interface 
    double*        hzref;            // h(x,y) defines the interface 
    double*        hdot;             // Fictitious velocities of h(x,y)
    double*        modksq;           // Modulus square of the wave-vector
    fftw_complex*  htwiddle;         // Fourier transform of h(x,y)
    fftw_complex*  htwiddlecopy;     // Fourier transform of h(x,y)
    fftw_complex*  hdottwiddle;      // Fourier transform of hdot(x,y)
    fftw_complex*  hdottwiddlecopy;  // Fourier transform of hdot(x,y)
    fftw_complex*  psitwiddle;       // Random noise term
    fftw_complex*  Gt;               // Force vector at time t 
    fftw_complex*  Gtdt;             // Force vector at time t+dt
    
    // Terms to deal with external forces
    double*        fpin;
    fftw_complex*  fpintwiddle;       
    t_Field        fd;
   
    // Diffusion
    int elstart;

    // Open the log file 
    if ( ir->nsteps > 0 )
    {
	logfile.open ( ir->log.c_str() );
	gridfile.open ( ir->grf.c_str() );
    }
    else
    {
	std::cout << "Nothing to run\n.";
	std::cout << "If you have not specified the number of steps \n"
		  << "in the parameter file the simulation will not run.\n";
	exit (EXIT_SUCCESS);
    }

    // Print the grid information and close grid_file
    if ( gridfile.is_open() )
    {
	gridfile << ir->grid.box[XX] << "   "
		 << ir->grid.box[YY] << "   " 
		 << ir->grid.box[ZZ] << "\n";
	gridfile << ir->grid.dx << "   " <<  ir->grid.dy << "\n";
	gridfile << ir->grid.nx << "   " <<  ir->grid.ny << "\n";
	gridfile.close();
    }

    if ( ir->bdiffusion )
	elstart = 0;
    else
	elstart = 1;

    std::cout << "elstart = " << elstart << "\n";
    
    // Declare some variables based on input
    nx = ir->grid.nx;
    ny = ir->grid.ny;
    nyh = ir->grid.nyh; 

    beta = 1.0 / ( k_b * ir->T );
    m = ir->mu * ir->grid.dx * ir->grid.dy;
    eta = m / ir->tau;     
    cutoffk = ( TWOPI / ir->l ) * ( TWOPI / ir->l );
    mean = 0.0;
    var = eta * k_b * ir->T / ir->dt;
    vxcoeff = ir->dt * ( 1 - ir->dt / ( 2.0 * ir->tau ) );
    fxcoeff = 0.5* ir->dt *ir->dt / m ;
    vvcoeff = 1 - ir->dt / ( 2.0 * ir->tau ) ;
    fvcoeff = 0.5 * ir->dt / m ;
    vvdenom = 1 + ir->dt / ( 2.0 * ir->tau ) ;

    // Write box information in a form that will be used by the xtc file
    for ( int cDIM = 0; cDIM < 3; cDIM++ )
    {
	for ( int dDIM = 0; dDIM < 3; dDIM++ )
	{
	    xtcbox[cDIM][dDIM] = 0;
	}
    }
    xtcbox[0][0] = ir->grid.box[XX];
    xtcbox[1][1] = ir->grid.box[YY];
    xtcbox[2][2] = ir->grid.box[ZZ];

    // Write stuff to log file
    print_input_parameters ( ir );
    printlog_input_parameters ( ir, logfile );

    // Allocate arrays
    int ngridpoints = ir->grid.nx * ir->grid.ny;
    int nfftwpoints = ir->grid.nx * ir->grid.nyh;

    h      = new double[ngridpoints];
    hdot   = new double[ngridpoints];
    fpin   = new double[ngridpoints];
    modksq = new double[nfftwpoints];
    htwiddle        = new_fftw ( nfftwpoints );
    hdottwiddle     = new_fftw ( nfftwpoints );
    fpintwiddle     = new_fftw ( nfftwpoints );
    htwiddlecopy    = new_fftw ( nfftwpoints );
    hdottwiddlecopy = new_fftw ( nfftwpoints );
    psitwiddle      = new_fftw ( nfftwpoints );
    Gt              = new_fftw ( nfftwpoints );
    Gtdt            = new_fftw ( nfftwpoints );
           
    // Initialize the coordinates and velocities
    srand((unsigned)time(NULL));
    for ( int i = 0; i < ngridpoints; i++ )
    {
        h[i] = ir->hinit;
	hdot[i] = 0.0;
    }
    //init_velocity(hdot,nx,ny,ir->T);

    // Open output file
    if ( ir->nsteps > 0 )
    {
	write_gro ( ir->gro.c_str(), h, &ir->grid );
	cnutraj = xdrfile_open(ir->traj.c_str(),"w");
	writeframe(h, &ir->grid,cnutraj,0,0.0,xtcbox);
    }
    
    // Preprocess some variables before entering loop
    // Calculate force due to constraint hamiltonian here 
    generate_modksq ( modksq, &ir->grid ); 
    
    if ( ir->bfield )
    {
        init_fields(ir->fieldfile.c_str(), &fd, ir->T, logfile );
	field_potential ( h, fpin, &ir->grid, &fd );
	dftreal ( fpin, fpintwiddle, nx, ny, nyh );
	normalize_fft_array ( fpintwiddle, nx, ny, nyh );
    }
    else
    {
	// Set the force array to zero
	for ( int el = 0; el < nfftwpoints; el++ )
	    for ( int cDIM = 0; cDIM < 2; cDIM++ )
		fpintwiddle[el][cDIM] = 0.0; 
    }
    dftreal ( h, htwiddle, nx, ny, nyh );
    dftreal ( hdot, hdottwiddle, nx, ny, nyh );
    normalize_fft_array ( htwiddle, nx, ny, nyh );
    normalize_fft_array ( hdottwiddle, nx, ny, nyh );    
    generate_psitwiddle ( psitwiddle, nx, ny, nyh, mean, var );
    calcforce(htwiddle,psitwiddle,fpintwiddle,Gt,modksq,ir->gamma,
	      ir->grid.dx,ir->grid.dy,elstart,nfftwpoints);    
    
    // Main loop
    std::cout << "Starting the integration\n";
    for ( step = 1; step <= ir->nsteps; step++ ) 
    {
	xtctime = step*ir->dt;

	if ( step % 1000 == 0 )
	{
	    std::cout << "\r" << "Finished " << step << " of " << ir->nsteps << " steps       " << std::flush;
	}
	// Generate the coordinates
	for ( int el = elstart; el < nfftwpoints; el++ )
	{
	    if ( modksq[el] < cutoffk )
	    {
		for ( int cDIM = 0; cDIM < 2; cDIM++ ) 
		{
		    htwiddle[el][cDIM] += 
			vxcoeff * hdottwiddle[el][cDIM] + 
			fxcoeff * Gt[el][cDIM];
		}
	    }
	    else
	    {
		for ( int cDIM = 0; cDIM < 2; cDIM++ ) 
		    htwiddle[el][cDIM] = 0.0;
	    }
	}
		
	// Duplicate arrays since the inverse transform destroys the original arrays
	duplicate_fftw_array ( htwiddle, htwiddlecopy, nx, nyh );
	idftreal( htwiddlecopy, h, nx, ny, nyh );
	normalize_ifft_array ( h, nx, ny );

	// Evaluate force due to constrain hamiltonian terms in real space
	if ( ir->bfield )
	{	
	    field_potential ( h, fpin, &ir->grid, &fd );
	    dftreal ( fpin, fpintwiddle, nx, ny, nyh );
	    normalize_fft_array ( fpintwiddle, nx, ny, nyh );
	}

	// Generate psitwiddle(t+dt)
	generate_psitwiddle(psitwiddle,nx,ny,nyh,mean,var);
	
	// Compute the force in the current configuration
	calcforce( htwiddle,psitwiddle,fpintwiddle,Gtdt,modksq,ir->gamma,
		   ir->grid.dx,ir->grid.dy,elstart,nfftwpoints);

	// Update the velocities
	for ( int el = elstart; el < nfftwpoints; el++ ) 
        {         
	    if ( modksq[el] < cutoffk )
	    {
		for ( int cDIM = 0; cDIM < 2; cDIM++ ) 
		{ 
		    hdottwiddle[el][cDIM] = 
			vvcoeff * hdottwiddle[el][cDIM] + 
			fvcoeff * ( Gt[el][cDIM] + Gtdt[el][cDIM] );
		    
		    hdottwiddle[el][cDIM] = 
			hdottwiddle[el][cDIM] / vvdenom;
		}
		
	    }
	    else
	    {
		for ( int cDIM = 0; cDIM < 2; cDIM++ ) 
		    hdottwiddle[el][cDIM] = 0.0;
	    }
	}
        	
	// Duplicate arrays since the inverse transform destroys the original arrays
	duplicate_fftw_array ( hdottwiddle, hdottwiddlecopy, nx, nyh );
	idftreal( hdottwiddlecopy, hdot, nx, ny, nyh );
	normalize_ifft_array ( hdot, nx, ny );
	
	// Write output
	if ( step % ir->nxtc == 0 )
	    writeframe(h,&ir->grid,cnutraj,step,xtctime,xtcbox);
        
        // Initialize stuff for next step
	duplicate_fftw_array ( Gtdt, Gt, nx, nyh );
    }
    
    // Deallocate arrays and close files
    delete [] h;
    delete [] hdot;
    delete [] fpin;
    delete [] modksq;
    fftw_free ( htwiddle );
    fftw_free ( hdottwiddle );
    fftw_free ( psitwiddle ); 
    fftw_free ( htwiddlecopy );
    fftw_free ( hdottwiddlecopy );
    fftw_free ( Gt );
    fftw_free ( Gtdt );
    fftw_free ( fpintwiddle );
   
    if ( ir->bfield )
	finish_fields( &fd );
 
    if ( cnutraj )
	xdrfile_close(cnutraj);
    
    if ( logfile.is_open() )
	logfile.close();
    
    std::cout << "\n";

}


