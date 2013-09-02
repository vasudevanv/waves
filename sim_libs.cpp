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

  sim_libs.cpp

  Functions for carrying out various calculations 

*/


// Headers
#include <vector>
#include <cmath>
#include <fftw3.h>
#include "constants.h"
#include "grid.h"
#include "random.h"


void init_velocity ( double *hdot, int nx, int ny, double temp )
{
    /*
      Velocity initialization: Maxwell distribution with rescaling
    */
    int    num;
    double v_sum = 0;
    double v2_sum = 0;
    double sf;
    
    num = nx * ny;
    for ( int i = 0; i < num; i++ )
    {
	hdot[i] = marsaglia();
	v_sum = v_sum + hdot[i];
	v2_sum = v2_sum + hdot[i]*hdot[i];
    }
    v_sum = v_sum / (double) num;    //average velocity
    v2_sum = v2_sum / (double) num;  //average velocity square
    sf = sqrt ( 1.0 * temp / v2_sum );     //scaling factor to desired temperature
    
    //Rescale the velocity to set the initial temperature
    for ( int i = 0; i < num; i++ )
    {
	hdot[i] = ( hdot[i] - v_sum ) *sf;
    } 
}


void dftreal ( double *hfunc, fftw_complex *dfthfunc, 
	       int nx, int ny, int nh ) 
{
    /* 
       Fourier transform of a 2D array using fftw3 library
    */
    fftw_plan  plan_forward;
    
    plan_forward = fftw_plan_dft_r2c_2d ( nx, ny, hfunc, 
					  dfthfunc, FFTW_ESTIMATE );
    fftw_execute ( plan_forward );
    fftw_destroy_plan ( plan_forward );
}


void idftreal ( fftw_complex *dfthfunc, double *hfunc,
		int nx, int ny, int nh ) 
{
    /* 
       Inverse Fourier transform of a 2D array using fftw3 library
    */
    fftw_plan plan_backward;
    
    plan_backward =  fftw_plan_dft_c2r_2d ( nx, ny, dfthfunc,
					    hfunc, FFTW_ESTIMATE );
    fftw_execute ( plan_backward );
    fftw_destroy_plan( plan_backward );
}

void normalize_fft_array ( fftw_complex *h, int nx, int ny, int nyh )
{
    /* 
       Normalize a fftw_complex array using the symmetric normalization
       condition
    */
    int  kx, ky;
    int  el, cDIM;
    
    for ( kx = 0; kx < nx; kx++ )
    {
	for ( ky = 0; ky < nyh; ky++ )
	{
	    el = kx*nyh + ky;
	    for ( cDIM = 0; cDIM < 2; cDIM++ )
	    {
		h[el][cDIM] = h[el][cDIM] / sqrt( (double) (nx * ny) );
	    }
	}
    }
}


void normalize_ifft_array ( double *h, int nx, int ny )
{
    /* 
       Normalize a double array returned by inverse Fourier transform
       using the symmetric normalization condition
    */
    int kx;

    for ( kx = 0; kx < nx * ny; kx++ )
    {
	h[kx] = h[kx] / sqrt( (double) (nx * ny) );
    }
}


void duplicate_fftw_array ( fftw_complex *h, fftw_complex* hc, int nx, int ny )
{
    /* 
       Duplicate a fftw_complex datatype
    */
    int kx, ky;
    int el, cDIM;
    
    for ( kx = 0; kx < nx; kx++ )
    {
	for ( ky = 0; ky < ny; ky++ )
	{
	    el = kx*ny + ky;
	    for ( cDIM = 0; cDIM < 2; cDIM++ )
	    {
		hc[el][cDIM] = h[el][cDIM];
	    }
	}
    }
}


void generate_psitwiddle(fftw_complex *psitwiddle, int nx, int ny, int nh, 
		 double mean, double var)
{
    /*
      Generate the Langevin noise term in fourier space
    */
    int                   i,kx,ky,n,p;
    std::vector<double>   x,y;
    int                   halfmode;
    double                sqrtvar;
    double                twosqrtvar;

    sqrtvar = sqrt ( var );
    twosqrtvar = sqrt ( 2.0 * var );
    n = nx * nh;
    
    x.reserve(n);
    y.reserve(n);
    
    // Generate the random array
    marsaglia_array(x,n);
    marsaglia_array(y,n);
    
    // Need to check for the pure real modes
    if ( nx % 2 == 0 ) 
    {
	halfmode = nx/2 - 1;
    }
    else
    {
	halfmode = (int) ( nx/2 );
    }
    for ( kx = 0; kx < nx; kx++ ) 
    {
	for ( ky = 0; ky < nh; ky++ ) 
	{
	    int el = kx*nh + ky;
	    if ( el == 0 )  
	    {
		x[el] = twosqrtvar*x[el] + mean;
		y[el] = 0.0;
	    }
	    else
	    {
		x[el] = sqrtvar*x[el] + mean;
		y[el] = sqrtvar*y[el] + mean;
	    }
	    psitwiddle[el][0] = x[el];
	    psitwiddle[el][1] = y[el];
	}
    }
    x.clear();
    y.clear();
} 


void calcforce(fftw_complex *htwiddle, fftw_complex *psitwiddle,
	       fftw_complex *fpintwiddle, fftw_complex* Gt, 
	       double *modksq, double gamma, double dx, double dy, 
	       int elstart, int nfftwpoints)
{
    /* 
       Calculate the force due a particular configuration specified
       by htwiddle. Add the contribution of Langevin damping terms
       to the force
    */
    int     el;
    double  gd;
    
    gd = gamma * dx * dy;
    for ( el = elstart; el < nfftwpoints; el++ )
    {
	for ( int cDIM = 0; cDIM < 2; cDIM++ )
	{
	    Gt[el][cDIM] = psitwiddle[el][cDIM] -
		gd * modksq[el] * htwiddle[el][cDIM]
		- fpintwiddle[el][cDIM];
	}
    }
    
    if ( elstart != 0 )
    {
	Gt[0][0] = 0.0;
	Gt[0][1] = 0.0;
    }
}


// void calcforce(fftw_complex *htwiddle, fftw_complex *psitwiddle,
// 	       fftw_complex *fpintwiddle, fftw_complex* Gt, double *modksq, 
// 	       double gamma, double dx, double dy, int nx, int ny, int nh)
// {
//     /* 
//        Calculate the force due a particular configuration specified
//        by htwiddle. Add the contribution of Langevin damping terms
//        to the force
//     */
//     // int     kx,ky;
//     int el;
//     double  gd;
//     
//     gd = gamma * dx * dy;
//     
//     for ( el = elstart; el < nfftwpoints; el++ )
//     {
// 	for ( int cDIM = 0; cDIM < 2; cDIM++ )
// 	{
// 	    Gt[el][cDIM] = psitwiddle[el][cDIM] -
// 		gd * modksq[el] * htwiddle[el][cDIM]
// 		- fpintwiddle[el][cDIM];
// 	}
//     }
// 
//     // for ( kx = 0; kx < nx; kx++ )
//     // {
//     // 	for ( ky = 0; ky < nh; ky++ ) 
//     // 	{
//     // 	    int el = kx * nh + ky;
//     // 	    for ( int cDIM = 0; cDIM < 2; cDIM++ )
//     // 	    {
//     // 		Gt[el][cDIM] = psitwiddle[el][cDIM] -
//     // 		    gd * modksq[el] * htwiddle[el][cDIM]
//     // 		    - fpintwiddle[el][cDIM];
//     // 	    }
//     // 	}
//     // }
// }


void generate_modksq ( double *modksq, t_Grid* grid )
{
    /* 
       Generate the modulus-squared of the k space wave-vectors
    */
    int    kx, ky;
    double tplx, tply;
    
    tplx = TWOPI / grid->LX;
    tply = TWOPI / grid->LY;
    
    for ( kx = 0; kx < grid->nx; kx++ )
    {
	for ( ky = 0; ky < grid->nyh; ky++ ) 
	{ 
	    int el = kx*grid->nyh+ky;
	    if ( grid->nx % 2 == 0 ) 
	    {
		if ( kx > grid->nx/2 ) 
		{
		    modksq[el] = ( (tplx*(kx-grid->nx)) * 
				   (tplx*(kx-grid->nx)) + 
				   (tply*ky) * (tply*ky) );
		}
		else
		{
		    modksq[el] = ( (tplx*kx) * (tplx*kx) + 
				   (tply*ky) * (tply*ky) );
		}
	    }
	    else
	    {
		if ( kx > (grid->nx-1)/2 ) 
		{
		    modksq[el] = ( (tplx*(kx-grid->nx)) * 
				   (tplx*(kx-grid->nx)) +
				   (tply*ky) * (tply*ky) );
		}
		else
		{
		    modksq[el] = ( (tplx*kx) * (tplx*kx) + 
				   (tply*ky) * (tply*ky) );
		}
	    }
	}
    }
} 


double calculate_com ( double *h, int nx, int ny )
{
    /*
      Calculate the center of mass of a given configuration
    */
    double zinit = 0;
    int    num;
    num = nx * ny;
    for ( int i = 0; i < num; i++ )
    {
	zinit = zinit + h[i];
    }
    zinit = zinit/(double)(num);    //center of mass
    
    return zinit;
}


void recenter_comm ( double *h, int nx, int ny, double zinit )
{
    /*
      Recenter the configuration to a given z location
    */
    int     num;
    double  sum = 0;
    
    sum = calculate_com ( h, nx, ny );
    
    // Rescale the coordinates 
    for ( int i = 0; i < num; i++ )
    {
	h[i] = h[i] - sum + zinit;
    }
}

