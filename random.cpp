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

  random.cpp

  Routines to generate random deviates

*/
  
  
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "random.h"


// Uniform random variables
float ran(int iseed)
{
    /* Returns a uniform random deviate between 0.0 and 1.0. ;
       Based on: Park and Miller's "Minimal Standard" random number;
       generator (Comm. ACM, 31, 1192, 1988);
    */
    int  IM=2147483647, IA=16807, IQ=127773, IR= 2836;
    float AM=128.0/IM;
    int K;
    
    K = iseed/IQ;
    iseed = IA*(iseed-K*IQ) - IR*K;
    if (iseed < 0) iseed = iseed+IM;
    return AM*(iseed/128);
}


// Gaussian random variables
double marsaglia()
{
    /* 
       Function to return a zero mean, unit variance normally 
       distributed deviate using the Marsaglia polar method
    */
    double x, y;
    double r;
    int   flag;
    
    // Generate uniform random points inside a unit circle;
    do 
    {
	x = ((double) rand() / ((double)RAND_MAX+1)) ;
	y = ((double) rand() / ((double)RAND_MAX+1)) ;
	x = 2.0*x - 1.0;
	y = 2.0*y - 1.0;
	r = x*x + y*y;
    } while (r > 1.0);
    x = x*sqrt(-2.0*log(r)/r);
    y = y*sqrt(-2.0*log(r)/r);
    return x;
}


void marsaglia_array(std::vector<double> &x, int n)
{
    /* 
       Function to return an array of size n containing zero mean, 
       unit variance normally distributed deviates using the Marsaglia 
       polar method
    */
    double  r,tx,ty;
    int    i;
    
    // Generate uniform random points inside a unit circle;
    for(i = 1;i<=n;i+=2) 
    {
	do 
	{ 
	    tx = ((double) rand() / ((double)RAND_MAX+1)) ;
	    ty = ((double) rand() / ((double)RAND_MAX+1)) ;
	    tx = 2.0*tx - 1.0;
	    ty = 2.0*ty - 1.0;
	    r = tx*tx + ty*ty;
	} while ( r > 1.0 );
	
	tx = tx*sqrt(-2.0*log(r)/r);
	ty = ty*sqrt(-2.0*log(r)/r);
	
	if (i+1 <= n) 
	{
	    x[i-1] = tx;
	    x[i] = ty;
	}
	else 
	{
	    x[i-1] = tx;
	}
    }
}

