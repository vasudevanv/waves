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

  grid.cpp

  Defines the member functions of the t_Grid class and other routines need 
  for manipulating the grid

*/


#include "constants.h"
#include "grid.h"


t_Grid :: t_Grid ( )
{
    /* 
       Constructor: sets the default values
    */
    dx      = 0.1;     
    dy      = 0.1;     
    box[XX] = 6.0;    
    box[YY] = 6.0;    
    box[ZZ] = 5.0;    
    LX      = box[XX];
    LY      = box[YY];    
    nx      = (int) ( box[XX] / dx );
    ny      = (int) ( box[YY] / dy );
    nyh     = ( ny / 2 ) + 1;
}


t_Grid& t_Grid :: operator = ( const t_Grid& param )
{
    /*
      = operator
    */
    dx      = param.dx;
    dy      = param.dy;
    box[XX] = param.box[XX];
    box[YY] = param.box[YY];
    box[2]  = param.box[2];
    LX      = param.LX;
    LY      = param.LY;    
    nx      = param.nx;
    ny      = param.ny;
    nyh     = param.nyh;
    
    return *this;
}

void  t_Grid :: grid_reinit ( )
{
    /*
      Recalculate some of the variables after user input
    */
    nx  = ( box[XX] / dx );
    ny  = ( box[YY] / dy );
    nyh = ( ny / 2 ) + 1;

    if ( box[XX] - nx*dx < dx && box[XX] - nx*dx > 0 )
    {
	nx = nx + 1;
	box[XX] = nx * dx;
    }

    if ( box[YY] - ny*dy < dy  && box[YY] - ny*dy > 0 )
    {
	ny = ny + 1;
	box[YY] = ny * dy;
    }
    
    LX = box[XX];
    LY = box[YY];
}
