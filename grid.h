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

  grid.h

  Defines a class which represents the grid.

*/


#ifndef _grid_h
#define _grid_h


class t_Grid
{
    /* 
       Grid class : defines variables and routines to deal with the grid
    */
public:
    double dx;
    double dy;
    double box[3];
    double LX, LY;
    int    nx;
    int    ny;
    int    nyh;
    
    t_Grid();
    t_Grid& operator = ( const t_Grid& param );
    void  grid_reinit ( );
};

#endif /* _grid_h */
