/*

  Copyright 2013, Vasudevan Vekateshwaran, Srivathsan Ranganathan
  
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

  main.cpp

  The main program.

*/

#include "main.h"
#include "sim_io.h"
#include "sim_field.h"


int main ( int argc, char *argv[] )
{
    int       opt,prev_ind,flag;  // Option parser variable 
    char      *parmfile;          // Filename of parameter file
    inputrec   ir;

    // Parse the command line options
    flag = 0;
    if ( argc < 2 ) 
    {
	print_usage ( ) ;
	exit ( EXIT_FAILURE );
    }
    while ( prev_ind = optind, ( opt = getopt ( argc, argv, "f:h" ) ) != -1 )
    {
	if ( optind == prev_ind + 2 && *optarg == '-' ) 
	{
	    opt = ':';
	    -- optind;
	}
	switch (opt)
	{
	case 'f': 
	    if ( optarg[0] != '-' )
	    {
		parmfile = optarg;
		flag = 1;
	    }
	    break;
	case 'h':
	case ':':
	case '?': 
	    print_usage ( );
	    exit(EXIT_FAILURE);
	}
    }

    if ( ! flag )
    {
	print_usage ( );
	exit(EXIT_FAILURE);
    }

    parse_parameter_file ( parmfile, &ir );

    run_simulation ( &ir );

    return 0;
}


