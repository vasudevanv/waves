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

  sim_ion.cpp

  Functions used to parse the input parameter file 

*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cstring>
#include <map>
#include <algorithm>
#include <xdrfile.h>
#include <xdrfile_xtc.h>

#include "constants.h"
#include "grid.h"
#include "sim_io.h"


static std::map< std::string, ParameterValue > mapParameterValues;


static void initialize_map ( )
{
    /* This routine maps some string values to the enum values.
       The enum values will be used later to ascertain which 
       inputs were specified in the parameter file
    */
    mapParameterValues["temperature"] = evT;
    mapParameterValues["nsteps"] = evNsteps;
    mapParameterValues["gamma"] = evGamma;
    mapParameterValues["tau"] = evTau;
    mapParameterValues["dt"] = evdt;
    mapParameterValues["trajectory_file"] = evTraj;
    mapParameterValues["coordinate_file"] = evGro;
    mapParameterValues["grid_file"] = evGrf;
    mapParameterValues["log_file"] = evLog;
    mapParameterValues["mu"] = evMu;
    mapParameterValues["hinit"] = evHinit;
    mapParameterValues["trajectory_writing_frequency"] = evNxtc;
    mapParameterValues["l"] = evl;
    mapParameterValues["grid"] = evGrid;
    mapParameterValues["box"] = evBox;
    mapParameterValues["field"] = evField;
    mapParameterValues["diffusion"] = evDiffusion;
}


inputrec :: inputrec ( )
{
    nsteps = 2000;
    gamma  = 72.0;                    // mN/m 
    T      = 300;                     // K
    mu     = 100;                     // amu / nm**2
    l      = 0.9;                     // nm
    dt     = 0.05;                    // ps
    tau    = 100*dt;                  // ps
    hinit  = 2.5;                     // nm (half of default box[ZZ])
    nxtc   = 1;
    gro    = "system.gro";
    traj   = "traj.xtc";
    log    = "md.log";
    grf    = "interface.grid";
    fieldfile   = "";
    bhinit = false;
    bfield  = false; 
    bgrid  = false;
    bdiffusion = true;
}

void inputrec :: get_internal_units ( )
{
    /* 
       Converts the input to consitent set of internal units.
       Refer main.h for the input units.
    */
    gamma = gamma * MJPM2TOAMUPPS2; // mN/m to amu/ps**2
}


// String routines 
double string_to_double( const std::string& s )
{
    std::istringstream i(s);
    double x;
    if (!(i >> x))
	return 0;
    return x;
} 

int string_to_int( const std::string& s )
{
    std::istringstream i(s);
    int x;
    if (!(i >> x))
	return 0;
    return x;
} 

long string_to_long( const std::string& s )
{
    std::istringstream i(s);
    long x;
    if (!(i >> x))
	return 0;
    return x;
} 

char *fgets2(char *line, int n, FILE *stream)
{
    /* 
       This routine reads a string from stream of max length n
       and zero terminated, without newlines line should be long 
       enough (>= n)
    */
    char *c;
    if (fgets(line,n,stream)==NULL) return NULL;
    if ((c=strchr(line,'\n'))!=NULL) *c=0;
    if ((c=strchr(line,'\r'))!=NULL) *c=0;
    return line;
}

void ltrim (char *str)
{
    /* 
       Removes all leading spaces from a string 
    */
    char *tr;
    int c;
    
    if (!str)
	return;
    
    tr = strdup (str);
    c  = 0;
    while ((tr[c] == ' ') || (tr[c] == '\t'))
	c++;
    
    strcpy (str,tr+c);
    free (tr);
}

void rtrim (char *str)
{
    /* 
       Remove all trailing spaces from a string 
    */
    int nul;
    
    if (!str)
	return;
    
    nul = strlen(str)-1;
    while ((nul > 0) && ((str[nul] == ' ') || (str[nul] == '\t')) ) {
	str[nul] = '\0';
	nul--;
    }
}

void trim (char *str)
{
    /* 
       Remove leading and trailing spaces from a string 
    */
    ltrim (str);
    rtrim (str);
}

std::string trim_right_copy ( const std::string& s,
				     const std::string& delimiters )
{
    /* 
       Remove trailing spaces and tabs from a string 
    */
    if ( s.length() > 0 )
    {
	unsigned pos = s.find_last_not_of( delimiters ) + 1;
	if ( pos > 0 )
	{ 
	    return s.substr( 0, pos );
	}
	else 
	{
	    return "";
	}
    }
    else
    {
	return "";
    }
}

std::string trim_left_copy ( const std::string& s, 
				    const std::string& delimiters )
{ 
    /* 
       Remove leading spaces and tabs from a string 
    */
    if ( s.length() > 0 )
    {
	unsigned pos = s.find_first_not_of( delimiters );
	if ( pos < s.length ( ) )
	{
	    return s.substr( s.find_first_not_of( delimiters ) );
	}
	else
	{
	    return "";
	}
    }
    else
    {
	return "";
    }
}

std::string trim_copy( const std::string& s,
			      const std::string& delimiters )
{
    /* 
       Remove leading and trailing spaces and tabs from a string 
    */
    return trim_left_copy( trim_right_copy( s, delimiters ), delimiters );
}


// Parsing routines
void parse_parameter_file ( const char *fn, inputrec* inp )
{
    /* 
       Parse the input file 
    */
    std::ifstream   in;
    std::string     line;
    int             lc;

    // Open the file for reading
    in.open(fn);
    if (in.is_open()) 
    { 
	std::cout << "Opened " << fn << " for reading\n"; 
    }
    else
    {
	std::cout << "Cannot open " << fn << " for reading\n"; 
	exit( EXIT_FAILURE );
    }
    
    // Initialize the mapping
    initialize_map();
    std::cout << "Initialized map\n";

    // Parse the file
    lc = 0;
    while( std::getline(in, line) )
    {
	lc++;
	line = line.substr(0,line.find_first_of("#;"));
	line = trim_copy(line);
	
	if ( line.length() > 0 )
	{
	    // line store in args structure
	    std::string             lbuf,rbuf,buf;
	   	    
	    lbuf = trim_copy(line.substr(0,line.find_first_of("=")));
	    if ( lbuf.length() <= 0 ) 
	    {
		std::cout << "Error on line " << lc << " in paramter file\n";
		exit( EXIT_FAILURE );	
	    } 
	    else
	    {
		rbuf = line.substr(line.find_first_of("=")+1);
		if (rbuf.length() <= 0 )
		{
		    std::cout << "Error on line " << lc << " in paramter file\n";
		    exit( EXIT_FAILURE );  
		} 
		else
		{
		    std::vector<std::string>  rbufData;
		    std::string               value;
		    std::stringstream         rbufStream(trim_copy(rbuf));
		    int nval = 0;
		    while(rbufStream >> value)
		    {
			rbufData.push_back(value);
			nval++;
		    }
		    if (nval == 1)
		    {
			switch ( mapParameterValues[ strdup ( lbuf.c_str() ) ] )
			{
			case evT:
			    inp->T = string_to_double(rbufData[0]);
			    break;
			case evdt:
			    inp->dt = string_to_double(rbufData[0]);
			    break;
			case evNsteps:
			    inp->nsteps = string_to_long(rbufData[0]); 
			    break;
			case evGamma:
			    inp->gamma = string_to_double(rbufData[0]); 
			    break;
			case evTau:
			    inp->tau = string_to_double(rbufData[0]); 
			    break;
			case evMu:
			    inp->mu = string_to_double(rbufData[0]); 
			    break;
			case evHinit:
			    inp->hinit = string_to_double(rbufData[0]); 
			    inp->bhinit = true;
			    break;
			case evTraj:
			    inp->traj = rbufData[0];
			    break;
			case evGro:
			    inp->gro = rbufData[0];
			    break;
			case evGrf:
			    inp->grf = rbufData[0];
			    break;
			case evLog:
			    inp->log = rbufData[0];
			    break;
			case evl:
			    inp->l = string_to_double(rbufData[0]);
			    break;
			case evNxtc:
			    inp->nxtc = string_to_int(rbufData[0]);
			    break;
			case evField:
			    inp->fieldfile = rbufData[0];
			    inp->bfield = true;
			    std::cout << "Field file = " << inp->fieldfile <<"\n";
			    break;
			case evDiffusion:
			    std::transform(rbufData[0].begin(), rbufData[0].end(), rbufData[0].begin(), ::tolower);
			    if ( ! rbufData[0].compare("yes")) 
			    {
				inp->bdiffusion = true;
			    }
			    else if ( ! rbufData[0].compare("no"))
			    {
				inp->bdiffusion = false;
			    } 
			    else
			    {
				std::cout << "Error on line "<< lc << " in paramter file.\n"; 
				std::cout << line <<"\n";
				exit( EXIT_FAILURE );
			    }
			    break;
			default:
			    std::cout << "Error on line "<< lc << " in paramter file.\n"; 
			    std::cout << line <<"\n";
			    exit( EXIT_FAILURE );
			}
		    }
		    else if ( nval == 2 )
		    {
			switch ( mapParameterValues[ strdup ( lbuf.c_str() ) ] )
			{
			case evGrid:
			    inp->bgrid = true;
			    inp->grid.dx = string_to_double(rbufData[0]);
			    inp->grid.dy = string_to_double(rbufData[1]);
			    break;
			default: 
			    std::cout << "Error on line "<< lc << " in paramter file.\n";
			    std::cout << line <<"\n";
			    exit( EXIT_FAILURE );
			}
		    }
		    else if (nval == 3)
		    {
			switch ( mapParameterValues[ strdup ( lbuf.c_str() ) ] )
			{
			case evBox:
			    inp->bgrid = true;
			    inp->grid.box[0] = string_to_double(rbufData[0]);
			    inp->grid.box[1] = string_to_double(rbufData[1]);
			    inp->grid.box[2] = string_to_double(rbufData[2]);
			    break;
			default: 
			    std::cout << "Error on line "<< lc << " in paramter file.\n";
			    std::cout << line <<"\n";
			    exit( EXIT_FAILURE );
			}
		    }
		    else
		    {
			std::cout << "Error on line " << lc << " in paramter file.\n";
			std::cout << line <<"\n";
			exit( EXIT_FAILURE );
		    }
		    rbufData.clear();
		}
	    }
	}
    }

    // Perform some checks on the inputs
    // Set the grid variable inside the inputrec
    if ( inp->bgrid )
    {
	inp->grid.grid_reinit ( ); 
    }
    if ( inp->bhinit )
    {
	if ( inp->hinit < 0.0 || inp->hinit > inp->grid.box[ZZ] )
	{
	    std::cout << "Initial height needs to between the box limits of "
		      << 0.0 << " and " << inp->grid.box[ZZ] << "\n";
	    exit ( EXIT_FAILURE );
	}
    }
    inp->get_internal_units ( ); 
}


// Output writing routines
void print_usage ( void )
{
    /* 
       Routine to print help information about the program
    */
    const char *line[] = { " Elastic membrane model simulation\n",
			   " Usage: ./em_interface -f inputs.parm \n",
			   " Please refer the README.md document for the\n",
			   " parameter file (.parm) format.\n" };
    
    for ( int i = 0; i < 4; i++ )
    {
	std::cout << line[i] ;
    }
}

void print_input_parameters ( inputrec* inp )
{
    /*
      Print the input parameters on screen
    */
    std::cout << "Temperature = " << inp->T      << "\n";
    std::cout << "Timestep    = " << inp->dt     << "\n";
    std::cout << "Tau         = " << inp->tau    << "\n";
    std::cout << "Log         = " << inp->log    << "\n";
    std::cout << "Gro         = " << inp->gro    << "\n";
    std::cout << "Xtc         = " << inp->traj   << "\n";
    std::cout << "nstep       = " << inp->nsteps << "\n";
    std::cout << "gamma       = " << inp->gamma  << "\n";
    std::cout << "mu          = " << inp->mu     << "\n";
    std::cout << "l           = " << inp->l      << "\n";
    std::cout << "GridL       = " << inp->grid.box[XX] << " x " 
	      << inp->grid.box[YY] << "\n";
    std::cout << "GridD       = " << inp->grid.nx      << " x " 
	      << inp->grid.ny      << "\n";
}

void printlog_input_parameters ( inputrec* inp, std::ofstream& fp )
{
    /*
      Print input parameters to the log file
    */
    fp  << "Temperature = " << inp->T      << "\n";
    fp  << "Timestep    = " << inp->dt     << "\n";
    fp  << "Tau         = " << inp->tau    << "\n";
    fp  << "Log         = " << inp->log    << "\n";
    fp  << "Gro         = " << inp->gro    << "\n";
    fp  << "Xtc         = " << inp->traj   << "\n";
    fp  << "nstep       = " << inp->nsteps << "\n";
    fp  << "gamma       = " << inp->gamma  << "\n";
    fp  << "mu          = " << inp->mu     << "\n";
    fp  << "l           = " << inp->l      << "\n";
    fp  << "GridL       = " << inp->grid.box[XX] << " x " 
	<< inp->grid.box[YY] << "\n";
    fp  << "GridD       = " << inp->grid.nx      << " x " 
	<< inp->grid.ny      << "\n";
}

void write_gro ( const char *gro, double* h, t_Grid* grid ) 
{
    /* 
       Write a gro file with coordinates
    */
    FILE* cnugro;    
    int   kx, ky, i;
    float x, y;
    
    cnugro = fopen ( gro, "w" );
    fprintf (cnugro, "Test output gro\n");
    fprintf (cnugro, "%d\n",grid->nx*grid->ny);
    for ( kx = 0; kx < grid->nx; kx++ )
    {
	for ( ky = 0; ky < grid->ny; ky++ )
	{
	    i = kx*grid->ny + ky;
	    x = kx*grid->dx;
	    y = ky*grid->dy;
	    fprintf ( cnugro, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", 
		      1,"PT","PT",i,x,y,h[i] ); 
	}
    }
    fprintf ( cnugro, "%f %f %f\n", 
	      grid->box[XX], grid->box[YY], grid->box[ZZ] );
    fclose ( cnugro );
}

int writeframe(double *h, t_Grid* grid, XDRFILE* xd, int step, 
	       float xtctime, matrix box) 
{
    /*
      Write a single frame of data to the XTC file pointed to by
      the variable xd.
    */
    int     npt,kx,ky;
    rvec    *coord;
    float   prec;
    int     ret;
    
    // Set precision for XTC
    prec = 1000000;
    
    // Create coordinate array
    coord = new rvec[grid->nx*grid->ny]; 
    npt = 0;
    for( kx = 0; kx < grid->nx; kx++ ) 
    {
	for( ky = 0; ky < grid->ny; ky++ ) 
	{
	    coord[npt][0] = kx*grid->dx;
	    coord[npt][1] = ky*grid->dy;
	    coord[npt][2] = h[kx*grid->ny+ky];
	    npt = npt + 1;
	}
    }
  
    // Write XTC frame
    ret = write_xtc(xd,npt,step,xtctime,box,coord,prec);
    
    delete [] coord;
    return ret;
} 

