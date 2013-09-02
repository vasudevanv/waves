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

  sim_harmonicpin.cpp

  Routines to deal with 9-3 surface. 

*/

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <cstring>
#include <cmath>
#include <algorithm> 
#include "constants.h"
#include "grid.h"
#include "sim_io.h"
#include "sim_field.h"
#include "sim_surface93.h"
#include "sim_harmonicpin.h"
#include "sim_samsurface.h"

using std::vector;

Field:: Field()
{

}
Field:: ~Field()
{

}
void Field:: init ( int argc, char **argv, std::ofstream& fp ) 
{
}

void Field:: calc_forces( double*h, double*f, t_Grid* grid )
{
}

Field* createField ( const char *type )
{
    if( strcmp( type, "sam" ) == 0 )
    {
       return new samsurface();
    }
    else if( strcmp( type, "surface93" ) == 0 )
    {
	return new surface93();
    }
    else if( strcmp(type, "harmonicpin") == 0 )
    {
	return new harmonicpin();
    }
    else
    {
	return NULL;
    }
}


vector<FieldArgs>* parse_field_file ( const char *fn )
{
   
    std::ifstream              in;
    std::string                line;
    int                        lc;
    vector<FieldArgs> *result = new vector<FieldArgs>();
    
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
    
    lc = 0;
    while( std::getline(in, line) )
    {
	line = line.substr(0,line.find_first_of("#;"));

	if ( line.length() > 0 && trim_copy(line).length() > 0 )
	{
	    
	    // line store in args structure
	    vector<std::string>     lineData;
	    std::string             value;
	    std::stringstream       lineStream(trim_copy(line));
	    
	    int nval = 0;
	    while(lineStream >> value)
	    {
		lineData.push_back(value);
		nval++;
	    }
	    if ( nval < 1 )
	    {
		std::cout << "Error on line " << lc << " in field file\n";
		exit ( EXIT_FAILURE );
	    } 
	    else
	    {
		FieldArgs args;
		args.argc = (int)lineData.size();
		args.argv = new char*[ args.argc ];
		args.argv = new char*[ line.length() + 1 ];
		for( int i = 0; i < args.argc; i++ )
		{
		    args.argv[i] = strdup(lineData[i].c_str());
		}
		args.line = strdup(line.c_str());
		std::cout << "\n";
		lineData.clear();
		result->push_back( args );
		lc++;
	    }
	}
    }

    if ( lc == 0 )
    {
	std::cout << "Error in field file: No fields found \n";
	exit ( EXIT_FAILURE );
    }
    return result;
}

void dispose_parsed_fields( vector<FieldArgs> *parsed )
{
    vector<FieldArgs>::iterator iter, begin, end;
    begin = parsed->begin();
    end = parsed->end();
    for( iter = begin; iter != end; ++iter )
    {
	FieldArgs& args = *iter;
	for( int i = 0; i < args.argc; i++ )
	    delete [] args.argv[i];
	delete [] args.argv;
    }
    delete parsed;
} 

int check_surface ( const char* typeName )
{
    std::vector<std:: string> mylist;
    mylist.push_back("surface93");
    mylist.push_back("samsurface");
    std::string myinput(typeName);
    if (std::find(mylist.begin(), mylist.end(), myinput) != 
	mylist.end())
	return 1;
    else
	return 0;
} 

void init_fields ( const char* fn, t_Field* fd, double T, std::ofstream& fp )
{
    // Check flag to make sure only one surface is given
    int    surface_checkflag = 0;

    // Routine to parse the umbrella 
    vector<FieldArgs> *parsed = parse_field_file( fn );
    
    vector<Field*> *fieldList = new vector<Field*>();
    fd->sFieldList = fieldList;
    
    // Go through each external force in turn, decode its
    // type, request the main group, and give a chance
    // to the umbrella to request other groups
    for( vector<FieldArgs>::iterator iter = parsed->begin();
	 iter != parsed->end(); ++iter )
    {
	FieldArgs& args = *iter;
	
	if( args.argc < 3 )
	{
	    std:: cout << "Invalid field description: "<< args.line << "\n";
	    exit ( EXIT_FAILURE );
	}
	else
	{
	    const char *typeName = args.argv[0];
	    
	    std::cout << "Found field with type name " << typeName << "\n";
	    Field *U = createField( typeName );
	    if( U == NULL )
	    {
		std::cout << "Invalid field type " << typeName << "\n";
		exit ( EXIT_FAILURE );
	    }
	    else
	    {
		if ( check_surface(typeName) && surface_checkflag == 1 )
		{
		    std::cout << "Can have only one surface in simulation\n"; 
		    exit ( EXIT_FAILURE );
		}
		else if ( check_surface(typeName) && surface_checkflag == 0 )
		{
		    surface_checkflag = 1;
		}
		U->T = T;
		U->init ( args.argc, args.argv, fp );
		fieldList->push_back( U );
	    }
	}
    }
    dispose_parsed_fields( parsed );
}

void finish_fields( t_Field* fd )
{
    vector<Field*> *fieldList = static_cast< vector<Field*>* >( fd->sFieldList );
    for( vector<Field*>::iterator iter = fieldList->begin();
	 iter != fieldList->end(); ++iter )
    {
	Field *U = *iter;
	delete U;
	*iter = NULL;
    }
    delete fieldList;
    fd->sFieldList = NULL;
} 


void field_potential ( double* h, double* f, t_Grid* grid, t_Field* fd )
{
    /*
     Function to calculate the force due to all fields 
    */
    vector<Field*> *fieldList = static_cast< vector<Field*>* >( fd->sFieldList );
    int kx,ky,el;
    for ( kx = 0; kx < grid->nx; kx++ )
    {
        for ( ky = 0; ky < grid->ny; ky++ )
        {
            el = kx* grid->ny + ky;
            f[el] = 0.0;
        }
    }
    for( vector<Field*>::iterator iter = fieldList->begin();
	 iter != fieldList->end(); ++iter )
    {
        Field *U = *iter;
	U->calc_forces( h, f, grid );
    }
}


