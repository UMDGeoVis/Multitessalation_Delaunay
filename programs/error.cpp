/*****************************************************************************
Delaunay Triangulator and MT constructor, version 1.0, 1999.
Copyright (C) 1999 DISI - University of Genova, Italy.
Group of Geometric Modeling and Computer Graphics DISI.
DISI - University of Genova, Via Dodecaneso 35, 16146 Genova - ITALY.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*****************************************************************************/

// ---------------------------------------------------------------------------
//
//   file   : error.cpp
//   author : Christian Melchiorre
// 
//   Routines for the management of error messages
//

#include <iostream>
#include <stdlib.h>

#include "defs.h"
#include "error.h"

using namespace std;

//
// Show an error message, stop execution and ask for confirmation to
// resume exectution.
//

void Pause( char *msg )
{
   cerr << "\n" << msg << "\n";
   char c;
   do 
   {
        cerr << "Proseguire ? (S/N) " << endl;
	     cin >> c;
	}  while( c != 's' && c != 'S' && c != 'n' && c != 'N' );

   if ( c=='n' || c=='N' ) exit(-1);
}


//
// Show an error message, terminate execution.
//

void error( char *msg )
{
   cerr << "\nERROR: " << msg << "\n";
   exit(-1);
}


//
// If condition is true, show an error message and terminate execution.
//

void check( boolean cond, char *msg )
{
   if (cond) error(msg);
}


//
// Show an error message on cerr.
//

void message( char *msg )
{
   cerr << "\nMESSAGE: " << msg << "\n";
}
