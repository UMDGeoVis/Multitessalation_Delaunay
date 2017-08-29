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
//   file   : error.h
//   author : Christian Melchiorre
//
//   Routines for the management of error messages
//

#ifndef _ERROR_H
#define _ERROR_H

#include "defs.h"

// Show an error message, stop execution and ask for confirmation to
// resume exectution.
void Pause( char *msg );

// Show an error message, terminate execution.
void error( char * );

// If condition is true, show an error message and terminate execution.
void check( boolean, char * ); 
  
// Show an error message.
void message( char * );

#endif // _ERROR_H
