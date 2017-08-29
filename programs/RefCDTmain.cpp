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

/* ------------------------------------------------------------------------ */
/*           CONSTRUCTION OF A CONSTRAINED DELAUNAY TRIANGULATION           */
/* ------------------------------------------------------------------------ */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "error.h"
#include "geom.h"

#include "mttracer.h"

#include "builddel.h"
#include "refdel.h"
#include "decCDT.h"
#include "refCDT.h"

int main( int argc, char **argv )
{
   PTTriangulation T;
   MTTracer MT;

   cerr.precision(16);
   cerr << endl;

   if ( argc != 3 )
   {
      cerr << "usage: " << argv[0] << " infile outfile" << endl;
      exit(-1);
   }
 
   //
   // Create triangulation
   //
       
   T = new TRefCDT( &MT );
   check( (T == NULL), "INSUFFICIENT MEMORY" );
   cerr << endl;

   //
   // Execute triangolation
   //

   T->BuildTriangulation( argv[1], argv[2] );
   cerr << "triangulation completed" << endl << endl;
   return 0;
}
