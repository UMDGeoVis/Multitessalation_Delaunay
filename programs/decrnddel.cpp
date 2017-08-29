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

// -----------------------------------------------------------------------------
//  
//   file   : decrnddel.h
//   author : Christian Melchiorre
//
//   Implementation of class TDecRndDelaunay, sub-class of TDecimDelaunay,
//   for decimation of a Delaunay triangulation with a random choice of
//   the vertex to be removed.
//   This class re-implements method InitialTriangulation in this way:
//   before inserting the removable vertices
//   in ElimVtxTree, it sets a random key in their Error field
//   (function compare for PTVertex used in TBTree uses this field).
//

#include <stdlib.h>
#include <time.h>

#include "defs.h"
#include "error.h"
#include "ttriang.h"
#include "decrnddel.h"


#ifdef CC_VISUAL5
#define srand48 srand
#define drand48() ((double) rand() / RAND_MAX)
#endif


// -----------------------------------------------------------------------------
//  
//  void TDecRndDelaunay::InitialTriangulation()
//
//  Similar to TDecimDelaunay::InitialTriangulation. The difference is
//  that, before inserting the removable vertices
//  in ElimVtxTree, it sets a random key in their Error field 
//  (function compare for PTVertex used in TBTree uses this field).
//


void TDecRndDelaunay::InitialTriangulation()
{
    #ifdef DEBUG
       DEBUG << "TDecRndDelaunay::InitialTriangulation()" << endl;
    #endif
       
    srand48( (int) time(NULL) );

    for( int v=0; v<nPts; v++ )
	 Points[v]->Error = drand48();
	
       
    TDecimDelaunay::InitialTriangulation();
    
}
