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

// ---------------------------------------------------------------------------------
//
//  file   : decrnddb.cpp
//  author : Christian Melchiorre
//
//  Implementation of class TDecRndDeBerg, sub-class of TDecRndDelaunay, for
//  decimation of a Delaunay triangulation with removal of independent
//  vertices and a random choice of vertices to be removed.
//  To this aim, procedure DeBerg::SelectVertices of class DeBerg is used
//  (files: deberg.h/deberg.cpp).
//  The only re-implemented functions are NoMoreUpdates (implemented
//  inline in file decrnddeb.h) and NextVertex().
//



#include "defs.h"
#include "error.h"
#include "ttriang.h"

#include "tbtree.h"
#include "tdoublelist.h"

#include "deberg.h"
#include "decrnddb.h"


// ---------------------------------------------------------------------------------
//
//  void TDecRndDeBerg::NextVertex()
//
//  Choose the next vertex in the following way: if there is still some
//  vertex in list DeBergVertices, take the head of such list; 
//  otherwise, by calling DeBerg::SelectVertices(), create a new list
//  of independent removable vertices and take the head of such new list.
//

void TDecRndDeBerg::NextVertex()
{

    #ifdef DEBUG
       DEBUG << "\nTDecRndDeBerg::NextVertex()" << endl;
    #endif

    #ifdef ROBUST
       check( (DeBergVertices.IsEmpty() && ElimVtxTree.IsEmpty()), 
              "TDecRndDeBerg::NextVertex(), <1> no more points to remove" );
    #endif
 
    if ( DeBergVertices.IsEmpty() )
    {
 
        DeBerg::SelectVertices( ElimVtxTree, DeBergVertices );
 
	#ifdef DEBUG
        
	   // show the content of DeBergVertices
	
	   DEBUG << "DeBergVtx : ";
	   TDoubleListIterator<PTVertex> I( &DeBergVertices );
	   while( !I.EndOfList() )
	   {
	      DEBUG << " -> V" << (I.Current()->object)->VID;
	      I.GoNext();
	   }
	   DEBUG  << endl << endl;
	
	#endif // DEBUG
	
        #ifdef ROBUST
           check( (DeBergVertices.IsEmpty()),
                "TDestroyDelaunay::NextVertex(), <2> no more points to remove" );
        #endif
 
    } // end ...if ( DeBergVertices.IsEmpty() )
 
    
    VertexToRemove = DeBergVertices.RemoveHead();

    #ifdef DEBUG    
       DEBUG << "VertexToRemove: V" << VertexToRemove->VID  << endl;
    #endif

}
