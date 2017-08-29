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
//  file   : decrnddbCDT.cpp
//  author : Alessio Calcagno
//
//  Implementation of class TDecRndDeBergCDT, sub-class of TDecRndCDT
//  and TDecRndDeBerg, for
//  the decimation of a constrained Delaunay triangulation (CDT) with
//  removal of independent vertices and a random choice of vertices to be
//  removed.
//  To this aim, procedure DeBerg::SelectVertices of class DeBerg is used
//  (files: deberg.h/deberg.cpp).
//  The only re-implemented functions are NoMoreUpdates (implemented
//  inline in file decrnddbCDT.h) and NextVertex().
//



#include "defs.h"
#include "error.h"
#include "ttriang.h"

#include "tbtree.h"
#include "tdoublelist.h"

#include "deberg.h"
#include "decrnddbCDT.h"



// ---------------------------------------------------------------------------------
//
//  void TDecRndDeBergCDT::NextVertex()
//
//  Choose the next vertex in the following way: if there is still some
//  vertex in list DeBergVertices, take the head of such list;
//  otherwise, by calling DeBerg::SelectVertices(), create a new list
//  of independent removable vertices and take the head of such new list.
//


/* PAOLA: In .h it is commented
void TDecRndDeBergCDT::NextVertex()
{

    #ifdef DEBUG
       DEBUG << "\nTDecRndDeBergCDT::NextVertex()" << endl;
    #endif

    #ifdef ROBUST
       check( (DeBergVertices.IsEmpty() && ElimVtxTree.IsEmpty()), 
              "TDecRndDeBergCDT::NextVertex(), <1> no more points to remove" );
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

*/
