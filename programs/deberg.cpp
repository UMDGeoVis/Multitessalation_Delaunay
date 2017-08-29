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

// ----------------------------------------------------------------------------------------
//
//  file   : deberg.cpp
//  author : Christian Melchiorre
//
//  Implementation of class DeBerg, which exports function SelectVertices.
//  Such function takes a set (structured as a tree) of removable vertices
//  in a triangulation, and returns a list of a subset of independent
//  vertices in it.
//  Such function will be used by classes TDecRndDeBerg (sub-class of
//  TDecRndDelaunay) and TDecErrDeBerg (sub-class of TDecErrDelaunay).
//  rispettivamente, da TDecRndDelaunay e TDecErrDelaunay.
//  The solution of a static class avoids using multiple inheritance
//  with all problems connected to virtual base classes etc...
//

#include "defs.h"
#include "geom.h"
#include "error.h"
#include "ttriang.h"
#include "tbtree.h"
#include "tlist.h"
#include "tdoublelist.h"
#include "deberg.h"


// ----------------------------------------------------------------------------------------
//
//  static void DeBerg::SelectVertices( TBTree<PTVertex>& ElimVtxTree,
//                                      TDoubleList<PTVertex>& DeBergVertices )
//
//  This is the only function exported by class DeBerg, and performs the
//  task described in the header of this file. Vertices of tree 
//  ElimVtxTree must be removable vertices of a triangulation which they
//  belong to (i.e., they are consistently connected with other entities of
//  the triangulation through topological relations). Starting from one
//  vertex chosen as the minimum element of the tree (minimality will have
//  a different meaning when SelectVertices is called from 
//  TDecRndDeBerg or TDecErrDeBerg), extract a set of such vertices which
//  are independent from each other (i.e., for each pair of vertices 
//  Vi, Vj, there is no edge (Vi, Vj) in the triangulation). Extracted
//  vertices are returned in list DeBergVertices, passed as a parameter by
//  reference.
//

void DeBerg::SelectVertices( TBTree<PTVertex>& ElimVtxTree, TDoubleList<PTVertex>& DeBergVertices )
{
   
    int vid;

    #ifdef DEBUG
       DEBUG << "DeBerg::SelectVertices()" << endl;
    #endif

    //
    // Transform the tree in a sorted list.
    // At the same time, take note of the maximum VIS among the vertices
    // in the tree, useful to dimension the array of marks 
    // (which is indexed on VIDs).
    //
    
    TList<PTVertex> ElimVtxList;
    
    int MaxVID = 0;
    
    while ( !ElimVtxTree.IsEmpty() )
    {
       PTVertex NextV = ElimVtxTree.RemoveMax();
       if ( NextV->VID > MaxVID ) MaxVID = NextV->VID;
       ElimVtxList.AddHead( NextV );
    }
       
    //
    // Creare array of marks. Such array is indexed on vertex VIDs and
    // contains a boolean mark: Marked[vid] is TRUE if some adjacent vertex
    // of vid has been selected as part of list DeBergVertices.
    //   
        
    boolean *Marked = new boolean[MaxVID+1];
    
    check( (Marked==NULL), "DeBerg::SelectVertices(), insufficient memory");
    
    for( vid=0; vid<=MaxVID; vid++ )
       Marked[vid] = FALSE;
    
    //
    // Main loop of function SelectVertices: extract first vertex from
    // list ElimVtxList and, if not marked, mark all its neighbors, and
    // insert the vertex in the ouput list DeBergVertices. If it is 
    // marked, inert it back into the tree.
    //

    while( !ElimVtxList.IsEmpty() )
    {
		PTVertex NextV = ElimVtxList.RemoveHead();

		if ( Marked[ NextV->VID ] )
		{
			ElimVtxTree.Insert( NextV );
		}
		else
		{
			MarkAllNeighbours( NextV, Marked, MaxVID );
			DeBergVertices.AddTail( NextV );
		}
	} // end ...while( !ElimVtxList.IsEmpty() )

    
    //
    // End of traversal. Free the array of marks.
    //
    
    delete[] Marked;
    
}


// ----------------------------------------------------------------------------------------
//
//  static void DeBerg::MarkAllNeighbours( PTVertex V, boolean *Marked )
//
//  Given vertex V, mark (by setting a TRUE value in the corresponging
//  positions of array Marked[] passed as a parameter) all vertices
//  adjacent to V (i.e., endpoints of the same edge).
//

void DeBerg::MarkAllNeighbours( PTVertex V, boolean *Marked, int MaxVID )
{

   int i, j;

   #ifdef DEBUG
     DEBUG << "DeBerg::MarkAllNeighbours( V" << V->VID << " )" << endl;
   #endif

   //
   // Usual traversal of the vertices of the polygon of influence of V,
   // as done, for instance, in TDecErrDelaunay::CalcVertexErrorXXX().
   //

   PTTriangle TFirst;
   
   PTEdge   EFirst = V->VE[0];
   PTVertex VFirst = ( EFirst->EV[0] != V ? EFirst->EV[0] : EFirst->EV[1] );

   if ( EFirst->OnConvexHull() )
   {
       
       TFirst = ( EFirst->ET[0] != NULL ? EFirst->ET[0] : EFirst->ET[1] );
  
       #ifdef ROBUST
          check( (TFirst == NULL), "TDecErrDelaunay::RecalcVertexErrorApprox(), <0> inconsistency detected");
       #endif
            
       PTVertex v[3];
       TFirst->GetTV( v[0], v[1], v[2] );
       
       for( i=0; i<3; i++ )
          if ( v[i] != V && v[i] != VFirst ) break;
	 
       #ifdef ROBUST
          check( (i>=3), "TDecErrDelaunay::RecalcVertexErrorApprox(), <1> inconsistency detected" );
       #endif
	  
       if ( Geom::Turnxy( V, VFirst, v[i] ) != TURN_LEFT )
          EFirst = V->VE[1];
      
   } // end ...if( EFirst->OnConvexHull() )
   
   
   // triangle after EFirst, in counterclockwise order
   
   TFirst = ( EFirst->EV[0] == V ? EFirst->ET[0] : EFirst->ET[1] );
   
   #ifdef ROBUST
      check( ( TFirst==NULL), "TDecErrDelaunay::RecalcVertexErrorApprox(), <2> inconsistency detected" );
   #endif
   
   PTEdge ENext = EFirst;
   PTTriangle TNext;

   do   
   {
      // the endpoint of ENext different from V belongs to the boundary
      // of the influence region of V
      
      PTVertex VNext = ( ENext->EV[0] != V ? ENext->EV[0] : ENext->EV[1] );

      // mark the vertex
            
      if (VNext->VID <= MaxVID)
		Marked[ VNext->VID ] = TRUE;
      
      // go to the next one
      
      TNext = ( ENext->EV[0] == V ? ENext->ET[0] : ENext->ET[1] );
      
      if ( TNext != NULL )
      {
         // search for index of E in the TE relation of TNext
	 
	 for( j=0; j<3; j++ )
	   if ( TNext->TE[j] == ENext ) break;
	   
	 #ifdef ROBUST
	    check( (j>=3), "TDecErrDelaunay::RecalcVertexErrorApprox(), <3> inconsistency detected" );
	 #endif

         ENext = TNext->TE[(j+2)%3]; // edge preceding E in TE of TNext
	 
      }
           
   } while ( TNext != NULL && ENext != EFirst );
   

   // end of traversal
     
}
