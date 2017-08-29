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

// -------------------------------------------------------------------------
//
//   file   : referrdel.cpp
//   author : Christian Melchiorre
//
//   Implementation of class TRefErrDelaunay, sub-class of TRefineDelaunay,
//   for building a Delaunay triangulation through refinement, with an
//   error-based choice of the point to be inserted.
//   With respect to the base class, it redefines method NextPoint():
//   the next point to be inserted is the one that causes the maximum
//   error w.r.t. the triangle in which it falls.
//   Moreover, it redefines method UpdateStep(): new triangles created
//   in the current update step are inserted in a balanced binary tree
//   of triangles (in order to find in logarithmic time the triangle
//   containing the point of maximum error).
//


#include <iostream>

#include "defs.h"
#include "error.h"
#include "geom.h"

#include "tbtree.h"

#include "builddel.h"
#include "refdel.h"
#include "referrdel.h"

#include "mttracer.h"


// -------------------------------------------------------------------------
//
//  Constructor of class TRefErrDelaunay.
//


TRefErrDelaunay::TRefErrDelaunay( PMTTracer iMT ) : TRefineDelaunay( iMT ), PtsErrTree()
{
   #ifdef DEBUG
       DEBUG << "TRefErrDelaunay() Constructor" << endl;
   #endif
}


// -------------------------------------------------------------------------
//
//  void TRefErrDelaunay::InitialTriangulation()
//
//  Redefine TRefineDelaunay::InitialTriangulation() in such a way to
//  call AddPointsToTree() after building the initial triangulation.
//

void TRefErrDelaunay::InitialTriangulation()
{
   #ifdef DEBUG
      DEBUG << "TRefErrDelaunay::InitialTriangulation()" << endl;
   #endif
   
   TRefineDelaunay::InitialTriangulation();
   
   AddAllPointsToTree();

}

      
// -------------------------------------------------------------------------
//
//  void TRefErrDelaunay::AddAllPointsToTree()
//
//  Add to PtsErrTree, the balanced binary tree of PTPoint, the remaining
//  points after the initial triangulation of the convex hull.
//
    
void TRefErrDelaunay::AddAllPointsToTree()
{    
    #ifdef DEBUG
       DEBUG << "TRefErrDelaunay::AddAllPointsToTree()" << endl;
    #endif
    
    for( int i=nChPts; i<nPts; i++ )
       PtsErrTree.Insert( Points[i] );
    
}


// -------------------------------------------------------------------------
//
//  void TRefErrDelaunay::NextPoint()
//
//  In the base classes, the new vertex was created starting from the next
//  point in array Points (the point of index iNextPoint). 
//  Here, the next vertex is found as the one causing the maximum error
//  in the triangles of the current triangulation.
//  To this aim, tree PtsErrTree is used. In such tree, points not yet
//  belonging to the triangulation are sorted based on their error.
//

void TRefErrDelaunay::NextPoint()
{

   #ifdef DEBUG
       DEBUG << "TRefErrDelaunay::NextPoint()" << endl;
   #endif

   #ifdef ROBUST
       check( ( iNextPoint < nChPts ),
            "TRefErrDelaunay::NextPoint(), this should not happen" );
       check( (PtsErrTree.IsEmpty()),
            "TRefErrDelaunay::NextPoint(), point tree empty" );
   #endif
   
   PTPoint PointToIns = PtsErrTree.RemoveMax();

   #ifdef DEBUG
      DEBUG << "Max point from the TREE : " << (*PointToIns) << endl;
   #endif

   // PointToIns is the point causing the max error
   
   Points[ iNextPoint ] = VertexToIns = new TVertex( PointToIns );
   
   check( (VertexToIns == NULL), "TRefErrDelaunay::NextPoint(), insufficient memory");
   
   iNextPoint++;
  
   //
   // In the base classes, variable iNextPoint was the index in array
   // points of the next point to be inserted, and also the array entry
   // in which to put VertexToIns as soon as it is created.
   // In TRefErrDelaunay, variable iNextPoint is used only for the latter
   // purpose.
   // Maintaining new vertices sorted in order of creation in array Points
   // is useful to write them in order of insertion in the triangulation.
   //
         
}




// -----------------------------------------------------------------------------
//  
//  void TRefErrDelaunay::RepositionPoint ( PTPoint );
//
//  Redefined in such a way that, after placing DetachedPoint in the 
//  appropriate triangle/edge (and recomputing its error through AddPoint),
//  we also insert such point in PtsErrTree
//

void TRefErrDelaunay::RepositionPoint( PTPoint P )
{
   TRefineDelaunay::RepositionPoint( P );
   PtsErrTree.Insert( P );
}


// -----------------------------------------------------------------------------
//  
//  void TRefErrDelaunay::DetachTriangle()
//  void TRefErrDelaunay::DetachEdge()
//
//  The code of the following functions is nearly equal to their code
//  TRefineDelaunay. The only differene is that points in the PointLists
//  are taken away from PointList and from tree PtsErrTrg.
//  They will be inserted again in the tree and in the PointList of
//  the new edge/triangle in which they fall (by RepositionPoint), after
//  computing their new error (the error has changes since, after the
//  update, the points fall in another triangle/edge).
//


void TRefErrDelaunay::DetachTriangle( PTTriangle T )
{
    PTPoint P = NULL;
    
    TListIterator<PTPoint> PLIter( &T->PointList );
    while( !PLIter.EndOfList() )
    {
       P = PLIter.Current()->object;
       PtsErrTree.Remove( P );
       PLIter.GoNext();
    }

    TRefineDelaunay::DetachTriangle(T);

    //
    // notify MT Tracer that triangle has been deleted 
    //

}


void TRefErrDelaunay::DetachEdge( PTEdge E )
{
    PTPoint P = NULL;
  
   TListIterator<PTPoint> PLIter( &E->PointList );

    while( !PLIter.EndOfList() )
    {
       P = PLIter.Current()->object;
       PtsErrTree.Remove( P );
       PLIter.GoNext();
    }

    TRefineDelaunay::DetachEdge(E);
    
}

