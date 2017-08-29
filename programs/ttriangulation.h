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
//   file   : ttriangulation.h
//   author : Christian Melchiorre
//
//   Definition of a generic class for triangulations.
//   Used as an abstract base class for building triangulations according
//   to various criteria.
//


#ifndef _TTRIANGULATION_H
#define _TTRIANGULATION_H

#include "defs.h"
#include "ttriang.h"
#include "stepbystep.h"

// -----------------------------------------------------------------------------
//
//   class TTriangulation 
//
//   Class definition
//


class TTriangulation;

typedef class TTriangulation *PTTriangulation;
typedef class TTriangulation &RTTriangulation;


// Definition of a generic class for triangulations.
// Used as an abstract base class for building triangulations according
// to various criteria.
// For a description of the data structure that implements the
// triangulation, see file ttriang.h 

class TTriangulation
{

   friend class StepByStep;
	
   public:

     TTriangulation();


     // Main procedure, it performs the loop of uodates to the triangulation.
     virtual void BuildTriangulation( const char *, const char * );
     

   protected:

     //
     // Status variables
     //
      
     // Array containing the initial indices of the points in array Points.
     // Used during refinement, it stores the correspondence to the
     // initial position of the points in array Points.
     int *OrderInitial;
      
     // Array of pointers to the input points in the order in which they
     // have been read.
     // During refinement, function InitialTriangulation(), which creates
     // the initial triangulation of the convex hull, resorts the
     // array Points in such a way that the first elements point to
     // the vertices of the convex hull, sorted counterclockwise.
     // After the execution of TRefineDelaunay::InitialTriangulation() 
     // the pointers to the points of the convex hull are in positions
     // [0...TRefineDelaunay::nChPts-1] dell'array Points, while the
     // other points (not vertices of the convex hull) are pointed by the
     // remaining positions [nChPts...nPts-1].
     // During refinement with point insertion in random order, function
     // RandomSortPoints() sorts positions [nChPts...nPts-1] based on a
     // random key.
     PTPoint *Points;
      
     // Number of input points, and size of array Points.
     // It is not the number of vertices currently present in the
     // triangulation.
     int nPts;


     // Point to be inserted in the triangulation as a new vertex.
     // In generic and random refinement, points are inserted in the same
     // order in which they are sorted in array Points; thus, at each step
     // VertexToIns is the point of index iNextPoint in array Points.
     // In error-driven refinement, VertexToIns is extracted from the
     // sorted set of points to be inserted: TRefErrDelaunay::PtsErrTree,
     // i.e., the point of maximum (real) error.
     PTVertex VertexToIns;

     // Vertex to be removed from the triangulation.
     // It is extracted from the set of removable vertices
     // TDestroyDelaunay::ElimVtxTree, i.e., it is the point of
     // minimum (real) error.
     // For edge-collapse, it is the first endpoint of the edge to
     // be collapsed.
     PTVertex VertexToRemove;

     // Index in array Points which divides the already-iserted points
     // from the still-to-be-inserted ones.
     // Used to scan the array Points during point insertion:
     //   for all i < iNextPoint: Points[i] refer to vertices
     //                     (points already inserted in the triangulation);
     //   for all i >= iNextPoint: Points[i] refer to points
     //                     (not yet inserted in the triangulation);
     // It is the entry of array Points in which to put vertex VertexToIns
     // when it is created.
     // In generic and random refinement,  points are inserted in the same
     // order in which they are sorted in array Points; thus, it also
     // denotes the next point to be inserted in the triangulation.
     int iNextPoint;


     // Entry point to the set of triangles of the triangulation.
     // It is an arbitrary triangle of the triangulation, from which we
     // start visiting the triangulation. 
     // Usually, it is the last triangle added to the triangulation.
     PTTriangle FirstTriangle;

     // Number of triangles present in the current triangulation.
     int nTrg;

     //
     // Functions and variables used in the point location process.
     //
      
     // Locates the triangle or edge which contains point p
     virtual void PointLocation( PTPoint p );

     PTVertex PLVertex;	// questa e' un'aggiunta di Paola!
     int PLLocation;
     PTTriangle PLTriangle;
     PTEdge PLEdge;


     //
     // Methods
     //

     //
     // Read/write data from/to disk.
     //
     
     virtual void ReadData( const char * ) = 0;
     virtual void WriteData( const char * );
     
     // Return the vertices, triangles and (if present) constraint edges
     // in indexed format, by putting them into three arrays that are passed
     // as parameters to the function.
     // It also returns the number of vertices, triangles and constraint edges.
     // Arrays are de-allocated (if not null) and re-allocated inside
     // the function.
     virtual void ConvertData(int *vNum, int *tNum, int *eNum,
                              float **vData, int **tData, int **eData);
     
     // Preliminary work before starting the updates on the initial
     // triangulation.
     // In refinement, it computes an initial triangulation
     // (triangulation of the vertices of the convex hull),
     // adjusts the 'array Points; and in TRefErrDelaunay it puts 
     // into PtsErrTree the points still to be inserted.
     // In decimation, it inserts into ElimVtxTree the removable vertices
     // which satisfy the bound on the degree.
     virtual void InitialTriangulation() = 0;
     
     
     
     // Boolean function, it checks if further updates are needed to
     // obtain the desired precision or the desired number of updates.
     virtual boolean NoMoreUpdates() = 0;
     
     
     // Function that performs one update step.
     // In refinement, it selects a vertex to be inserted and inserts it
     // in the triangulation.
     // In decimation, it selects a vertex to be removed and 
     // removed it from  the triangulation. 
     virtual void UpdateStep() = 0;
     

      //
      // Add and remove one entity (vertex, edge, triangle)
      //
      
      virtual void AddTriangle( PTTriangle T )	
      { 
          FirstTriangle = T; nTrg++; 
          #ifdef DEBUG
	    DEBUG << endl << "ATTACHED : T" << T->TID << "(ntrg = " << nTrg << ")" << endl;
	  #endif	    
      };
 
      virtual void DetachVertex( PTVertex );
      virtual void DetachEdge( PTEdge );
      virtual void DetachTriangle( PTTriangle );
          
     
      //
      // Function performing various work at the end of the construction
      // process. By default, it does nothing.
      virtual void EndTriangulation() {};
      virtual void PrepareToEnd() {};
};

#endif // _TTRIANGULATION_H
