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
//   file   : ttriang.h
//   author : Christian Melchiorre
//
//   Definition of classes for the main data structure.
//   The data structure used to encode the triangulation is the symmetric
//   data structure for triangulations.
//   The basic triangulation elements are vertices, edges and triangles.
//   Each of these elements is individually stored and is related with
//   other elements by means of topological relations.
//   This file defines classes that implement the basic elements.
//   It also defines a class which denotes a point in 3D space, used to
//   store input points that have not yet been inserted in the triangulation.


#ifndef _TTRIANG_H
#define _TTRIANG_H

using namespace std;


#include <iostream>
#include <fstream>
#include "defs.h"
#include "error.h"
#include "markable.h"
#include "tlist.h"

#ifdef MT_TRACER
#include "fieldmt.h"
#endif


// -----------------------------------------------------------------------------
//
//  Classes defined in this file
//


class TPoint;
class TVertex;
class TEdge;
class TTriangle;



// -----------------------------------------------------------------------------
//
//   pointers and references to the defined classes
//



typedef class TPoint *PTPoint;
typedef class TPoint &RTPoint;

typedef class TVertex *PTVertex;
typedef class TVertex &RTVertex;

typedef class TEdge *PTEdge;
typedef class TEdge &RTEdge;

typedef class TTriangle *PTTriangle;
typedef class TTriangle &RTTriangle;



// -----------------------------------------------------------------------------
//
//   class TPoint
//
//   An input point that has not yet been inserted in the triangulation.
//   It contains the 3D coordinates of the point and its error, but no
//   relation with edges or triangles.
//   All points of the dataset are read from the input file and stored
//   individually as objects of class TPoint.
//   When such points are inserted in the triangulation, and thus they
//   become vertices of the current triangulation, they are stored as
//   objects of class TVertex.

class TPoint
{
  public:
   
    // Unique identifier of this point.
    int PID;

    // Point coordinates.
    double x,y,z;
    
    // Point error: vertical distance from the triangulation entity
    // in which it lies / it whould lie.
    // It has a different meaning if the point is a vertex of the 
    // current triangulation or not (i.e., if the point has not yet been
    // inserted during refinement or has already been removed during 
    // decimation. Thus, it may represent:
    // - real error: if the point is not a vertex of the triangulation 
    // - potential error: error that the point would have if it were
    //   inserted / removed as vertex during refinement / decimation.
    double Error;

    // Default constructor.
    TPoint( double xi=0, double yi=0, double zi = 0 )
      : PID(-1), x(xi), y(yi), z(zi), Error(0.0){};
  

    // Compare the coordinates of two points and check if they are equal.    
    boolean Equals( PTPoint p );

    // Version that considers only x any y coordinates.
    boolean Equalsxy( PTPoint p );
    
    // Input/output operators.
    friend ostream& operator<< ( ostream&, RTPoint p );
    friend istream& operator>> ( istream&, RTPoint p );    
};

int compare( PTPoint, PTPoint );


// -----------------------------------------------------------------------------

//
// class TVertex
// 
//  A vertex of the triangulation.
//  It implements a point that has been inserted in the triangulation.
//  It inherits coordinates and error from the point, and, in addition, it
//  contains topological relations, the unique identifier of the vertex 
//  within the MT and the number of constraint edges incident in it.


class TVertex : public TPoint, public Markable
{
      
     public:

     // Number of constraint edges incident in this vertex.
     int nIncConstr;

     // Next free identifier (static variable).
     static int NextVID;
   
     // Unique vertex identifier.
     int VID;
   
     // Partial Vertex-Edge relation.
     // For internal vertices, VE[0] and VE[1] are two arbitrary edges 
     // incident in this vertex. For vertices lying on the convex hull,
     // VE[0] and VE[1] are the two incident edges that lie on the
     // convex hull.
     PTEdge VE[2];

     // List of edges incident in this vertex (total Vertex-Edge relation)
     // and function to compute such list.
     TList<PTEdge> EdgeList;
     void GetVE();

     // For MT tracer.

#ifdef MT_TRACER
      MT_INDEX MTIdx;
#endif
     
     //
     // Constructors
     //
      
     TVertex( double, double, double );
      
     TVertex( PTPoint );
     TVertex( RTPoint );

     //
     // Destructor
     //
      
     ~TVertex();
      
     // output operator
     friend ostream& operator<< ( ostream&, RTVertex );
      
};


 

// -----------------------------------------------------------------------------
//
//   class TEdge
//
//   An edge of the triangulation.
//


class TEdge : public Markable
{
    public:
   
      // Next free vertex identifier (static variable).      
      static int NextEID;
       
      // Unique identifier of this edge.
      int EID;
   
      // Edge-Vertex relation.
      PTVertex EV[2];
      
      // Edge-Triangle relation.      
      PTTriangle ET[2];
//----------------------------------------------------------------------------
	  //EDGE COLLAPSE
      double Error;
      // Vertex towards which this edge must be collapsed in case
      // it is on the domain boundary (NULL otherwise).
      PTVertex ToVertex;
//----------------------------------------------------------------------------
      //
      // Construstors
      //
      
      TEdge( PTVertex, PTVertex );

      TEdge( TEdge &);

      TEdge();

      // Is this edge on the convex hull?
      boolean OnConvexHull()
      {
	 // cerr << "edge: " << this->EID << endl;
         return (
	    (ET[0]==NULL && ET[1] != NULL)
	     ||
	    (ET[1]==NULL && ET[0] != NULL)
	 );
      };
      
      // Return TRUE if this edge is matching with filter or
      // if they point to the same edge object.
      boolean Match( PTEdge filter );

      boolean OldMatch( PTEdge filter );

      // Points lying on this edge.
      // Points that do not belong to the current triangulation (they
      // have not yet been inserted, or have already been removed) and
      // their vertical projection falls on this edge.
      // The head of this list contains the point with maximum error.
      TList<PTPoint> PointList;

      
      // Insert the given point into the PointList of this edge.
      void AddPoint( PTPoint );
      

      //
      // Destructor
      //
      
      ~TEdge();
      
      friend ostream& operator<< ( ostream&, RTEdge );
 
};



// -----------------------------------------------------------------------------
//
//   class TTriangle
//
//   A triangle of a triangulation
//


class TTriangle : public Markable
{

   public:
   
      // Next free identifier (static variable).
      static int NextTID;

      // Unique identifier of this triangle.      
      int TID;

      // Coordinates x, y and radius of the circum-circle of this triangle.
      double InCircleX, InCircleY, InCircleRad;
   

      // Triangle-Edge relation.
      // The three edges that bound this triangle must be stored
      // in counterclockwise order.
      PTEdge TE[3];
      
      // For the MT tracer.
#ifdef MT_TRACER
      MT_INDEX MTIdx;
#endif

      //
      // ConstruCtor
      // 
      
      TTriangle( PTEdge, PTEdge, PTEdge );
      
      
      // Compute Triangle-Vertex relation. 
      // The computation exploits the Triangle-Edge and Edge-Vertex
      // adjacency relations. The three vertex pointers are returned
      // (by reference, as parameters) sorted in counterclockwise order.      
      void GetTV( PTVertex &v0, PTVertex &v1, PTVertex &v2 );
      
      // Compute Triangle-Triangle relation.
      // The computation exploits the Triangle-Edge and Edge-Triangle
      // adjacency relations. As for GetTV(), the three triangle pointers
      // are returned (by reference, as parameters) sorted in
      // counterclockwise order.
      void GetTT( PTTriangle &t0, PTTriangle &t1, PTTriangle &t2 );
      PTTriangle GetTT( PTEdge e );
      PTTriangle GetTT( int i );
      
      
      // Compute the coordinates and the radius of the circum-circle.
      void CalcCircle();

      // Check whether a point is inside the circle or not.
      boolean InCircle( PTPoint );
      

      // Points lying in the triangle.
      // Points that do not belong to the current triangulation (they
      // have not yet been inserted, or have already been removed) and
      // their vertical projection falls inside the triangle.
      // The head of this list contains the point with maximum error.
      TList<PTPoint> PointList;
      
      void AddPoint( PTPoint );

      double GetError();
      
      //
      // Destructor
      //
      
      ~TTriangle();
      
      // Output operator.
      friend ostream& operator<< ( ostream&, RTTriangle );
  
};

#endif // _TTRIANG_H
