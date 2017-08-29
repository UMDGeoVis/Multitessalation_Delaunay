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
//   file   : refdel.h
//   author : Christian Melchiorre
//
//   Definition of class TRefineDelaunay, subclass of TBuildDelaunay.
//   Unlike its base class, this class builds an initial
//   triangulation by triangulating the convex hull of the input points
//   (instead of taking the first three non-aligned points).
//   Then, it inserts all remaining points.
//   In summary, it redefined method InitialTriangulation() and adds
//   method CalcConvexHull() (and other auxiliary functions).
//   Unlike TBuildDelaunay, it uses an object of class MTTracer for defining
//   the various termination criteria and for tracing the history of the 
//   updates. Such history will be used to build a multiresolution model,
//   the MT.
//   In order to allow using an error-based termination criterion,
//   it adds the namagement of PointLists associated with edges and
//   triangles in the parent class (bucketing technique).
//

#ifndef _REFDEL_H
#define _REFDEL_H

#include "defs.h"
#include "ttriang.h"
#include "builddel.h"
#include "mttracer.h"


class TRefineDelaunay;

class TRefineDelaunay : public TBuildDelaunay
{
   protected:

      //
      // Status variables
      //
      
      // Number of points on the convex hull.
      int nChPts; 

      // Static variable used for sorting the points in CalcConvexHull().
      // WHY STATIC?
      static PTPoint CenterPoint;

      // Points that belonged to the PointLists of triangles/edges
      // of the part of triangulation that has been modified.
      // When we remove an edge/triangle W from the triangulation, the
      // points that were in the PointList of W (points, that are not
      // vertices of the triangulation, falling inside W) need to find
      // their place in the PointLists of the new edges/triangles in
      // which they will fall. Such points are temporarily placed in
      // list DetachedPoints, before performing such repositioning
      // operation.
      // The old triagles/edges of the influence region of a vertex are
      // deleted by function DeleteInfluenceRegion. Such function, for each
      // deleted edge/triangle, moves the points from the PointList of it
      // into list DetachedPoints.
      // After the retriangulation and optimization of the influence region,
      // DetachedPoints is scanned, and each point of it is added to the
      // PointList of the triangle/edge that contains the vertical projection
      TList<PTPoint> DetachedPoints;

      // True if we are in the stage of building the initial triangulation.
      boolean InitialPhase;

      // Used for the various termination criteria and for tracing the
      // history of the updates, in order to build the MT (multi-resolution
      // model).
      PMTTracer MT;

      //
      // methods
      //

      virtual void InitialTriangulation();

      virtual boolean NoMoreUpdates()
        { return( iNextPoint >= nPts || MT->TerminateCondition() ); };

           
      void CalcConvexHull();

      void GrahamSortPoints();
         void HeapSortSift( int, int, PTPoint );
         int GrahamSortFunc( PTPoint, PTPoint, PTPoint );

      virtual void NextPoint();

      virtual void DeleteInfluenceRegion();
        
      virtual void RepositionPoint( PTPoint );
      
      virtual void AddTriangle( PTTriangle );
      
      virtual void DetachTriangle( PTTriangle );
      virtual void DetachEdge( PTEdge );

      virtual void EndTriangulation();
 
      virtual void MT_Initial();
      virtual void MT_KillInterference();
      virtual void MT_AddComponent();
      virtual void MT_Terminate();

      virtual void WriteData( const char * );
      
   public:

      TRefineDelaunay( PMTTracer );

};

typedef class TRefineDelaunay *PTRefineDelaunay;
typedef class TRefineDelaunay &RTRefineDelaunay;

#endif // _REFDEL_H

