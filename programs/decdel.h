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
//   file   : decdel.h
//   author : Christian Melchiorre
//
//   Class TDecimDelaunay inherits from class TDestroyDelaunay and adds
//   an object of type MTTracer. Such object defines several termination
//   criteria and allows tracing the history of updates (the history will
//   be used to build a multi-resolution model, the MT).
//   In order to allow using an error-based termination criterion,
//   it adds the namagement of PointLists associated with edges and 
//   triangles in the parent class (bucketing technique).
//


#ifndef _DECDEL_H
#define _DECDEL_H

#include "defs.h"
#include "ttriang.h"
#include "destrdel.h"
#include "mttracer.h"


class TDecimDelaunay;

typedef class TDecimDelaunay *PTDecimDelaunay;
typedef class TDecimDelaunay &RTDecimDelaunay;

class TDecimDelaunay : public TDestroyDelaunay
{
   protected:
   
      //
      // Status variables
      //

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
      // of such point.
      TList<PTPoint> DetachedPoints;
      
      // True if we are in the stage of building the initial triangulation.
      boolean InitialPhase;

      // Used for the various termination criteria and for tracing the
      // history of the updates, in order to build the MT (multi-resolution
      // model).
      PMTTracer MT;

      //
      // Methods
      //
      
      virtual void InitialTriangulation();

      virtual boolean NoMoreUpdates()
        { return( ElimVtxTree.IsEmpty() || MT->TerminateCondition() ); };

      virtual void RepositionPoint( PTPoint );
      
      virtual void AddTriangle( PTTriangle );
      virtual void DetachTriangle( PTTriangle );
      virtual void DetachEdge( PTEdge );
 
      virtual void DeleteInfluenceRegion();

      virtual void EndTriangulation();
 
      virtual void MT_Initial();
      virtual void MT_KillInterference();
      virtual void MT_AddComponent();
      virtual void MT_Terminate();

      virtual void WriteData( const char * );
          
   public:
   
      TDecimDelaunay( int, PMTTracer );

};

#endif // _DECDEL_H
