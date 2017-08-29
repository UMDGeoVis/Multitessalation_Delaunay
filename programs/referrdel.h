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
//   file   : referrdel.h
//   author : Christian Melchiorre
//
//   Definition of class TRefErrDelaunay, sub-class of TRefineDelaunay,
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

#ifndef _REFERRDEL_H
#define _REFERRDEL_H


#include "defs.h"
#include "tbtree.h"
#include "ttriang.h"

#include "builddel.h"
#include "refdel.h"
#include "mttracer.h"


class TRefErrDelaunay : public TRefineDelaunay
{
   protected:
   
       // Set of points not yet inserted in the triangulation, sorted based
       // on their Error field.
       // Implemented as a balanced binary tree in order to extract the
       // point of maximum error efficielntly.
       TBTree<PTPoint> PtsErrTree;
       
       virtual void InitialTriangulation();

       virtual void NextPoint();

       virtual void AddAllPointsToTree();

       virtual void RepositionPoint( PTPoint );

       virtual void DetachTriangle( PTTriangle );
       virtual void DetachEdge( PTEdge );
       
   public:
   
       TRefErrDelaunay( PMTTracer );

};


#endif // _REFERRDEL_H
