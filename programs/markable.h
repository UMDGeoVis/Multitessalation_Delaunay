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
//  file   : markable.h
//  author : Christian Melchiorre
//
//  Implementation of class Markable used to define the behaviour of
//  objects which can be marked.
//



#ifndef _MARKABLE_H
#define _MARKABLE_H

#include "defs.h"

// The data type used to mark objects.
// Remark: the mark field is used as a bit-vector. We have to make sure
// that the number of bits is sufficient with respect to the number of
// different marking values. For instance, it is possible to define it
// of "long" type without any other modifications to the implementation 
// of class Markable or any code using such class.

typedef int MARKTYPE;


//
// Some constants used to mark edges, triangles...
//
// REMARK: If you add or modify any constant, ALWAYS remember to update
// constants FIRST_MARK and LAST_MARK, in such a way that they correspond
// to the first and last defined constant.
//

// This is the value for a non-marked entity
const MARKTYPE NULLMARK     = 0x00;

// This mark is used during the computation of the influence region,
// to mark entities that must be deleted from the triangulation
const MARKTYPE TO_DELETE    = 0x01;         // bit 0

// This mark is used during the computation of the influence region,
// to mark an edge belonging to the influence region
const MARKTYPE INFL_BORDER  = 0x02;         // bit 1

// This mark is used during the traversal of the whole triangulation,
// or a part of the triangulation, to mark entities that have already
// been visited
const MARKTYPE VISITED      = 0x04;         // bit 2

// This mark is used during the retriangulation or Delaunay optimization,
// to mark new triangles, and means that such triangles must be added in
// MT_Hist.
const MARKTYPE NEW_TRIANGLE = 0x08;         // bit 3

// This mark is used by function MT_KillInterference to mark triangles
// previously marked as TO_DELETE, and means that such triangles have been
// registered as deleted by MT_Hist
const MARKTYPE MT_DELETED   = 0x10;         // bit 4

// Other constants used by function TDecimDelaunay::RetriangulateInflRegn()
const MARKTYPE INFL_BORDER_AUX = 0x20;      // bit 5

// This mark is used during optimization of the retriangulation of the
// influence region, to mark edges that need to be checked if they are
// optimal
const MARKTYPE SWAP_EDGE_QUEUE = 0x40;      // bit 6

// This mark is used in a constrained triangulation, to mark edges
// that are constraints
const MARKTYPE CONSTRAINED = 0x80;          // bit 7

// This mark is used to mark new edges inserted in the triangulation,
// during the retriangulation after the removal of a vertex with 
// some constraint edges indicent in it.
const MARKTYPE NEW_EDGE = 0x100;          // bit 8

// This mark is used to mark triangles and vertices that have already been
// visited, during the recomputation of degree, constr-degree, potential
// error and possibility of adding in ElimVtxTree of the vertices of the
// triangulation that are involved in the retriangulation after the removal
// of a vertex with some constraint edges indicent in it.
const MARKTYPE RECHECKED = 0x200;          // bit 9

// This mark is used to mark original edges that have a copy of themselves 
// in TDecCDT::OrigEdgList. The copy contains all edge attributes as they
// were before starting the retriangulation process.
const MARKTYPE COPIED = 0x400;          // bit 10

// REMARK: If you add or modify any constant, ALWAYS remember to update
// constants FIRST_MARK and LAST_MARK, in such a way that they correspond
// to the first and last defined constant.

// These costants must be the same as the first and last mark, they are
// used in method TEdge::Match( PTEdge )
const MARKTYPE FIRST_MARK = 0x01; 
const MARKTYPE LAST_MARK  = 0x400;

//
// Class Markable allows associating one or more marks with an object
//

class Markable
{
   private:

      // Attribute containing all marking information (bit vector).
      MARKTYPE mark;
      
   public:

      Markable() : mark( NULLMARK ) {};

      void Mark( MARKTYPE MarkValue );
      boolean Marked( MARKTYPE MarkValue );
      void UnMark( MARKTYPE MarkValue );
      void MarkReset();

      MARKTYPE MarkValue();    
};

#endif // _MARKABLE_H
