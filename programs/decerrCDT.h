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

// ----------------------------------------------------------------------------
//
//  file   : decerrCDT.cpp
//  author : Alessio Calcagno
//
//  Definition of class TDecErrCDT, sub-class of TDecCDT and
//   TDecErrDelaunay, for decimation of a Delaunay triangulation with an
//   error-based choice of the next vertex to be removed.
//


#ifndef _DECERRCDT_H
#define _DECERRCDT_H

/*
#include "defs.h"
#include "geom.h"

#include "tdoublelist.h"

#include "ttriang.h"
#include "mttracer.h"

#include "decdel.h"
*/
#include "decCDT.h"
#include "decerrdel.h"


class TDecErrCDT;

typedef class TDecErrCDT *PTDecErrCDT;
typedef class TDecErrCDT &RTDecErrCDT;


// Decimation of a constrained Delaunay triangulation with an error-driven
// choice of the vertex to be removed. The initial CDT is read from input.
// Class TDecErrCDT, sub-class of TDecCDT and TDecErrDelaunay.
// The only modificaton  with respect to the super-classes is in the
// re-implementation of method InitialTriangulation: before inserting
// ramovable vertices in ElimVtxTree, set their Error field to the potential
// error (field Error is needed for function compare of PTVertex used in the
// binary tree TBTree).

class TDecErrCDT : public TDecCDT, public TDecErrDelaunay
{

   protected:
/*
      //
      // Attributes.
      //
      
      int RecalcError;
      
      //
      // Moethods.
      //
*/
      virtual void InitialTriangulation();

      virtual boolean ReCheckVertex( PTVertex );
/*
      virtual void RecalcVertexError( PTVertex );
      virtual void RecalcVertexErrorApprox( PTVertex );
      virtual void RecalcVertexErrorExact( PTVertex );
      virtual boolean OkTriangle( PTVertex, PTVertex, PTVertex, TDoubleList<PTVertex>& );   
*/
   public:

     TDecErrCDT( int, PMTTracer, int iRecalcError = RECALC_APPROX );
     TDecErrCDT( int iK, PMTTracer iMT, boolean EXTActive,
                 boolean ALLOWFeaturesDel, boolean ALLOWChainBrk, int iRecalcError = RECALC_APPROX );

};


#endif // _DECERRCDT_H
