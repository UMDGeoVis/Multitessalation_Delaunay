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
//  file   : decerrdel.cpp
//  author : Christian Melchiorre
//
//  Definition of class TDecErrDelaunay, sub-class of TDecimDelaunay,
//  for decimation of a Delaunay triangulation with an error-based choice 
//  of the next vertex to be removed.
//


#ifndef _DECERRDEL_H
#define _DECERRDEL_H


#include "defs.h"
#include "geom.h"
#include "tdoublelist.h"
#include "ttriang.h"
#include "mttracer.h"
#include "decdel.h"


class TDecErrDelaunay;

typedef class TDecErrDelaunay *PTDecErrDelaunay;
typedef class TDecErrDelaunay &RTDecErrDelaunay;


class TDecErrDelaunay : virtual public TDecimDelaunay
{

   protected:

      //
      // Attributes
      //
      
      int RecalcError;
      
      //
      // Methods
      //

      virtual void InitialTriangulation();

      virtual boolean ReCheckVertex( PTVertex );

      virtual void RecalcVertexError( PTVertex );
      virtual void RecalcVertexErrorApprox( PTVertex );
      virtual void RecalcVertexErrorExact( PTVertex );
      virtual boolean OkTriangle( PTVertex, PTVertex, PTVertex, TDoubleList<PTVertex>& );

   public:

     TDecErrDelaunay( int, PMTTracer, int iRecalcError = RECALC_APPROX );

};


#endif // _DECERRDEL_H
