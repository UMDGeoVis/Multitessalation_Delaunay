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

// ---------------------------------------------------------------------------------
//
//  file   : decrndb.h
//  author : Christian Melchiorre
//
//  Definition of class TDecRndDeBerg, sub-class of TDecRndDelaunay, for
//  decimation of a Delaunay triangulation with removal of independent 
//  vertices and a random choice of vertices to be removed.
//  To this aim, procedure DeBerg::SelectVertices of class DeBerg is used
//  (files: deberg.h/deberg.cpp).
//


#ifndef _DECRNDDEB_H
#define _DECRNDDEB_H


#include "defs.h"
#include "ttriang.h"

#include "tdoublelist.h"
#include "tbtree.h"

#include "decdel.h"
#include "decrnddel.h"


class TDecRndDeBerg;

typedef class TDecRndDeBerg *PTDecRndDeBerg;
typedef class TDecRndDeBerg &RTDecRndDeBerg;


class TDecRndDeBerg : public TDecRndDelaunay
{
   protected:

      TDoubleList<PTVertex> DeBergVertices;
      
   public:
   
      virtual boolean NoMoreUpdates()
        { return( ( ElimVtxTree.IsEmpty() && DeBergVertices.IsEmpty() )
	          || MT->TerminateCondition() ); };

      virtual void NextVertex();		  
		  
      TDecRndDeBerg( int iK, PMTTracer iMT ) 
          : TDecRndDelaunay( iK, iMT ), DeBergVertices(), TDecimDelaunay( iK, iMT ) {};

};


#endif // _DECRNDDEB_H
