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
//  file   : decrnddbCDT.h
//  author : Alessio Calcagno
//
//  Definition of class TDecRndDeBergCDT, sub-class of TDecRndCDT
//  and TDecRndDeBerg, for
//  the decimation of a constrained Delaunay triangulation (CDT) with
//  removal of independent vertices and a random choice of vertices to be
//  removed.
//  To this aim, procedure DeBerg::SelectVertices of class DeBerg is used
//  (files: deberg.h/deberg.cpp).
//



#ifndef _DECRNDDEBCDT_H
#define _DECRNDDEBCDT_H


#include "defs.h"
#include "ttriang.h"

#include "tdoublelist.h"
#include "tbtree.h"

#include "decdel.h"
#include "decrndCDT.h"


class TDecRndDeBergCDT;

typedef class TDecRndDeBergCDT *PTDecRndDeBergCDT;
typedef class TDecRndDeBergCDT &RTDecRndDeBergCDT;


#include "decrnddb.h"
class TDecRndDeBergCDT : public TDecRndCDT, public TDecRndDeBerg
{
   public:
      TDecRndDeBergCDT( int iK, PMTTracer iMT ) 
	      : TDecRndCDT( iK, iMT ), TDecRndDeBerg( iK, iMT ), TDecimDelaunay( iK, iMT )
	{
		// #ifdef DEBUG
		 cout << "TDecRndDeBergCDT Constructor" << endl;
		// #endif
	};

       TDecRndDeBergCDT( int iK, PMTTracer iMT, boolean EXTActive, boolean ALLOWFeaturesDel, boolean ALLOWChainBrk  )
          : TDecRndCDT( iK, iMT, EXTActive, ALLOWFeaturesDel, ALLOWChainBrk ),
            TDecRndDeBerg( iK, iMT ),
            TDecimDelaunay( iK, iMT )
        {
		// #ifdef DEBUG
		 cout << "TDecRndDeBergCDT(Options) Constructor" << endl;
		// #endif
	};

	protected:
		virtual void InitialTriangulation() { TDecRndCDT::InitialTriangulation(); };
		virtual boolean NoMoreUpdates() { return TDecRndDeBerg::NoMoreUpdates(); }
};

/*
class TDecRndDeBergCDT : public TDecRndCDT
{
   protected:

      TDoubleList<PTVertex> DeBergVertices;
      
   public:
   
      virtual boolean NoMoreUpdates()
        { return( ( ElimVtxTree.IsEmpty() && DeBergVertices.IsEmpty() )
	          || MT->TerminateCondition() ); };

      virtual void NextVertex();		  
		  
      TDecRndDeBergCDT( int iK, PMTTracer iMT ) 
	 : TDecRndCDT( iK, iMT ), DeBergVertices()
	 {
	     #ifdef DEBUG
		 DEBUG << "TDecRndDeBergCDT Constructor" << endl;
	     #endif
	 };

};
*/

#endif // _DECRNDDEBCDT_H
