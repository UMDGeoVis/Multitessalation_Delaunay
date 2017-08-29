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
//   file   : decerrCDT.cpp
//   author : Alessio Calcagno
//
//   Implementation of class TDecErrCDT, sub-class of TDecCDT and
//   TDecErrDelaunay, for decimation of a Delaunay triangulation with an
//   error-based choice of the next vertex to be removed.
//   This class re-implements method InitialTriangulation in this way:
//   before inserting the removable vertices
//   in ElimVtxTree, it sets their Error field with the potential error
//   (function compare for PTVertex used in TBTree uses this field).
//

/*
#include <time.h>

#include "defs.h"
#include "error.h"

#include "ttriang.h"

#include "decCDT.h"
*/
#include "decerrCDT.h"

// ----------------------------------------------------------------------------
// 
//  Constructor of class TDecErrCDT
//


TDecErrCDT::TDecErrCDT( int iK, PMTTracer iMT, int iRecalcError ) 
   : TDecCDT( iK, iMT ),
     TDecErrDelaunay( iK, iMT, iRecalcError ),
     TDecimDelaunay( iK, iMT )
{

  #ifdef DEBUG
      DEBUG << "TDecErrCDT Constructor" << endl;
  #endif 
   /*
   check( (iRecalcError != RECALC_APPROX && iRecalcError != RECALC_EXACT ),
      "TDecErrDelaunay constructor: invalid value for parameter RecalcError" );
      
   RecalcError = iRecalcError;
   */
}


TDecErrCDT::TDecErrCDT( int iK, PMTTracer iMT, boolean EXTActive,
                        boolean ALLOWFeaturesDel, boolean ALLOWChainBrk, int iRecalcError ) 
   : TDecCDT( iK, iMT, EXTActive, ALLOWFeaturesDel, ALLOWChainBrk ),
     TDecErrDelaunay( iK, iMT, iRecalcError ),
     TDecimDelaunay( iK, iMT )
{

  #ifdef DEBUG
      DEBUG << "TDecErrCDT(Options) Constructor" << endl;
  #endif 
   /*
   check( (iRecalcError != RECALC_APPROX && iRecalcError != RECALC_EXACT ),
      "TDecErrDelaunay constructor: invalid value for parameter RecalcError" );
      
   RecalcError = iRecalcError;
   */
}


// -----------------------------------------------------------------------------
//  
//  void TDecErrCDT::InitialTriangulation()
//
//  Similar to the implementation given in TDecRndCDT.
//

void TDecErrCDT::InitialTriangulation()
{
   #ifdef DEBUG
      DEBUG << "\nTDecErrCDT::InitialTriangulation()" << endl;
   #endif
       
   for( int v=0; v<nPts; v++ )
   {

      if( ((PTVertex) Points[v])->VE[0] != NULL )
      // in order to read also files which contain points not present
      // in the triangulation (e.g., output of a decimation): for such
      // points,  only coordinates are initialized (e.g.,  0.0  0.0  0.0 ).
      {
         PTVertex V = (PTVertex)(Points[v]);
                
         if ( IsVtxElim( V ) && OkDegree( CalcDegree(V) ) && OkConstrDegree( V ) )
         {
            RecalcVertexError( V );
            ElimVtxTree.Insert( V );
         }
       
      } // end ...if(  Points[v]->VE[0] != NULL
   } // end ...for(v)

   // Since TDecCDT inherits from both TDecimDelaunay and TRefCDTDelaunay
   // it is necessary to specify the base class

   TDecimDelaunay::MT_Initial();

   TDecimDelaunay::InitialPhase = FALSE;

}


// ----------------------------------------------------------------------------
// 
//  void TDecErrCDT::ReCheckVertex( PTVertex V )
//
//  Re-implemented in such a way that, before inserting a vertex in tree
//  ElimVtxTree in TDestroyDelaunay::DeleteInfluenceRegion(), the error
//  of such a vertex is re-computed.
//

boolean TDecErrCDT::ReCheckVertex( PTVertex V )
{

   if ( TDecCDT::ReCheckVertex(V) )
   {
      RecalcVertexError( V );
      return(TRUE);
   }
   else
      return(FALSE);
}
