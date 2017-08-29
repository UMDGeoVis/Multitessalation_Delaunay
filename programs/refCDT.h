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
//   file   : refCDT.h
//   author : Alessio Calcagno
//

#ifndef _REFCDT_H
#define _REFCDT_H

#include "utils.h"
#include "refdel.h"
#include "decCDT.h"


class TRefCDT;

typedef class TRefCDT *PTRefCDT;
typedef class TRefCDT &RTRefCDT;

class TRefCDT : public TRefineDelaunay, public TDecCDT
{
   protected:

      // Array of input constraints
      // Used for CDT refinement. It contains pairs of point indices
      // in array Points.
      int *constr_array;

      // Number of input constraints
      int nConstr;

/**********************                   *********************************/

      // Re-defined functions

      virtual void Add_Constrain();
//      virtual boolean ReCheckVertex( PTVertex V )  { return TRUE; }
      virtual void Delete_Constr_InfluenceRegion( PTEdge newConstr );

      // da RefCDTDelaunay
      virtual void ReadData(const char *infname);
      virtual void InitialTriangulation();
      virtual void NextPoint();
      virtual void BuildTriangulation( const char *infname, const char *outfname );

      // Redefine ttriangulation::WriteData to write constraints
      virtual void WriteData( const char * outfname );
      // Redefine ttriangulation::ConvertData to convert constraints
      virtual void ConvertData(int *vNum, int *tNum, int *eNum,
                               float **vData, int **tData, int **eData);
                    
      // Redefined to remove ambiguity
      // (functions inherited from both parent classes)
      virtual boolean NoMoreUpdates() { return TRefineDelaunay::NoMoreUpdates(); }
      virtual void UpdateStep() { TRefineDelaunay::UpdateStep(); }
      virtual void AddTriangle( PTTriangle T ) { TRefineDelaunay::AddTriangle( T ); }
      virtual void DetachEdge( PTEdge E ) { TRefineDelaunay::DetachEdge( E ); }
      virtual void DetachTriangle( PTTriangle T ) { TRefineDelaunay::DetachTriangle( T ); }
      virtual void DeleteInfluenceRegion() { TRefineDelaunay::DeleteInfluenceRegion(); }
      virtual void EndTriangulation() { TRefineDelaunay::EndTriangulation(); }
      virtual void MT_AddComponent() { TRefineDelaunay::MT_AddComponent(); }
      virtual void MT_Initial() { TRefineDelaunay::MT_Initial(); }




/**********************                  *********************************/

      // New functions
      
      virtual void PrepareToEnd();

   public:
  
       TRefCDT( PMTTracer iMT ) : TRefineDelaunay( iMT ), TDelaunayBase(),
                                  TDecimDelaunay( 0, iMT ), TDecCDT( 0, iMT )
       {
	    #ifdef DEBUG
		 DEBUG << "TRefCDT Constructor" << endl;
  	    #endif
	    constr_array = NULL;
            nConstr = 0;
       };
};

#endif // _REFCDT_H
