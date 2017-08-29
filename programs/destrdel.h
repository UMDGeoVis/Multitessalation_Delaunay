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
//   file   : destrdel.h
//   author : Christian Melchiorre
//
//   Definition of class TDestroyDelaunay, for the decimation of an
//   initial triangulation that has been read from an input file.
//



#ifndef _DESTRDEL_H
#define _DESTRDEL_H

#include <fstream>

#include "tbtree.h"
#include "basedel.h"

class TDecimDelauay;

typedef class TDestroyDelaunay *PTDestroyDelaunay;
typedef class TDestroyDelaunay &RTDestroyDelaunay;


class TDestroyDelaunay : virtual public TDelaunayBase
{
   protected:


     // ------------------------------------------------------------------
     //
     // Status variables
     //
      
     
     //
     // Variables used for determinign the influence region and 
     // retriangulating it.
     //
               
     // Queue containing the edges of the triangulation to be optimized.
     // The retriangulation process initially creates an arbitrary 
     // triangulation of the influence region. Such triangulation may
     // not be a Delaunay one. Then, it is optimized. The optimization
     // process works on queue SwapEdgeQueue. The queue contains the
     // edges that are potentially non-optimal. It is initialized 
     // with the new internal edges of the retriangulation. 
     TDoubleList<PTEdge> SwapEdgeQueue;

     // Set of the removable vertices that satisfy the degree constraint.
     // Such set is sorted on the Error field.
     // It is implemented as a btree in order to make more efficient 
     // the extraction of the vertex of minimum potential error.
     TBTree<PTVertex> ElimVtxTree;

     // Maximum degree of a vertex in order to be removable.
     int KDegree;

     // ------------------------------------------------------------------
     //
     // Methods
     //
      
     virtual void ReadData( const char * );
       virtual void ReadVertices( ifstream & );
       virtual void ReadTriangles( ifstream & );
 

     virtual void InitialTriangulation();
       virtual boolean IsVtxElim( PTVertex );
       virtual boolean OkDegree( int );
       virtual int  CalcDegree( PTVertex );
              
     virtual boolean NoMoreUpdates()
       { return( ElimVtxTree.IsEmpty() ); };

      
     //
     // Update step.
     //
      
     virtual void UpdateStep();
        virtual void NextVertex();
        virtual void RemoveVertex();
      
      
     //
     // Computation of the influence region of a point to be inserted.
     //
       
     virtual boolean CalcInfluenceRegion();
         virtual void InitInflRegnHull();
	 virtual void InitInflRegnInternal();
	 virtual void CalcInflRegnMain();
      
     //
     // Retriangulation of the influence region.
     //
     virtual void RetriangulateInfluenceRegion( boolean );

     // Compute an arbitrary triangulation of the influence region.
     virtual void RIRInitialTrg( TDoubleList<PTEdge> & InflRegnBorder,
                                TDoubleList<PTEdge> & SwapEdgeQueue   );

     // Optimize the arbitrary triangulation of the influence region
     // (i.e., make it a Delaunay triangulation).
     virtual void RIRDelOptTrg();

	 virtual boolean EdgeToSwap( PTEdge );
	 virtual PTEdge SwapEdge( PTEdge & );

         virtual boolean OkTriangle(PTVertex , PTVertex, PTVertex, TDoubleList<PTEdge> &);
 
         virtual void GetQuadBorder( PTEdge, PTEdge&, PTEdge&, PTEdge&, PTEdge& ); 
	 
	 
     //
     // Delete the influence region.
     // Temporarily made public in order to make it visible to function
     // Delete_Constr_InfluenceRegion under the RH52 compiler.
     //

     public:   

         virtual void DeleteInfluenceRegion();   

     virtual boolean ReCheckVertex( PTVertex );
	 
	 
   public:

      TDestroyDelaunay( int );

};

#endif // _DESTRDEL_H
