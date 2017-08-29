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
//   file   : decCDT.h
//   author : Alessio Calcagno
//
//   Definition of Class TDecCDT, sub-class of TDecimDelaunay, for the
//   decimation of a constrained Delaunay triangulation (CDT) read from
//   file.
//   It re-implements methods ReadData, WriteData, InitialTriangulation...
//

#ifndef _DECCDT_H
#define _DECCDT_H

#include "decdel.h"


// Macros for the result of function
// TRA_intersect( PTEdge e1, PTEdge e2 ): intersection test 
// between two edges

#define NO_INTER      0   /* no intersection */
#define ONLY_1CV      1   /* intersection just in a vertex */
#define NOT_PROPER    2   /* a vertex of one falls inside the other one */
#define PROPER_INTER  3   /* proper intersection */
#define UPON_NO_CV    4   /* partially coincident without common vertices */
#define UPON_AND_1CV  5   /* partially coincident with one common vertex */
#define UPON_AND_2CV  6   /* coincident edges */

// Macros for the result of function InitCalcEdgeWideInflRegn

#define EDGE_ALREADY_EXISTS 1
#define QUAD_REGION         2

// Macros for the result of function IsEdgeAdmissible

#define  NOT_EXISTENT       3
#define  EXIST              2
#define  EXIST_CONSTRAINED  1
#define  NOT_ADMISSIBLE     0

class TDecCDT;

typedef class TDecCDT *PTDecCDT;
typedef class TDecCDT &RTDecCDT;

// PF1b2PTEdge is the type pointer to a function that takes two 
// PTEdge's and return a boolean value
typedef boolean ( * PF1b2PTEdge )( PTEdge, PTEdge );


//
// Base class for decimation of a Constrained Delaunay Triangulation (CDT)
// read from file. Sub-class of TDecimDelaunay.
//

class TDecCDT : virtual public TDecimDelaunay
{
   protected:

      // Number of constraint edges incident into vertex VertexToRemove 
      //when it is selected.
      // Set by TDecCDT::RemoveVertex() to the number of constraint 
      // edges incident in VertexToRemove just after its selection in
      // function NextVertex() and never updated.
      // Thus, it stores ConstrDegree(VetexToRemove) before any modification
      // of the triangulation for the removal of VertexToRemove.
      // Used to drive the execution of the decimation steps:
      // - if nIncConstr == 0  apply functions for normal decimation;
      // - if nIncConstr == 1  apply functions that manage extension of
      //   optimization beyond the influence region of VertexToRemove;
      // - if nIncConstr == 2  as previous case, in addition activate
      //   function Tra_Add_Constraint that adds a constraint;
      // - otherwise       error.
      int nIncConstr; 

      // True if we want to remove vertex VertexToRemove with incident 
      // constraints by using functions that manage extension of
      // optimization beyond the influence region of VertexToRemove.
      boolean ExtActive;

      // True if we allow eliminating features, i.e., eliminating
      // vertices with constrDegree == 1
      boolean AllowFeaturesDel;
                                
      // True if we allow breaking closed chains of constraints,
      // i.e., if we allow eliminating vertices taking part of 
      // a chain of three edges.
      boolean AllowChainBrk;                   

      // Right and left boundary w.r.t. to the edge directed from v0 to v1.
      // Used during the insertion of a constraint.
      // Set by function CalcEdgInflRegn. They store, in CCW order,
      // the boundary edges of the regions of influence of an edge.
      TDoubleList<PTEdge> Right_Border, Left_Border;
      
      // List containing a copy of the original edges involved in the 
      // optimization.
      // Such copies contain tha attributes of the original edges as
      // they were before optimization.
      // When the copy of an edge is inserted in this list, the edge is
      // marked as COPIED. When the copy is deleted from the list, the
      // COPIED mark is removed from the edge.
      // The list is initialized just after the computation of the
      // influence region by TDedCDT::ExtRIRPrepare, with the edges of
      // InflRegnBorder, that may be later involved in the optimization,
      // and be removed from the triangulation.
      // It is not correct to eliminate the original edges directly:
      // we must reconstruct the original topological relations of such
      // edges and mark them as TO_DELETE.
      // The aim is obtaining a consistent region formed by original
      // entities (marked as TO_DELETE) and a region formed by new and
      // final entities (marked as NEW_TRIANGLE/EDGE) that occupy the
      // same area.
      // When an original edge is marked as TO_DELETE during the 
      // optimization, certainly it is no more considered for the
      // optimization, thus we can remove its copy from OrigEdgList and
      // remove the COPIED mark from the corresponding original edge.
      // At the end of optimization, list OrigEdgList contains all original
      // edges considered during the optimization but not swapped. Such
      // edges represent the boundary of the region involved in the 
      // optimization. Thus we empty the list of copies and remove the
      // COPIED mark from original edges. This last operation can be
      // done while traversing the influence region marked as TO_DELETE
      // in function ExtDelBaseDelInflRegn.
      // The optimization may coniser an original edge outside the region 
      // of influence. It is necessary to save a copy of original 
      // topological relations before its two original adjacent triangles
      // are modified. During optimization, before an edge is swapped,
      // we compute the quadrilateral and, if the edges of the quadrilateral
      // are original and not in OrigEdgList, we insert them in OrigEdgList
      // and mark them as COPIED.
      // See also comments in function ExtSwapEdge.
      TDoubleList<TEdge> OrigEdgList;

      // VIDs of already removed vertices
      TDoubleList<int> RemovedVertexIndex;

      // Number of vertices already removed from the triangulation
      int nRemovedVertex;
      
      // Array of pointers to the constraints
      PTEdge *Constraints;
      
      // Number of constraints present in input file (.tac)
      int nConstrInFile;
      
      // Number of constraints currently present in the triangulation
      int nConstrInTRI;  
							      
      // NOTE: the input file may cointain "wrong" constraints,
      // not corresponding to an edge of the triangulation.
      // In this case nConstrInTRI is initialized to the number
      // of "correct" constraints present in the input file.

      // Used only in FindConstraintInVE to decide whether to
      // continue in case of failure
      boolean ReadingConstr; 

      // Functions required for the insertion of a constraint
      // (Add_Constraint) 
      boolean OkConstrDegree( PTVertex V );
      virtual boolean ReCheckVertex( PTVertex );
      int IsEdgeAdmissible( PTVertex V1, PTVertex V2 );
      void Edge2Constraint( PTEdge E );
      void AdjustBordersVertices_VE( TDoubleList<PTEdge> Border );
      virtual void Delete_Constr_InfluenceRegion( PTEdge newConstr );
      void TTra_Add_Constraint( PTVertex v0, PTVertex v1 );
      virtual void UpdateStep();

      // Functions required for the extension of optimization
      // beyond the region of influence of the vertex, necessaty
      // to manage the removal of constraints in case of vertices
      // with constraints incident in them.
      void ExtRIRIPrepare(TDoubleList<PTEdge> & SwapEdgeQueue, PTEdge constr1, PTEdge constr2,
                          boolean & closed, TDoubleList<TEdge> &OrigEdgList);
      void ExtRIRIPrepare(TDoubleList<PTEdge> & SwapEdgeQueue, PTEdge constraint,
                          TDoubleList<TEdge> &OrigEdgList );
      void Constraint2Edge( PTEdge E );

      void ExtRIRInitialTrg( TDoubleList<PTEdge> & InflRegnBorder,
                             TDoubleList<PTEdge> & SwapEdgeQueue   );

      void ExtRIRDelOptTrg(TDoubleList<PTEdge> & SwapEdgeQueue, TDoubleList<TEdge> &OrigEdgList);
      PTEdge ExtSwapEdge( PTEdge &E, TDoubleList<TEdge> &OrigEdgList );

      void ExtDeleteInfluenceRegion( TDoubleList<PTEdge> & InflRegnBorder );
      void RecheckOptimizedRegion();
      void ExtMT_AddComponent();
      void Del_InflRegnBorder( TDoubleList<PTEdge> & InflRegnBorder );
      void ExtDelBaseDelInflRegn();
      void ReconstructET( PTEdge E, TDoubleList<TEdge> &OrigEdgList );


      void CheckEdgeInflBorder( PTEdge E, TDoubleList<PTEdge> &border);
      void CalcEdgeWideInflRegn( PTEdge E,  TDoubleList<PTTriangle> &InternalTrgs,
                    TDoubleList<PTEdge> &Right_border, TDoubleList<PTEdge> &Left_border );
      int InitCalcEdgeWideInflRegn( PTEdge E,  TDoubleList<PTTriangle> &InternalTrgs,
                    TDoubleList<PTEdge> &Right_border, TDoubleList<PTEdge> &Left_border );

 
      int  Tra_intersect( PTEdge e1, PTEdge e2 );
      void Tra_intersect_test();

      void  Tra_Add_Constraint( PTVertex V1, PTVertex V2 );
      boolean CalcEdgInflRegn( PTEdge newConstr, TDoubleList<PTEdge> &Right_Border,
                    TDoubleList<PTEdge> &Left_Border, PTEdge & origEdg  );

      boolean IsNewConstrOK( PTVertex );
      boolean OldIsNewConstrOK( PTVertex );

      virtual void RemoveVertex();


      virtual void DetachEdge( PTEdge );
 //   virtual boolean IsVtxElim( PTVertex );  replaced by OkConstrDegree

      virtual void ReadData( const char *infname );

      void ReadConstraints( ifstream & );

       
      // Since the filtering function used more often is EqualsEV,
      // this is used as a default argument for the 
      // function pointer parameter filterFun.
      int FindEdgesInVE( PTVertex V, PTEdge filter, TDoubleList<PTEdge> & filtered, PF1b2PTEdge filterFun /*= EqualsEV*/);
      int FindEdgesInVE( PTVertex V, PTEdge filter, TDoubleList<PTEdge> & filtered);

      int FindConstraintsInVE( PTVertex V, TDoubleList<PTEdge> & IncConstr);

      static boolean (EqualsEV)( PTEdge PTE1, PTEdge PTE2 );
      static boolean (EqualsMarkConstrained)( PTEdge PTE1, PTEdge PTE2 )
      { return PTE1->Marked(CONSTRAINED) && PTE2->Marked(CONSTRAINED); }

      PTEdge FindEdgeInVE( PTVertex V, PTEdge Edge2find );
      PTEdge FindConstraintInVE( PTVertex V, PTEdge constraint );

   #ifdef array_Constraints
       virtual void WriteData( const char *outfname  ) 
       {  // WriteData1: rename in WriteData to use this one
          #ifdef DEBUG3
             DEBUG3 << "\nTDecCDT::WiteData1(" << outfname << ")" << endl;
          #endif
			 TDecimDelaunay::WriteData(outfname);
          cerr << "output " << nPts << " vertices, " << nPts-nRemovedVertex << " vertices still present in triangulation.\n";
			 WriteConstraints( outfname );
		 };
   #else
       virtual void WriteData( const char *outfname  );  
       // WriteData2: rename in WriteData to use this one
   #endif

      virtual void ConvertData(int *vNum, int *tNum, int *eNum,
                               float **vData, int **tData, int **eData);

      void WriteConstraints( const char *outfname  );


   public:
  
       TDecCDT( int iK, PMTTracer iMT );
       TDecCDT( int iK, PMTTracer iMT, boolean EXTActive, boolean ALLOWFeaturesDel, boolean ALLOWChainBrk );

};


#endif // _DECCDT_H
