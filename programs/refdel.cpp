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
//   file   : refdel.cpp
//   author : Christian Melchiorre
//
//   Implementation of class TRefineDelaunay, subclass of TBuildDelaunay.
//   Unlike the base class, it builds the initial triangulation as a 
//   triangulation of the convex hull of the input points
//   (instead of taking the first three non-aligned points).
//   Then, it inserts all remaining points.
//   In summary, it re-implements method InitialTriangulation(), adds
//   method CalcConvexHull(), and some auxiliary functions.
//

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "defs.h"
#include "error.h"
#include "geom.h"

#include "ttriang.h"

#include "builddel.h"
#include "refdel.h"

#include "mttracer.h"


// -------------------------------------------------------------------------
//
//   Constructor of class TRefineDelaunay
//


TRefineDelaunay::TRefineDelaunay( PMTTracer iMT ) : TBuildDelaunay(), DetachedPoints()
{

   #ifdef DEBUG
     DEBUG << "TRefineDelaunay Constructor" << endl;
   #endif
   
   //
   // variables used in CalcConvexHull():
   //

   // number of points on convex hull
   nChPts = 0;

   // flag if we are in phase of initial triangulation
   InitialPhase = FALSE;

   // MT Tracer
   MT = iMT;
   
}


// -------------------------------------------------------------------------
//
//   void TRefineDelaunay::InitialTriangulation()
//
//   Re-implement TBuildDelaunay::InitialTriangulation() in such a way 
//   that the initial triangulation is a triangulation of the convex hull
//   rather than a triangle formed by the first three non-aligned points.
//
//   INPUT: array Points containing the pointers to the input points.
//
//   OUTPUT: a consistent status corresponding to the triangulation of the
//   input points lying on the convex hull.
//

void TRefineDelaunay::InitialTriangulation()
{

    int i;

    #ifdef DEBUG
      DEBUG << "TRefineDelaunay::InitialTriangulation()" << endl;
    #endif
  
    InitialPhase = TRUE;

    //
    // compute convex hull of input points
    //

    CalcConvexHull();

    #ifdef DEBUG
       DEBUG << "computed ch" << endl;
    #endif
    
    
    //
    // After CalcConvexHull(), the points of the convex hull are in 
    // Points[0..nChPts-1]. We insertd them in the initial triangulation.
    // Method InitialTriangulation() of the base class TBuildDelaunay,
    // builds an initial triangulation with the first three non-aligned
    // points in Points. Because of the disposition of poimnts in such 
    // array (first the convex hull points, then the others), the first
    // three non-aligned points certainly belong to the convex hull.
    //

    TBuildDelaunay::InitialTriangulation();

    #ifdef DEBUG
      DEBUG << "After TBuildDelaunay::InitialTriangulation()" << endl;
    #endif

    for ( i=3; i<nChPts; i++ )
    { 
       VertexToIns = new TVertex( Points[i] );
       check( (VertexToIns == NULL), "TRefineDelaunay::InitialTriangulation(), insufficient memory");
       
       delete (Points[i]);
       Points[i] = VertexToIns;
       
       InsertVertex();

       if ( i % 1000 == 0 ) cerr << "inserted " << i << " points in the convex hull   \r" << flush;
    }

    cout << "inserted " << nChPts << " points in the convex hull" << endl;
    

    //
    // The remaining points after the triangulation of the convex hull,
    // i.e., Points[nChPts]...Points[nPts-1], are inserted in the pointlists
    // of edges/triangles in which they fall.
    //
    
    for ( i=nChPts; i<nPts; i++ )
       RepositionPoint( Points[i] );

    //
    // the initial triangulation has been built, the next point to
    // be inserted is that of index nChPts in Points
    //

    iNextPoint = nChPts;

    cout << "triangulation of convex hull completed" << endl;
    
    MT_Initial();
    
    InitialPhase = FALSE;
 
}



// -------------------------------------------------------------------------
//
//   void TRefineDelaunay::CalcConvexHull()
//
//   Compute convex hull of input points by using Graham's algorithm.
//
//   INPUT: array Points
//
//   OUTPUT: array Points is re-ordered in such a way that, if there are
//   nChPts = N points on the convex hull, then the first N positions
//   contain such points, while the other positions contain the remaining
//   points.
//

void TRefineDelaunay::CalcConvexHull()
{
    
    int i;

    #ifdef DEBUG
      DEBUG << "TRefineDelaunay::CalcConvexHull()" << endl;
    #endif

    //
    // sort points in Points in counterclockwise order w.r.t. a point
    // internal to the set (the baricenter of the first three
    // non-aligned points).
    //
    
    GrahamSortPoints();

    //
    // search the point with maximum x coordinate (if not unique, 
    // the point of maximum y coordinate among them), that will be
    // the starting point for Graham's scan.
    //

    // index of point of max x
    int iMaxX = 0;

    // to compare x values
    double xi, yi, maxxi, maxyi;

    for( i=0; i<nPts; i++ )
    {

        xi = Points[i]->x;
        yi = Points[i]->y;

        maxxi = Points[iMaxX]->x;
        maxyi = Points[iMaxX]->y;


        //
        // index of point of max x (if not unique, of max y)
        //

        if ( (xi > maxxi) || ( (xi == maxxi) && (yi > maxyi) ) )
           iMaxX = i;

    }

    //
    // Use two temporary arrays, ChNext and ChPrev, to maintain the links
    // to a circular double list of points on the current hull.
    // The initial hull contains all points in the counterclockwise
    // order computed at the beginning.
    //

    int *ChNext = new int[nPts],
        *ChPrev = new int[nPts];
	
    check( (ChNext == NULL || ChPrev == NULL),
    	   "TRefineDelaunay::CalcConvexHull(), insufficient memory");	   


    for ( i=0; i<nPts; i++ )
    {
       ChPrev[i] = i-1;
       ChNext[i] = i+1;
    }

    // close circle

    ChPrev[ 0 ] = nPts-1;
    ChNext[ nPts-1 ] = 0;


    nChPts = nPts;
    // number of points on the hull


    //
    // Main loop of Graham's scan
    //

    int ic = iMaxX, icNext, icNext2;

    boolean complete = FALSE;
         
    do
    {


       icNext  = ChNext[ic];
       icNext2 = ChNext[icNext];

       //
       // the following test rejects vertices that are on the hull but 
       // are endpoints of two aligned edges
       //
        
       if ( Geom::Turnxy( Points[ic], Points[icNext], Points[icNext2] ) != TURN_LEFT )
       {
    
          //
          // Points[icNext] does not belong to the hull, remove it from
          // circular list ChPoints
          //

          ChNext[ic] = icNext2;
          ChPrev[icNext2] = ic;

          ChPrev[icNext] = ChNext[icNext] = -1;
	  // -1 means that the point is no longer on the hull

	  nChPts--;
          // one less point on the hull

          //
          // resume scan from point preceding ic. If ic = iMaxX
	  // since we know that Points[iMaxX] is certainly on the hull,
	  // it is not needed to step backward
          //

          if (ic != iMaxX) ic = ChPrev[ic];

       }

       else // ok, go to next step
       {
          ic = ChNext[ic];
	  if ( ic == iMaxX ) complete = TRUE;
       }


    } while( !complete );


    //
    // Now the points of the convex hull are all in the circular list
    // represented by the links in ChNext[], starting from index iMaxX. 
    // The remaining points are those with ChNext[i] = -1. 
    // Now we sort Points, in such a way that, if the points on the 
    // convex hull are N, then the pointers to them are in the first
    // N positions of Points, and the pointers to the remaining points are
    // in the following positions. We do all this in O( nPts ) time.
    // We use a new auxiliary array, PointsAux, that will replace Points.
    //

    delete( ChPrev ); ChPrev = NULL; // no longer needed


#ifndef PAOLO

    register PTPoint PointsAux;
    int iCh = 0;    // positions for the hull 0...nChPts-1
    int iNoCh = 0;  // for the other points nChPts...nPts-1

    for(;;iNoCh++)
    {
      for( ; ChNext[iNoCh]!=-1; iNoCh++ );  // find 1st point not on CH
      if ( iNoCh >= nChPts )  break;        // all points ok, exit loop
      for( iCh = iNoCh+1; ChNext[iCh]==-1 && iCh<nPts; iCh++ ); 
                                            // find first point on CH
      if ( iCh >= nPts )  break;            // no more points, exit loop
      PointsAux     = Points[iCh];
      Points[iCh]   = Points[iNoCh];
      Points[iNoCh] = PointsAux;
      ChNext[iCh]   = -1;
    }
#else
    PTPoint *PointsAux = new PTPoint[nPts];

    check( (PointsAux == NULL),
	"TRefineDelaunay::CalcConvexHull(), insufficient memory" );


    int iCh    = 0;     // posizions for convex hull 0...nChPts-1
    int iNoCh = nChPts; // for other points nChPts...nPts-1

    for( i=0; i<nPts; i++ )
    {
      if ( ChNext[i] != -1 ) // if point is on convex hull
	 PointsAux[iCh++] = Points[i];
      else
	 PointsAux[iNoCh++] = Points[i];
    }

    //
    // now replace Points with PointsAux
    //

    delete( Points );

    Points = PointsAux;
    PointsAux = NULL;
#endif

    delete( ChNext );

    #ifdef DEBUG
      DEBUG << "Points( end of CalcConvexHull ): " << endl;
      int kk;
      for( kk=0; kk<nPts; kk++ )
      {
        if ( Points[kk] == NULL ) DEBUG << "NULL"; else DEBUG << "V " << (*Points[kk]);
	if ( kk < nChPts ) DEBUG << "*";
	DEBUG << endl;
      }
      DEBUG << endl;
    #endif // DEBUG


}


// ---------------------------------------------------------------------------------
//
//   void TRefineDelaunay::GrahamSortPoints()
//
//   Sort array Points counterclockwise w.r.t. a point internal to the set.
//   Use heapsort, with an optimal computational cost of 
//   O( n log n ) and no need for auxiliary structures.
//
//   INPUT: array Points, containing the input points sorted in the order in
//   which they have been read.
//
//   OUTPUT: array Points sorted counterclockwise w.r.t. a point internal
//   to the convex hull (the baricenter of the first three non-aligned
//   points).
//

void TRefineDelaunay::GrahamSortPoints()
{
    
    int i, insI;

    #ifdef DEBUG
      DEBUG << "TRefineDelaunay::GrahamSortPoints()" << endl;
    #endif

    //
    // find baricenter of three non-aligned points in array Points
    //

    PTPoint CenterPoint = NULL;

    for( i=2; i<nPts && CenterPoint==NULL; i++ )
    {

       if ( ! Geom::Alignedxy( Points[0], Points[1], Points[i] ) )
       {

          //
          // found three non-aligned points, compute baricenter
          //

          CenterPoint = new TPoint(
                         ( Points[0]->x + Points[1]->x + Points[i]->x ) / 3.0,
                         ( Points[0]->y + Points[1]->y + Points[i]->y ) / 3.0,
                         0.0
                      );

          check( (CenterPoint == NULL),
             "TRefineDelaunay::CalcConvexHull(), insufficient memory" );

          //
          // since now CenterPoint != NULL, the loop ends
          //

       }

    } // end ...for

   
    //
    // if now CenterPoint is still NULL, then three non-aligned points
    // do not exist in the input set
    //

    check ( (CenterPoint == NULL),
            "TRefineDelaunay::CalcConvexHull(), all points are aligned");     
    
    //
    // start of real sorting algorithm
    //
        
    //
    // Elements in positions from [nPts/2] to [nPts-1] alread form a heap.
    // Insert in such heap the elements of previous positions.
    //

    for( insI=((int)nPts/2)-1; insI>=0; insI-- )
         HeapSortSift( insI, nPts-1, CenterPoint );

    //
    // now transform heap into sorted array
    //

    PTPoint pTmp;
    
    for( insI=nPts-1; insI>=1; insI-- )
    {
         // ...swap Points[0] and Points[insI]
	 pTmp = Points[0];
	 Points[0] = Points[insI];
	 Points[insI] = pTmp;
	 
	 HeapSortSift( 0, insI-1, CenterPoint );
    }
    
    // 
    // array has been sorted
    //
        
    delete( CenterPoint );
    CenterPoint = NULL;

 }


// -----------------------------------------------------------------------------------------
//
//  int  TRefineDelaunay::HeapSortSift( int first, int last, PTPoint CenterPoint )
//
//  Auxiliary function for heapsort. Merge the element of positipon
//  first into the heap.
//

void TRefineDelaunay::HeapSortSift( int first, int last, PTPoint CenterPoint )
{

   int i = first;
   int j = (2 * i) + 1; // left son of i
   
   PTPoint pTmp = NULL;
   
   while( j <= last )
   {   
   
      // find the larger of the two children of Points[i], and
      // go on into that branch
      if ( j+1 <= last ) 
         if ( GrahamSortFunc( Points[j], Points[j+1], CenterPoint ) < 0 )
	    j = j+1;
      
      if ( GrahamSortFunc( Points[i], Points[j], CenterPoint ) < 0 )
      {
          // ...swap Points[i] and Points[j]
	  pTmp = Points[i];
	  Points[i] = Points[j];
	  Points[j] = pTmp;	  

          i = j;
          j = (2 * i) + 1;	  
      }
      else break; 
   
   } // end ...while( j < nPts )

}


// -----------------------------------------------------------------------------------------
//
//  int  TRefineDelaunay::GrahamSortFunc( PTPoint, PTPoint, PTPoint )
//
//  Set the ordering of the points based on the angle formed by the segment
//  joining each point to CenterPoint (a point inside the convex hull).
//  Two points are compared in the following way:
//
//   a) a point of coordinates equal to CenterPoint is minimum 
//
//   b) if the quadrants of the two points are different, then the
//      points are in the same relation as their quadrants:
//
//         I quadrant < II quadrant < III quadrant < IV quadrant
//
//   c) if the quadrants of the two points are equal, then
//      P1 is less than P2 iff (CenterPoint, P1, P2) form a right turn.
//

int TRefineDelaunay::GrahamSortFunc( PTPoint P1, PTPoint P2, PTPoint CenterPoint )
{
     
    //
    // test if one point is the same as CenterPoint
    //
    
    boolean C1 = P1->Equalsxy( CenterPoint );
    boolean C2 = P2->Equalsxy( CenterPoint );

    if ( C1 && C2 )
       return 0;              // P1 = P2
    else if ( C1 && ! C2 )
       return -1;              // P1 < P2
    else if ( !C1 && C2 )
       return 1;             // P1 > P2

    //
    // find the quadrants of the two points (by convention, the positive 
    // x axis belongs to the I quadrant, the positive y axis belongs to 
    // the II quadrant, the negative x axis belongs to the III quadrant,
    // and the negative y axis belongs to the IV quadrant)
    //

    int QP1, QP2;
    
    // first P1
    if ( P1->x > CenterPoint->x && P1->y >= CenterPoint->y )
       QP1 = 1;
    else if ( P1->x <= CenterPoint->x && P1->y > CenterPoint->y )
       QP1 = 2;
    else if ( P1->x < CenterPoint->x && P1->y <= CenterPoint->y )
       QP1 = 3;
    else
       QP1 = 4;

    // then P2
    if ( P2->x >  CenterPoint->x && P2->y >= CenterPoint->y )
       QP2 = 1;
    else if ( P2->x <= CenterPoint->x && P2->y > CenterPoint->y )
       QP2 = 2;
    else if ( P2->x < CenterPoint->x && P2->y <= CenterPoint->y )
       QP2 = 3;
    else
       QP2 = 4;

    // compare the two quadrants

    if ( QP1 < QP2 )
       return -1;
    else if ( QP1 > QP2 )
       return 1;
    else // QP1 == QP2
    {
       switch( Geom::Turnxy( CenterPoint, P1, P2 ) )
       {
          case TURN_LEFT: return -1; 
          case TURN_RIGHT: return 1;
          case ALIGNED: return 0; 
       }
    }

    // the following code is not reached, but needed in order to
    // avoid warning from the compiler

    error( "TRefineDelaunay::GrahamSort(), this should not happen" );
    return(0);

}


// ---------------------------------------------------------------------------------
//
//   void TRefineDelaunay::NextPoint()
//
//   Put into variable VertexToIns the pointer to the next point to be
//   inserted in the triangulation.
//
//   INPUT: iNextPoint, index in array Points of the next
//   point to be inserted
//
//   OUTPUT: VertexToIns, the value of iNextPoint is incremented by one
//
//   The code is nearly equal to the one of TBuildDelaunay::NextPoint(),
//   the only difference is that PointToIns is not deleted, because it
//   could belong to the  PointList of some triangle/edge of the 
//   triangulation. It will be deleted during the 
//   DetachEdge/DetachTriangle of the entity containing it in its PointList.
//

void TRefineDelaunay::NextPoint()
{
   PTPoint PointToIns = Points[iNextPoint];
   
   #ifdef ROBUST
     check( (PointToIns == NULL), "TBuildDelaunay::NextPoint(), inconsistency detected" );
   #endif

   Points[iNextPoint] = VertexToIns = new TVertex( PointToIns );
   check( (VertexToIns == NULL), "TBuildDelaunay::NextPoint(), insufficient memory");

   //
   // Maintaining the vertices in array Points, in order of their creation,
   // is useful for function WriteData(), when we write them in order
   // of creation (in such a way that references from a triangle to the
   // VIDs of its vertices are coincident with the order in which
   // vertices are written).
   //

   iNextPoint++;
}


// -----------------------------------------------------------------------------
//
//  void TRefineDelaunay::RepositionPoint( PTPoint )
//
//  Find the triangle/edge containing the given point, and insert the
//  point into the  PointList of such triangle/edge.
//

void TRefineDelaunay::RepositionPoint( PTPoint PointToPos )
{

    #ifdef DEBUG
       cout << "TRefineDelaunay::RepositionPoint(" << (*PointToPos) << ")" << endl;
    #endif
           
    #ifdef ROBUST
       check( (PointToPos == NULL), "TRefineDelaunay::RepositionPoint(), NULL point" );
    #endif
    
    PointLocation( PointToPos );
       
    switch( PLLocation )
    {
       case PL_EDGE:
            {
	       #ifdef DEBUG
	          DEBUG << "   Insert point ( " << (*PointToPos) << " ) in E" << PLEdge->EID << endl;
	       #endif

               PLEdge->AddPoint( PointToPos ); 
	    }
	    break;
	  
       case PL_TRIANGLE:
            {
	       #ifdef DEBUG
	          DEBUG << "   Insert point ( " << (*PointToPos) << " ) in T" << PLTriangle->TID << endl;
	       #endif

  	       PLTriangle->AddPoint( PointToPos ); 
	    }
	    break;

       case PL_VERTEX:
            {
                cerr << "TRefineDelaunay::RepositionPoint(), duplicate point in input  " << (*PointToPos);
                break;
            }
       default: // PL_EXTERNAL
            {
	        cerr << (*PointToPos);
	        error( "  : TRefineDelaunay::RepositionPoint(), point outside convex hull");
	    }
	    
       }
       
}


// -----------------------------------------------------------------------------
//  
//  void TRefineDelaunay::AddTriangle()
//  void TRefineDelaunay::DetachTriangle()
//  void TRefineDelaunay::DetachEdge()
//
//  The following functions are re-implemented for the management of
//  PointList's.
//

void TRefineDelaunay::AddTriangle( PTTriangle T )
{
    TBuildDelaunay::AddTriangle(T);
    T->Mark( NEW_TRIANGLE );
}


void TRefineDelaunay::DetachTriangle( PTTriangle T )
{
    PTPoint P = NULL;
    
    TBuildDelaunay::DetachTriangle(T);
    
    //
    // Points belonging to the PointList of T must be stored separately,
    // in such a way that they can be re-inserted into the new triangles
    // created later in this update step.
    //
    
    while ( ! T->PointList.IsEmpty() )
    {
       P = T->PointList.RemoveHead();
       if ( P->Equalsxy( VertexToIns ) )
          delete(P);
       else
          DetachedPoints.AddHead( P );
    }
    
}


void TRefineDelaunay::DetachEdge( PTEdge E )
{
    PTPoint P = NULL;
    
    TBuildDelaunay::DetachEdge( E );
    
    //
    // Points belonging to the PointList of E must be stored separately,
    // in such a way that they can be re-inserted into the new triangles
    // created later in this update step.
    //
    
    while ( ! E->PointList.IsEmpty() )
    {
       P = E->PointList.RemoveHead();
       if ( P->Equalsxy( VertexToIns ) )
          delete(P);
       else
          DetachedPoints.AddHead( P );
    }
    
}


// -----------------------------------------------------------------------------
//  
//  void TRefineDelaunay::DeleteInfluenceRegion()
//
//  Re-implemented in order to pass information about deleted / added
//  triangles to the MT Tracer.
//

void TRefineDelaunay::DeleteInfluenceRegion()
{

   #ifdef DEBUG
     DEBUG << "TRefineDelaunay::DeleteInfluenceRegion()" << endl;
   #endif
    
   //
   // delete old triangles
   //
   
   if ( ! InitialPhase ) MT_KillInterference();
    
   //
   // call old version of  DeleteInfluenceRegion() to delete the triangles
   //
   
   TBuildDelaunay::DeleteInfluenceRegion();

   //
   // Reposition "Detached" points from deleted triangles to new triangles,
   // in such a way that we are able to compute their error correctly in
   // MT_AddComponent().
   //
   
   while( ! DetachedPoints.IsEmpty() )
   {
       PTPoint p = DetachedPoints.RemoveHead();
       RepositionPoint( p );       
   }

   //
   // add new triangles (thoise marked as NEW_TRIANGLE) and unmark them
   //
   
   if ( ! InitialPhase ) MT_AddComponent();

}


// ---------------------------------------------------------------------------------------------\
//
//   void TRefineDelaunay::EndTriangulation()
//
//   Re-implemented only for the call to MT_Terminate()
//

void TRefineDelaunay::EndTriangulation()
{
   TBuildDelaunay::EndTriangulation();

    //
    // create array to contain the triangles
    //
    
    PTTriangle *TrgArray;
    int it, i;

    TrgArray = new PTTriangle[ nTrg ];
        
    check( ( TrgArray == NULL ),
           "TRefineDelaunay::EndTriangulation(), insufficient memory" );
	  
    for( it=0; it<nTrg; it++ ) TrgArray[it] = NULL;

    TDoubleList<PTTriangle> Triangles;
   
    check( (FirstTriangle==NULL), "TRefineDelaunay::EndTriangulation(), No triangles?");
   
    Triangles.AddHead( FirstTriangle );
    FirstTriangle->Mark( VISITED );
   
    PTVertex v[3];
   
    PTTriangle NextTrg[3];
    PTTriangle CurTrg;
   
    it = 0;

    while( !Triangles.IsEmpty() )
    {
       
       CurTrg = Triangles.RemoveHead();
       TrgArray[it++] = CurTrg;

       // adjacent triangles
       
       CurTrg->GetTT( NextTrg[0], NextTrg[1], NextTrg[2] );

       for( i=0; i<3; i++ )
       {
	 //
         // enqueue adjacent triangles to CurTrg
	 //
	    
         if ( NextTrg[i] != NULL )
         {
            if (!( NextTrg[i]->Marked( VISITED ) ))
	    {
                 Triangles.AddTail( NextTrg[i] );
		 NextTrg[i]->Mark(VISITED);
	    }
         } //endif
       } // endfor
            
    } // end ...while( !Triangles.IsEmpty() )
    
    for( it=0; it<nTrg; it++ )
    {
       PTTriangle CT = TrgArray[it];

       check( (CT == NULL), "TRefineDelaunay::EndTriangulation(), <2> inconsistency detected");

       MT->KillTriangle( CT );
    }    
    MT->MeshOk();

   MT_Terminate();
}


// ---------------------------------------------------------------------------------------------\
//
//   Functions to call the MT Tracer functions
//
//     void TRefineDelaunay::MT_Initial()
//     void TRefineDelaunay::MT_KillInterference()
//     void TRefineDelaunay::MT_AddComponent()
//     void TRefineDelaunay::MT_Terminate()
// 


void TRefineDelaunay::MT_Initial()
{

    MT->StartHistory( MT_REFINING );
    
    //
    // Insert all triangles of the initial triangulation in the first
    // node of the MT
    //
    
    //
    // Triangles created during construction of the initial triangulation
    // are still marked as NEW_TRIANGLE (DeleteInfluenceRegion does not
    // unmark them). We can call MT_AddComponent().
    //
    
    MT_AddComponent();
    
}


void TRefineDelaunay::MT_KillInterference()
{

   #ifdef MT_DEBUG
     MT_DEBUG << endl << "MT: insert V" << VertexToIns->VID 
          << ": " << VertexToIns->x << ", " << VertexToIns->y 
          << " => ERR : " << VertexToIns->Error << endl;
   #endif

   if ( FirstTrgToDel != NULL )
   {
   
      #ifdef ROBUST
         check( (! FirstTrgToDel->Marked( TO_DELETE ) ),
            "TRefineDelaunay::MT_KillInterference(), inconsistency detected" );
      #endif

      TDoubleList<PTTriangle> Triangles;
      
      FirstTrgToDel->Mark( MT_DELETED );
      Triangles.AddHead( FirstTrgToDel );
      
      while( !Triangles.IsEmpty() )
      {
          int t;

          PTTriangle T = Triangles.RemoveHead();
          PTTriangle TT[3];
      
          // do not remove mark TO_DELETE, it is needed in
          //    TBuildDelaunay::DeleteInfluenceRegion()
      
          T->GetTT( TT[0], TT[1], TT[2] );
      
          MT->KillTriangle( T );
      
          for( t=0; t<3; t++ )
             if ( TT[t] != NULL && TT[t]->Marked( TO_DELETE ) 
                  && ! TT[t]->Marked( MT_DELETED ) )
                  {  
                     TT[t]->Mark( MT_DELETED );
	             Triangles.AddTail( TT[t] );
                  }
		     
      }
      
   }
   
}


void TRefineDelaunay::MT_AddComponent()
{

   int t;

   #ifdef ROBUST
      check( (FirstTriangle == NULL || !FirstTriangle->Marked( NEW_TRIANGLE )),
             "TRefineDelaunay::MT_AddComponent(), inconsistency detected" );
   #endif
     
   TDoubleList<PTTriangle> Triangles;
     
   FirstTriangle->UnMark( NEW_TRIANGLE );
   Triangles.AddHead( FirstTriangle );
   
   while ( ! Triangles.IsEmpty() )
   {
      PTTriangle T = Triangles.RemoveHead();
      PTTriangle TT[3];
 /*     PTVertex   TV[3];
      
      T->GetTV( TV[0], TV[1], TV[2] );
 */  
      T->GetTT( TT[0], TT[1], TT[2] );

 /*     double Error = 0.0;
      if ( ! T->PointList.IsEmpty() ) Error = T->PointList.GetHead()->Error;
 */
      MT->MakeTriangle( T );
      
      for( t=0; t<3; t++ )
         if ( TT[t] != NULL && TT[t]->Marked( NEW_TRIANGLE ) )
         {
            TT[t]->UnMark( NEW_TRIANGLE );
	    Triangles.AddTail( TT[t] );
         }
      
   }

   MT->MeshOk();
  
}


void TRefineDelaunay::MT_Terminate()
{
   MT->EndHistory();
}


// --------------------------------------------------------------------------------
//  
//  void TRefineDelaunay::WriteData( const char * )
//
//  Save triangulation to file.
//

void TRefineDelaunay::WriteData( const char *outfname )
{
    // Indeed, re-implementation is not needed.

    TTriangulation::WriteData(outfname);
}
