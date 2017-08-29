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
//  Implementation of class TDecErrDelaunay, sub-class of TDecimDelaunay,
//  for decimation of a Delaunay triangulation with an error-based choice
//  of the next vertex to be removed.
//


//#include <time.h>

#include "defs.h"
#include "geom.h"
#include "ttriang.h"
#include "tbtree.h"
#include "tdoublelist.h"
#include "mttracer.h"
#include "decdel.h"
#include "decerrdel.h"


// ----------------------------------------------------------------------------
// 
//  Constructor of class TDecErrDelaunay
//


TDecErrDelaunay::TDecErrDelaunay( int iK, PMTTracer iMT, int iRecalcError ) 
   : TDecimDelaunay( iK, iMT )
{

   #ifdef DEBUG
      DEBUG << "TDecErrDelaunay Constructor" << endl;
   #endif 
   
   check( (iRecalcError != RECALC_APPROX && iRecalcError != RECALC_EXACT ),
      "TDecErrDelaunay constructor: invalid value for parameter RecalcError" );
      
   RecalcError = iRecalcError;
   
}


// ----------------------------------------------------------------------------
// 
//  void TDecErrDelaunay::InitialTriangulation()
//
//  Redefine function inherited from base classes  TDestroyDelaunay and
//  TDecimDelaunay. The difference is that we insert removable vertices
//  in a balanced sorted binary tree TBTree<PTVertex> ElimVtxTree, 
//  rather than in a list ElimVtx.
//
//  INPUT: array Points of pointers to the vertices of the
//  input initial triangulation
//
//  OUTPUT: Tree ElimVtxTree, containing the pointers to the removable
//  vertices of the triangulation, sorted based on their errors, as
//  (approximatedly) computed by function
//  TDecErrDelaunay::RecalcVertexError()
//


void TDecErrDelaunay::InitialTriangulation()
{
    int v;

    #ifdef DEBUG
       DEBUG << "TDecErrDelaunay::InitialTriangulation()" << endl;
    #endif
    
    for( v=0; v<nPts; v++ )
    {
       PTVertex V = (PTVertex)(Points[v]);
       
       
       if ( v % 1000 == 0 ) cerr << "Initial error computed, " << v << " vertices\r";
       
       if ( IsVtxElim( V ) && OkDegree( CalcDegree(V) ) )  // TRADEBUG: aggiunta chiamata a OkDegree(...)
       {
            RecalcVertexError( V );
            ElimVtxTree.Insert( V );
       }
    }

    cerr << "Initial error computed, " << nPts << " vertices" << endl;
    
    MT_Initial();
    
    InitialPhase = FALSE;
    
}



// ----------------------------------------------------------------------------
// 
//  void TDecErrDelaunay::ReCheckVertex( PTVertex V )
//
//  Redefined in such a way that, before re-inserting a vertex in tree
//  ElimVtxTree, we re-compute its error, in
//  TDestroyDelaunay::DeleteInfluenceRegion()
//

boolean TDecErrDelaunay::ReCheckVertex( PTVertex V )
{

   if ( TDecimDelaunay::ReCheckVertex(V) )
   {
      RecalcVertexError( V );
      return(TRUE);
   }
   else
      return(FALSE);
}


// ----------------------------------------------------------------------------
//
//  void TDecErrDelaunay::RecalcVertexError( PTVertex V )
//
//  Recompute error associated with vertice V in an approximated or
//  exact way, depending on parameter RecalcError. Call one of the two
//  functions RecalcVertexErrorApprox(V) or RecalcVertexErrorExact(V).
//

void TDecErrDelaunay::RecalcVertexError( PTVertex V )
{

   #ifdef DEBUG
      DEBUG << "TDecErrDelaunay::RecalcVertexError( V" << V->VID << " )" << endl;      
   #endif
   
   if ( RecalcError == RECALC_APPROX )
      RecalcVertexErrorApprox( V );
   else
      RecalcVertexErrorExact( V );
}


// ----------------------------------------------------------------------------
//
//  void TDecErrDelaunay::RecalcVertexErrorApprox( PTVertex V )
//
//  Recompute approximated error associated with vertex V, as the difference
//  between the height of V and the mean of the heights of adjacent vertices
//  of V (the vertices of its influence region).
//

void TDecErrDelaunay::RecalcVertexErrorApprox( PTVertex V )
{

   int i, j;

   #ifdef DEBUG
     DEBUG << "TDecErrDelaunay::RecalcVertexErrorApprox( V" << V->VID << " )" << endl;
   #endif
   
   PTTriangle TFirst;
   
   PTEdge   EFirst = V->VE[0];
   PTVertex VFirst = ( EFirst->EV[0] != V ? EFirst->EV[0] : EFirst->EV[1] );
   
   if ( EFirst->OnConvexHull() )
   {
       
       //
       // Must decide if it is right to start from EFirst to traverse
       // moving in counterclockwise order, all the region of influence
       // of V. In case starting from EFirst in counterclockwise order
       // we get out of the convex hull, then we must start from the other
       // edge adjacent to V, V->VE[1]
       //  
       
       TFirst = ( EFirst->ET[0] != NULL ? EFirst->ET[0] : EFirst->ET[1] );
  
       #ifdef ROBUST
          check( (TFirst == NULL), 
	        "TDecErrDelaunay::RecalcVertexErrorApprox(), <0> inconsistency detected");
       #endif
            
       PTVertex v[3];
       TFirst->GetTV( v[0], v[1], v[2] );
       
       // search for the vertex of TFirst different from
       
       for( i=0; i<3; i++ )
          if ( v[i] != V && v[i] != VFirst ) break;
	 
	 
       #ifdef ROBUST
          check( (i>=3), "TDecErrDelaunay::RecalcVertexErrorApprox(), <1> inconsistency detected" );
       #endif

       if ( Geom::Turnxy( V, VFirst, v[i] ) != TURN_LEFT )
          EFirst = V->VE[1];
      
   } // end ...if( EFirst->OnConvexHull() )
   

   
   // triangle after EFirst, in counterclockwise order
   
   TFirst = ( EFirst->EV[0] == V ? EFirst->ET[0] : EFirst->ET[1] );
   
   #ifdef ROBUST
      check( ( TFirst==NULL), "TDecErrDelaunay::RecalcVertexErrorApprox(), <2> inconsistency detected" );
   #endif
   
   
   double zvtx = 0.0; 
   int    nvtx = 0;
   
   PTEdge ENext = EFirst;
   PTTriangle TNext;

   do   
   {
      // the endpoint of ENext different from V belongs to the boundary
      // of the region of influence of V
      
      PTVertex VNext = ( ENext->EV[0] != V ? ENext->EV[0] : ENext->EV[1] );
      
      zvtx += VNext->z;
      nvtx++;
      
      // go to next one
      
      TNext = ( ENext->EV[0] == V ? ENext->ET[0] : ENext->ET[1] );
      
      if ( TNext != NULL )
      {
         // search for index of E in adjacency relation TE of TNext
	 
	 for( j=0; j<3; j++ )
	   if ( TNext->TE[j] == ENext ) break;
	   
	 #ifdef ROBUST
	    check( (j>=3), "TDecErrDelaunay::RecalcVertexErrorApprox(), <3> inconsistency detected" );
	 #endif

         ENext = TNext->TE[(j+2)%3];  // edge preceding E in TE of TNext
	 
      }
           
   } while ( TNext != NULL && ENext != EFirst );
   
      
   // compute error of V as difference between the height of V and the
   // mean height of vertices adjacent to V
   
   V->Error = Geom::Abs( V->z - ( zvtx / nvtx ) );
      
   #ifdef DEBUG 
      DEBUG << "=> " << V->Error << endl;
   #endif   
}


// ----------------------------------------------------------------------------
//
//  void TDecErrDelaunay::RecalcVertexErrorExact( PTVertex V )
//
//  Recompute exact error associated with vertex V, by simulating the
//  removal of V from the triangulation: we find the triangle in which
//  V falls after the retriangulation, and compute its error as the 
//  difference between the height of V and the height of its vertical 
//  projection on such triangle.
//

void TDecErrDelaunay::RecalcVertexErrorExact( PTVertex V )
{
 
   int i, j;

   #ifdef DEBUG
      DEBUG << "TDecErrDelaunay::RecalcVertexErrorExact( V" << V->VID << " )" << endl;
   #endif
   
  
   PTTriangle TFirst;
   
   PTEdge   EFirst = V->VE[0];
   PTVertex VFirst = ( EFirst->EV[0] != V ? EFirst->EV[0] : EFirst->EV[1] );
   

   if ( EFirst->OnConvexHull() )
   {
       
       TFirst = ( EFirst->ET[0] != NULL ? EFirst->ET[0] : EFirst->ET[1] );
  
       #ifdef ROBUST
          check( (TFirst == NULL), "TDecErrDelaunay::RecalcVertexErrorApprox(), <0> inconsistency detected");
       #endif
            
       PTVertex v[3];
       TFirst->GetTV( v[0], v[1], v[2] );

       // search for vertex of TFirst different from
       
       for( i=0; i<3; i++ )
          if ( v[i] != V && v[i] != VFirst ) break;
	 
       #ifdef ROBUST
          check( (i>=3), "TDecErrDelaunay::RecalcVertexErrorApprox(), <1> inconsistency detected" );
       #endif
	  
       if ( Geom::Turnxy( V, VFirst, v[i] ) != TURN_LEFT )
          EFirst = V->VE[1];
      
   } // end ...if( EFirst->OnConvexHull() )
   

   // triangle after EFirst, in counterclockwise order
   
   TFirst = ( EFirst->EV[0] == V ? EFirst->ET[0] : EFirst->ET[1] );
   
   #ifdef ROBUST
      check( ( TFirst==NULL), "TDecErrDelaunay::RecalcVertexErrorApprox(), <2> inconsistency detected" );
   #endif
   
   
   PTEdge ENext = EFirst;
   PTTriangle TNext;   
   
   //
   // PHASE 1: find region of influence
   //
   // In the following loop, the vertices of the polygon of influence of
   // V are inserted, in counterclockwise order, in the doubly linked
   // list of PTVertex, InflVtxs.
   //

   
   TDoubleList<PTVertex> InflVtxs;   

   do   
   {
      // the endpoint of ENext different from V belongs to the boundary
      // of the region of influence of V
      
      PTVertex VNext = ( ENext->EV[0] != V ? ENext->EV[0] : ENext->EV[1] );
      
      #ifdef ROBUST
        check( (VNext == NULL), "TDecErrDelaunay::RecalcVertexErrorExact() <a> inconsistency detected" );
      #endif
 
      InflVtxs.AddTail( VNext );
      
      // go to the next one
      
      TNext = ( ENext->EV[0] == V ? ENext->ET[0] : ENext->ET[1] );
      
      if ( TNext != NULL )
      {
         // search for index of E in adjacency relation TE of TNext
	 
	 for( j=0; j<3; j++ )
	   if ( TNext->TE[j] == ENext ) break;
	   
	 #ifdef ROBUST
	    check( (j>=3), "TDecErrDelaunay::RecalcVertexErrorApprox(), <3> inconsistency detected" );
	 #endif

         ENext = TNext->TE[(j+2)%3]; // edge preceding E in TE of TNext

         #ifdef ROBUST
           check( (ENext == NULL), "TDecErrDelaunay::RecalcVertexErrorExact() <a> inconsistency detected" );
         #endif
	 
      }
           
   } while ( TNext != NULL && ENext != EFirst );
   

   //
   // PHASE 2: retriangulatione of the polygon of influence
   //
   // The 'ears' algorithm is used to find the triangle containing V...
   //
   
    
    // the next three vertices (in counterclockwise order) on InflRegnAux
    PTVertex v0, v1, v2;

    // iteratore for list of vertices
    TDoubleListIterator<PTVertex> IAuxIter( &InflVtxs );
    IAuxIter.Restart();
    
    // found triangle containing V
    boolean found = FALSE;
    
    
    #ifdef ROBUST
       check( (InflVtxs.Lenght() < 3), 
          "TDecErrDelaunay::RecalcVertexErrorExact(), <5> inconsistency detected" );
    #endif
        
    
    v0 = IAuxIter.Current()->object;
    IAuxIter.GoNext(); 
    v1 = IAuxIter.Current()->object;
    IAuxIter.GoNext(); 
    v2 = IAuxIter.Current()->object;
    
    IAuxIter.Restart();
    IAuxIter.GoNext();
    
    
    //
    // v1 is the current vertex in the list (pointed by IAuxIter), 
    // v0 the previous one, v2 the next one
    //

    while( InflVtxs.Lenght() >= 3 )
    {

       #ifdef DEBUG
          TDoubleListIterator<PTVertex> I( &InflVtxs );
	  I.Restart();
	  while( ! I.EndOfList() )
	  {
	     DEBUG << " -> V" << I.Current()->object->VID;
	     I.GoNext();
	  }
	  DEBUG << endl;
       #endif

 
       if ( Geom::Turnxy( v0, v1, v2 ) == TURN_LEFT && OkTriangle( v0, v1, v2, InflVtxs ) )
       {
 

           // triangle v0, v1, v2 is "virtually" inserted in the
           // retriangolation of the region of influence of V
	   
	   // if triangle v0, v1, v2 contains V => FOUND!!!
	   
	   if ( Geom::Turnxy( v0, v1, V ) != TURN_RIGHT &&
	        Geom::Turnxy( v1, v2, V ) != TURN_RIGHT &&
		Geom::Turnxy( v2, v0, V ) != TURN_RIGHT ) 
		{
		   found = TRUE;
		   break; // exit while
		}
		else // triangle does not contain V, remove v1 from list
		{
		   #ifdef ROBUST
		     check( (InflVtxs.Lenght() == 3), 
		        "TDecErrDelaunay::RecalcVertexErrorExact() <6> inconsistency detected" );
	           #endif
		   
		   // remove v1 from list InflVtxs, set v1 to
		   // v2 and v2 to the successor of v2.

		   IAuxIter.GoNext(); 
		   if ( IAuxIter.EndOfList() )
		   {
		      InflVtxs.RemoveLast();
		      IAuxIter.Restart(); 
		   }
		   else
		   {
		      InflVtxs.RemoveBefore( IAuxIter.Current() );
		   }
		   
		   v1 = IAuxIter.Current()->object;

		   IAuxIter.GoNext(); if ( IAuxIter.EndOfList() ) IAuxIter.Restart(); 
		   v2 = IAuxIter.Current()->object;
	   
                   // ...riporto IAuxIter a puntare al nuovo v1
   		   if ( IAuxIter.StartOfList() ) IAuxIter.GoLast(); else IAuxIter.GoPrev();
	       
	       }   
	 
	}
	else // ... not OkTriangle( v0, v1, v2, InflVtxs )
	{
	   // move to next vertex in the list
	
	   #ifdef ROBUST
	       check( (InflVtxs.Lenght() == 3), 
	           "TDecErrDelaunay::RecalcVertexErrorExact() <7> inconsistency detected" );
	   #endif
	
	   v0 = v1;
	   IAuxIter.GoNext(); if ( IAuxIter.EndOfList() ) IAuxIter.Restart();
	   v1 = IAuxIter.Current()->object;
	   IAuxIter.GoNext(); if ( IAuxIter.EndOfList() ) IAuxIter.Restart();
	   v2 = IAuxIter.Current()->object;
	   
	   // IAuxIter points to v2, reset it to v1
	   if ( IAuxIter.StartOfList() )
	      IAuxIter.GoLast();
	   else
	      IAuxIter.GoPrev();
	   
	}
	
	          
    } // end ...while( InflVtxs.Lenght() >= 3 )
    
   
    #ifdef ROBUST
       check( (!found), "TDecErrDelaunay::RecalcVertexErrorExact(), <8> inconsistency detected" );
    #endif
    
    //
    // v0, v1, v2 now form the triangle containing V... compute the
    // errore of V as difference between height of V and height of its
    // vertical projection on the triangle
    //
    
    
    V->Error = Geom::Abs( V->z - Geom::Trianglez( v0, v1, v2, V ) );
    
    #ifdef DEBUG
       DEBUG << "Errore di V" << V->VID << " : " << V->Error << endl;
    #endif

}


// -------------------------------------------------------------------------
//
//  boolean TDecErrDelaunay::OkTriangle( PTVertex... )
//
//  Almost equal to TDestroyDelaunay::OkTriangle, with two differences:
//
//   - the boundary of the influence region is taken as a list of vertices,
//     not of edges (in TDestroyDelaunay we passed the list of edges 
//     bacause it was available in InflRegnAux, and it was needed to
//     reconstruct adjacency relations TE/ET among triangles)
//
//   - we check if a triangle not only  contains a vertex inside it, but
//     also if it contains a vertex in its in-circle.
//   PAOLA: this seems to me a nonsense!!!!
//
//   Note: if it works, the method to do the incircle test in OkTriangle 
//   avoids the optimization of the new triangulation of the region of
//   influence of V after removing V. Thus, it could be used also in 
//   TDestroyDelaunay::RetriangulateInfluenceRegion for the actual 
//   elimination of V from the currennt triangulation.
//


boolean TDecErrDelaunay::OkTriangle( PTVertex v0, PTVertex v1, PTVertex v2,
                                     TDoubleList<PTVertex> &InflVtxs )
{
    #ifdef DEBUG
       DEBUG << "TDecErrDelaunay::OkTriangle( V"
             << v0->VID << ", V" << v1->VID << ", V" << v2->VID
             << " )" << endl;
       
    #endif       
  
  
    //
    // circum-circle of "virtual" triangle v0, v1, v2
    //
      
    double cx, cy, cr;
    Geom::CalcCirclexy( v0, v1, v2, cx, cy, cr );
      
    //
    // check if every vertex in InflVtxs (different from v0, v1 or v2)
    // is outside circle cx, cy, cr
    //
    
    TDoubleListIterator<PTVertex> Iter( &InflVtxs );
           
    Iter.Restart();   

    while ( !Iter.EndOfList() )
    {
  
       PTVertex VCurr = Iter.Current()->object;
       
       if ( VCurr!=v0 && VCurr!=v1 && VCurr!=v2 )
       {
           // test incircle: is VCurr inside circle through v0, v1, v2 ?

           double ds = Geom::Distancexy( cx, cy, VCurr->x, VCurr->y );

	   if ( Geom::LtDouble( ds, cr ) )
   	      return( FALSE );
       }
	     
       // else ...go to next vertex
       
       Iter.GoNext();
         
    }
    
    return( TRUE );
}
