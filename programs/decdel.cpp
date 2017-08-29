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
//   file   : decdel.cpp
//   author : Christian Melchiorre
//
//   Implementation of Class TDecimDelaunay, sub-class of TDestroyDelaunay,
//   for the decimation of an input triangulation.
//   As TRefineDelaunay sub-class of TBuildDelaunay, this class adds the
//   management of point lists to the base class, as well as 
//   the construction of the initial triangulation.
//


#include "defs.h"
#include "error.h"
#include "geom.h"
#include "destrdel.h"
#include "decdel.h"
#include "mttracer.h"

// -----------------------------------------------------------------------------
//  
//  Constructor of class TDecimDelaunay
//


TDecimDelaunay::TDecimDelaunay( int iK, PMTTracer iMT ) : TDestroyDelaunay( iK), DetachedPoints()
{
    #ifdef DEBUG
       DEBUG << "TDecimDelaunay Constructor" << endl;
    #endif
    
    InitialPhase = TRUE;

    MT = iMT;
}


// -----------------------------------------------------------------------------
//  
//  void TDecimDelaunay::InitialTriangulation()
//
//  This function is re-implemented to call MT_Initial() after the
//  construction of the initial triangulation.
//

void TDecimDelaunay::InitialTriangulation()
{
   
   #ifdef DEBUG
      DEBUG << "TDecimDelaunay::InitialTriangulation()" << endl;
   #endif
   
   TDestroyDelaunay::InitialTriangulation();
   
   MT_Initial();
   
   InitialPhase = FALSE;
   
}


// -----------------------------------------------------------------------------
//  
//  void TDecimDelaunay::RepositionPoint( PTPoint )
//
//  Find the triangle/edge containing the given point, and insert the point
//  in the  PointList of such triangle/edge.
//  The code is nearly equal to the one of the correspinding function
//  TRefineDelaunay::RepositionPoint()
//

void TDecimDelaunay::RepositionPoint( PTPoint PointToPos )
{

    #ifdef DEBUG
       DEBUG << "TDecimDelaunay::RepositionPoint()" << endl;
    #endif
    
    #ifdef ROBUST 
       check( (PointToPos == NULL), "TDecimDelaunay::RepositionPoint(), NULL point" );
    #endif
    
    PointLocation( PointToPos );
    
    switch( PLLocation )
    {
       case PL_EDGE:
            {
               PLEdge->AddPoint( PointToPos );
            }
	    break;
	    
       case PL_TRIANGLE:
            {
	       PLTriangle->AddPoint( PointToPos );
	    }
	    break;
	    
       case PL_VERTEX:
            {
	       error( "TDecimDelaunay::RepositionPoint(), duplicate point");
	    }
	    
       case PL_EXTERNAL:
            {
	       error( "TDecimDelaunay::RepositionPoint(), point outside convex hull" );
	    }
	    
    } // end ...switch( PLLocation )

}


// -----------------------------------------------------------------------------
//  
//  void TDecimDelaunay::DeleteInfluenceRegion()
//
//  Re-implemented in order to pass information about deleted / added
//  triangles to the MT Tracer.
//

void TDecimDelaunay::DeleteInfluenceRegion()
{

   #ifdef DEBUG
     DEBUG << "TDecimDelaunay::DeleteInfluenceRegion()" << endl;
   #endif
    
   //
   // delete old triangles
   //
   
   if ( ! InitialPhase ) MT_KillInterference();
    
   //
   // call old version of DeleteInfluenceRegion() to delete the triangles
   //
   
   TDestroyDelaunay::DeleteInfluenceRegion();
   
   //
   // Reposition "Detached" points from deleted triangles to new triangles,
   // in such a way that we are able to compute their error correctly in 
   // MT_AddComponent().
   //
   
   RepositionPoint( VertexToRemove );
   
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


// -----------------------------------------------------------------------------
//  
//  void TDecimDelaunay::AddTriangle()
//  void TDecimDelaunay::DetachTriangle()
//  void TDecimDelaunay::DetachEdge()
//
//  The following functions are re-implemented for the management of
//  PointList's.
//


void TDecimDelaunay::AddTriangle( PTTriangle T )
{
   TDestroyDelaunay::AddTriangle( T );
   T->Mark( NEW_TRIANGLE );   
}


void TDecimDelaunay::DetachTriangle( PTTriangle T )
{
   TDestroyDelaunay::DetachTriangle( T );
   
   while ( ! T->PointList.IsEmpty() )
      DetachedPoints.AddHead( T->PointList.RemoveHead() );
   
}

void TDecimDelaunay::DetachEdge( PTEdge E )
{
 
   TDestroyDelaunay::DetachEdge( E );
   
   while ( ! E->PointList.IsEmpty() )
      DetachedPoints.AddHead( E->PointList.RemoveHead() );
   
}


// ---------------------------------------------------------------------------------------------\
//
//   void TDecimDelaunay::EndTriangulation()
//
//   Re-implemented only for the call to MT_Terminate()
//

void TDecimDelaunay::EndTriangulation()
{
   TDestroyDelaunay::EndTriangulation();

    //
    // create array to contain the triangles
    //
    
    PTTriangle *TrgArray;
    int it, i;

    TrgArray = new PTTriangle[ nTrg ];
        
    check( ( TrgArray == NULL ),
           "TDecimDelaunay::EndTriangulation(), insufficient memory" );
	  
    for( it=0; it<nTrg; it++ ) TrgArray[it] = NULL;

    TDoubleList<PTTriangle> Triangles;
   
    check( (FirstTriangle==NULL), "TDecimDelaunay::EndTriangulation(), No triangles?");
   
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

       check( (CT == NULL), "TDecimDelaunay::EndTriangulation(), <2> inconsistency detected");

       MT->KillTriangle( CT );
    }    
    MT->MeshOk();

   MT_Terminate();
}


// ---------------------------------------------------------------------------------------------\
//
//   Functions to call the MT Tracer functions
//
//     void TDecimDelaunay::MT_Initial()
//     void TDecimDelaunay::MT_KillInterference()
//     void TDecimDelaunay::MT_AddComponent()
//     void TDecimDelaunay::MT_Terminate()
// 

void TDecimDelaunay::MT_Initial()
{
    
    MT->StartHistory( MT_COARSENING );
        
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


//
// INPUT: FirstTrgToDel points to an old triangle of the influence region,
// the old triangles are all marked as TO_DELETE.
//
// ALGORITHM: traverse the old triangles of the region of influence through 
// adjacency navigation, and call MT->KillTriangle(T) for each of such
// triangles; mark TO_DELETE is not removed since it is needed later in 
// TBuildDelaunay::DeleteInfluenceRegion().
//

void TDecimDelaunay::MT_KillInterference()
{

   #ifdef MT_DEBUG
     MT_DEBUG << endl << "MT_KillInterference: eliminate V" << VertexToRemove->VID 
          << ": " << VertexToRemove->x << ", " << VertexToRemove->y << ", " << VertexToRemove->z
          << " => Potential Error : " << VertexToRemove->Error << endl;
   #endif


   if ( FirstTrgToDel != NULL )
   {
   
      #ifdef ROBUST
         check( (! FirstTrgToDel->Marked( TO_DELETE ) ),
            "TDecimDelaunay::MT_KillInterference(), inconsistency detected" );
      #endif

      TDoubleList<PTTriangle> Triangles;
      
      FirstTrgToDel->Mark( MT_DELETED );
      Triangles.AddHead( FirstTrgToDel );
      
      while( !Triangles.IsEmpty() )
      {
          PTTriangle T = Triangles.RemoveHead();
          PTTriangle TT[3];
      
          // do not remove mark TO_DELETE, it is needed in 
          //    TBuildDelaunay::DeleteInfluenceRegion()
      
          T->GetTT( TT[0], TT[1], TT[2] );
      
          MT->KillTriangle( T );
      
          for( int t=0; t<3; t++ )
             if ( TT[t] != NULL && TT[t]->Marked( TO_DELETE ) 
                  && ! TT[t]->Marked( MT_DELETED ) )
                  {  
                     TT[t]->Mark( MT_DELETED );
	             Triangles.AddTail( TT[t] );
                  }
		     
      }
      
   }
   
}

//
// INPUT: FirstTriangle points to a new triangle of the just re-triangulated
// influence region; the new triangles of the influence region are all
// marked as NEW_TRIANGLE.
//
// ALGORITHM: traverse the new triangles of the region of influence through
// adjacency navigation, and call MT->MakeTriangle(T) for each of such
// triangles; at the end of traversal, call MT->MeshOk().
//

void TDecimDelaunay::MT_AddComponent()
{

   #ifdef MT_DEBUG
      if( VertexToRemove == NULL )
         MT_DEBUG << "TDecimDelaunay::MT_AddComponent() to record initial triangulation\n";
      else
         MT_DEBUG << "TDecimDelaunay::MT_AddComponent() to record re-triangulation after deleting V" << VertexToRemove->VID << endl;
   #endif


   #ifdef ROBUST
      check( (FirstTriangle == NULL || !FirstTriangle->Marked( NEW_TRIANGLE )),
             "TDecimDelaunay::MT_AddComponent(), inconsistency detected" );
   #endif
     
   TDoubleList<PTTriangle> Triangles;
     
   FirstTriangle->UnMark( NEW_TRIANGLE );
   Triangles.AddHead( FirstTriangle );
   
   while ( ! Triangles.IsEmpty() )
   {
      PTTriangle T = Triangles.RemoveHead();
      PTTriangle TT[3];

      T->GetTT( TT[0], TT[1], TT[2] );

      double Error = 0.0;  
      Error = Error;  //  to avoid a warning
      if ( ! T->PointList.IsEmpty() ) Error = T->PointList.GetHead()->Error;

      MT->MakeTriangle( T );
      
      for( int t=0; t<3; t++ )
         if ( TT[t] != NULL && TT[t]->Marked( NEW_TRIANGLE ) )
         {
            TT[t]->UnMark( NEW_TRIANGLE );
	    Triangles.AddTail( TT[t] );
         }
      
   }
      
   MT->MeshOk();
  
}


void TDecimDelaunay::MT_Terminate()
{
   MT->EndHistory();
}


// --------------------------------------------------------------------------------
//  
//  void TDecimDelaunay::WriteData( const char * )
//
//  Save triangulation to file.
//

void TDecimDelaunay::WriteData( const char *outfname )
{
    // Indeed, re-implementation is not needed.

    TTriangulation::WriteData(outfname);
}
