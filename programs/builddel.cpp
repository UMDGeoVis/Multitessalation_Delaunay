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
//   file:      builddel.cpp
//   author : Christian Melchiorre
//
//   Implementatione of class TBuildDelaunay, for the construcyion of a
//   Delaunay triangulation starting from a list of input points.
//


#include <iostream>
#include <iomanip>
#include <fstream>

#include "defs.h"
#include "error.h"
#include "geom.h"
#include "ttriang.h"
#include "ttrianggc.h"

#include "basedel.h"
#include "builddel.h"



// ---------------------------------------------------------------------------------
//
//   Constructor of class TBuildDealunay
//


TBuildDelaunay::TBuildDelaunay() : TDelaunayBase()
{   
    #ifdef DEBUG
     DEBUG << "TBuildDelaunay Constructor" << endl;
    #endif // DEBUG 
}



// ---------------------------------------------------------------------------------
//
//  void TBuildDelaunay::ReadData( char *fname )
//
//  Read input data from file. The input file consists of one integer n
//  specifying the number of points, floowed by a sequence of 3*n real 
//  numbers specifying the coordinates of the points.
//
//  OUTPUT: list Points, containing the input points
//

void TBuildDelaunay::ReadData( const char *infname )
{
   register int i;

   #ifdef DEBUG
    DEBUG << "\nTBuildDelaunay::ReadData()" << endl;
   #endif // DEBUG
    
   ifstream inFile;
   
   inFile.open( infname );
   
   check( !inFile, "TBuildDelaunay::ReadData(), cannot open input file");

   //
   // read number of points
   //
   
   inFile >> nPts;
   
   check( (nPts < 3), "input with less than 3 points" );
   
   //
   // read input points and put in array Points pointers to such points
   //
      
   Points = new PTPoint[nPts];
   check( (Points == NULL), "TBuildDelaunay::ReadData(), insufficient memory for Points");

   int ip = 0;
   PTPoint p = NULL;
   
   cerr << endl;

   for( i=0; i<nPts; i++ )
   {
       p = new TPoint();
      
       check( (p == NULL), "TBuildDelaunay::ReadData(), insufficient memory for TPoint");
       check( (inFile.eof()), "TBuildDelaunay::ReadData(), unexpected End Of File");
      
       //
       // read point (its three coordinates)
       //
     
       inFile >> (*p); 
       
       #ifdef DEBUG
        DEBUG << "read point: " << (*p) << endl;
       #endif // DEBUG
       
       //
       // add point to the end of list Points
       //
       
       Points[ip++] = p; 
       if ( ip % 1000 == 0 ) cerr << "\rread " << ip << " points" << flush;  
   }

   cerr << "\rread " << ip << " points" << endl;

   inFile.close();
   
}


// ---------------------------------------------------------------------------------
//
//   void TBuildDelaunay::InitialTriangulation()
//
//   Create initial triangulation containing one triangle formed by the
//   first three non-aligned input points.
//  
//   INPUT: lista Point
//
//   OUTPUT: a consistent status of TBuildDelaunay representing the 
//   triangulation of the first three non-aligned points in Points
//

void TBuildDelaunay::InitialTriangulation()
{

   #ifdef DEBUG
    DEBUG << "\nTBuildDelaunay::InitialTriangulation()" << endl;
   #endif // DEBUG

   //
   // search for the first three non-aligned points in Points[],
   // and store them in  Points[0], Points[1] e Points[2], in 
   // counterclockwise order
   //

   int pid = 2;

   while( Geom::Alignedxy( Points[0], Points[1], Points[pid] ) )
   {
      pid++;
      check( (pid >= nPts), "TBuildDelaunay::InitialTriangulation(), all points are aligned" );
   }

   // bring Points[pid] in 3rd position in the array (if not already there)
    
   if ( pid != 2 )
   {  
      PTPoint pTmp = Points[pid];
      Points[pid] = Points[2];
      Points[2] = pTmp;
   }

   //
   // make sure that Points[0], Points[1], Points[2] form a
   // counterclockwise sequence
   //

   if ( Geom::Turnxy( Points[0], Points[1], Points[2] ) != TURN_LEFT )
   {
      // swap!
      PTPoint pTmp = Points[1];
      Points[1]   = Points[2];
      Points[2]   = pTmp;
   }

   #ifdef DEBUG
     DEBUG << "Points of first triangle: " << (*Points[0]) << ", " << (*Points[1])
           << ", " << (*Points[2]) << endl;
   #endif // DEBUG
      
   //
   // create three vertices corresponding to the three points, and delete
   // the three points
   //
   
   PTVertex v0 = new TVertex( Points[0] ); 
   PTVertex v1 = new TVertex( Points[1] ); 
   PTVertex v2 = new TVertex( Points[2] ); 
   
   check( ( v0 == NULL || v1 == NULL || v2 == NULL ),
	   "TBuildDelaunay::InitialTriangulation(), insufficient memory" );
    
   pid = Points[0]->PID;  delete Points[0];
   Points[0] = v0;
   if(CDT_SPIRITO)
      OrderInitial[pid] = 0;

   pid = Points[1]->PID;  delete Points[1]; 
   Points[1] = v1;
   if(CDT_SPIRITO)
     OrderInitial[pid] = 1;

   pid = Points[2]->PID;  delete Points[2];
   Points[2] = v2;
   if(CDT_SPIRITO)
     OrderInitial[pid] = 2;
   
   iNextPoint = 3;   
   
   //
   // create three edges of triangle v0,v1,v2, in counterclockwise order
   //
   
   #ifdef _GC_ON
      PTEdge e0 = GC::NewEdge( v0, v1 );
      PTEdge e1 = GC::NewEdge( v1, v2 );
      PTEdge e2 = GC::NewEdge( v2, v0 );
   #else
      PTEdge e0 = new TEdge( v0, v1 );
      PTEdge e1 = new TEdge( v1, v2 );
      PTEdge e2 = new TEdge( v2, v0 );

      check( (e0 == NULL || e1 == NULL || e2 == NULL ),
         "TBuildDelaunay::InitialTriangulation(), insufficient memory" );
   #endif // _GC_ON
   
   //
   // now create the triangle
   //
   
   #ifdef _GC_ON
      PTTriangle t0 = GC::NewTriangle( e0, e1, e2 );
   #else
      PTTriangle t0 = new TTriangle( e0, e1, e2 );
      check( (t0 == NULL),  "TBuildDelaunay::InitialTriangulation(), insufficient memory" );
   #endif // _GC_ON   
   
   //
   // "inverse" relations, i.e., adjacency relations from one entity to
   // another entity of higher dimension, are still to be set.
   // Such relations are not initialized by the constructors.
   //
   
   //
   // partial Vertex-Edge relation
   //
   
   v0->VE[0] = e2;    v0->VE[1] = e0;
   v1->VE[0] = e0;    v1->VE[1] = e1;
   v2->VE[0] = e1;    v2->VE[1] = e2;
   
   //
   // Edge-Triangle relation. By now, the three edges are all adjacent to
   // the left (ET[0]) to t0, and lie on the convex hull ( ET[1] = NULL ).
   //

   e0->ET[0] = t0;    e0->ET[1] = NULL;   
   e1->ET[0] = t0;    e1->ET[1] = NULL;
   e2->ET[0] = t0;    e2->ET[1] = NULL;
   
   AddTriangle( t0 );

//DEBUG
//CheckIsolatedPoint(v0);
//CheckIsolatedPoint(v1);
//CheckIsolatedPoint(v2);
   
}


// ---------------------------------------------------------------------------------
//
//   void TBuildDealunay::UpdateStep()
//
//   Perform one update step on the current triangulation, by inserting a 
//   new point.
//   First find the new point to add (the first of list Points),
//   then locate such point in the current triangulation (by function
//   TTriangulation::PointLocation()), then find and re-triangulate
//   the region of influence of the point.
//

void TBuildDelaunay::UpdateStep()
{

   #ifdef DEBUG
    DEBUG << "\nTBuildDelaunay::UpdateStep()" << endl;
   #endif // DEBUG

   //
   // TBuildDelaunay::NextPoint() put into variable VertexToIns the next
   // point to be inserted as a vertex of the triangulation
   //
   
   NextPoint();
   
   //
   // TBuildDelaunay::InsertVertex() inserts vertex VertexToIns in the
   // current triangulation
   //
   
   InsertVertex();

}


// ---------------------------------------------------------------------------------
//
//   void TBuildDelaunay::NextPoint()
//
//  Put into variable VertexToIns the pointer to the next point to be
//  inserted in the triangulation.
//
//  INPUT: variable iNextPoint, index in array Points of the next point
//  to be inserted.
//
//  OUTPUT: VertexToIns. The value of iNextPoint is incremented by one.
//
//  NOTE: before calling this function, the point was stored as TPoint,
//  here it is stored as TVertex and the TPoint that previously 
//  contained it is deleted.
//

void TBuildDelaunay::NextPoint()
{
   PTPoint PointToIns = Points[iNextPoint];
   
   #ifdef ROBUST
     check( (PointToIns == NULL), "TBuildDelaunay::NextPoint(), inconsistency detected" );
   #endif

   Points[iNextPoint] = VertexToIns = new TVertex( PointToIns );
   check( (VertexToIns == NULL), "TBuildDelaunay::NextPoint(), insufficient memory");

   //
   // Maintaining new vertices sorted in order of creation in array Points
   // is useful to write them in order of insertion in the triangulation.
   //

   delete PointToIns;

   iNextPoint++;
}


// ---------------------------------------------------------------------------------
//
//   void TBuildDealunay::InsertVertex()
//
//   Insert the vertex pointed by variable VertexToIns into the current
//   triangulation.
//
//   INPUT: VertexToIns.
//
//   OUTPUT: a consistent status of TBuildDelaunay, representing the new
//   triangulation obtained by adding vertex VertexToIns.
//      
//   ALGORITHM:
//     - PointLocation ( VertexToIns );
//     - CalcInfluenceRegion ();
//     - RetriangulateInfluenceRegion ();
//

void TBuildDelaunay::InsertVertex()
{

   #ifdef DEBUG
    DEBUG << "\nTBuildDelaunay::InsertVertex()" << endl;
   #endif // DEBUG
     
   //
   // Locate point in the triangulation.
   // This call sets status variables PLLocation, PLTriangle, PLEdge
   // (the meaning of such variables is explained in the comments
   // of function TTriangulation::PointLocation)
   //

   #ifdef DEBUG
    DEBUG << "Insert point: " << (*((PTPoint)VertexToIns))  
         << " = V" << VertexToIns->VID << endl;
   #endif // DEBUG
   
   PointLocation ( VertexToIns );
   
   //
   // If the point VertexToIns coincides with an existing vertex is ignored.
   // Otherwise we compute its region of influence and we re-triangulate it.
   //
   
   if ( PLLocation != PL_VERTEX )
   {
      CalcInfluenceRegion ();
      RetriangulateInfluenceRegion ();
   }
   else 
   {
      cerr << "\n Duplicated point : " << *VertexToIns;
      cerr << "\n It is the same as : " << *PLVertex;
      cerr << "\n Not inserted." << endl;
   }
}


// ---------------------------------------------------------------------------------
//
//  void TBuildDelaunay::CalcInfluenceRegion()
//
//  Compute the region of influence of the new point to be inserted in
//  the current triangulation.
// 
//  INPUT: status variables VertexToIns (point to be inserted),
//         PLLocation, PLEdge, PLTriangle (result of point location
//         computed by PointLocation()).
//
//  OUTPUT: fill double list InflRegnBorder of PTEdges with pointers to
//         the boundary edges of the region of influence, sorted
//         counterclockwise, and the list InfluenceRegn of triangles 
//         inside the region  of influenceof point VertexToIns.
//


void TBuildDelaunay::CalcInfluenceRegion()
{

   #ifdef DEBUG
    DEBUG << "\nTBuildDelaunay::CalcInfluenceRegion()" << endl;
   #endif // DEBUG
     
   // ckeck, just for sure
   
   #ifdef ROBUST
     check( (!InflRegnBorder.IsEmpty()),
            "TBuildDelaunay::CalcInfluenceRegion(), inconsistency detected" );
   #endif 

   
   FirstTrgToDel = NULL;
  
   //
   // call a different procedure depending on whether the point 
   // is inside or outside the current triangulation
   //
      
   if ( PLLocation == PL_EXTERNAL )  
   {
       // outside
       InitInflRegnExternal();
   }
   else
   {
       // inside
       InitInflRegnInternal();
   }
   
   CalcInflRegnMain();
   
}


// ---------------------------------------------------------------------------------
//  
//  void TBuildDelaunay::InitInflRegnInternal()
//
//  Find the region of influence in case the point VertexToIns falls 
//  inside the domain of the current triangulation.
//
//  INPUT/OUTPUT: see description of TDelaunay::CalcInfluenceRegion()
//

void TBuildDelaunay::InitInflRegnInternal()
{
    int ie;

    #ifdef DEBUG
     DEBUG << "\nTBuildDelaunay::InitInflRegnInternal()" << endl;
    #endif // DEBUG
    

    PTTriangle T0=NULL, T1=NULL;

    //
    // initialize region of influence with the triangle containing the
    // point if PLLocation == PL_TRIANGLE, with the two triangles adjacent
    // to the edge containing the point if PLLocation == PL_EDGE.
    //
    // SPECIAL CASE: PLLocation == PL_EDGE with PLEdge on the convex hull!
    //


    if ( PLLocation == PL_TRIANGLE ) 
    {
    
    
        FirstTrgToDel = T0 = PLTriangle;
	
	//
	// insert edges of triangle in the double list representing the
	// boundary of the region of influence, and mark them as INFL_BORDER
	//
	
	T0->Mark( TO_DELETE );
	
	InflRegnBorder.AddTail( T0->TE[0] );
	T0->TE[0]->Mark( INFL_BORDER );
	
	InflRegnBorder.AddTail( T0->TE[1] );
	T0->TE[1]->Mark( INFL_BORDER );
        
	InflRegnBorder.AddTail( T0->TE[2] );
	T0->TE[2]->Mark( INFL_BORDER );
        	
    }
    else if ( PLLocation == PL_EDGE )
    {
    
        PLEdge->Mark( TO_DELETE );
	
        T0 = PLEdge->ET[0];
	T1 = PLEdge->ET[1];
	
	FirstTrgToDel = ( T0 != NULL ? T0 : T1 );
	
        #ifdef ROBUST
	   check( (FirstTrgToDel == NULL),
	          "TBuildDelaunay::InitInflRegnInternal(), <1> this should not happen!" );
	#endif
	
	
	//
	// put in list the edges of T0 different from PLEdge. Check
	// that T0 (as later T1) is different from NULL
	//
	
	if ( T0 != NULL )
	{
	    T0->Mark( TO_DELETE );
	    
	    for ( ie=0; ie<3; ie++ )
	        if ( T0->TE[ie] == PLEdge )
		   break;
	    
	    #ifdef ROBUST
               check( (ie>=3), "TBuildDelaunay::InitInflRegnInternal(), <2> this should not happen!" );
	    #endif   
	    	    
	    InflRegnBorder.AddTail( T0->TE[(ie+1)%3] );
	    T0->TE[(ie+1)%3]->Mark( INFL_BORDER );
	    
	    InflRegnBorder.AddTail( T0->TE[(ie+2)%3] );
	    T0->TE[(ie+2)%3]->Mark( INFL_BORDER );

	}
	    
	
	//
	// same thing for T1
	//    
	
	if ( T1 != NULL )
	{
	    T1->Mark( TO_DELETE );
	    
	    for ( ie=0; ie<3; ie++ )
	        if ( T1->TE[ie] == PLEdge )
		   break;
	    
	    #ifdef ROBUST
   	       check( (ie>=3), "TBuildDelaunay::InitInflRegnInternal(), <2> this should not happen!" );
	    #endif
	    
	    InflRegnBorder.AddTail( T1->TE[(ie+1)%3] );
	    T1->TE[(ie+1)%3]->Mark( INFL_BORDER );
	    
	    InflRegnBorder.AddTail( T1->TE[(ie+2)%3] );
	    T1->TE[(ie+2)%3]->Mark( INFL_BORDER );
	    
	}   
	
	
    }
    else // PLLocation != PL_TRIANGLE/PL_EDGE...
    {
        error( "TBuildDelaunay::InitInflRegnInternal(), <3> this should not happen!" );
    }

}


// ---------------------------------------------------------------------------------
//  
//  void TBuildDelaunay::InitInflRegnExternal()
//
//  Find the region of influence in case point VertexToIns falls outside
//  the domain of the current triangulation.
//
//  INPUT/OUTPUT: see description of TDelaunay::CalcInfluenceRegion()
//

void TBuildDelaunay::InitInflRegnExternal()
{

    #ifdef DEBUG
     DEBUG << "\nTBuildDelaunay::InitInflRegnExternal()" << endl;
    #endif // DEBUG


    //
    // Start from edge PLEdge, and find the edges of the convex hull that
    // are "visible" from VertexToIns. Then insert them, sorted in
    // counterclockwise order w.r.t.  VertexToIns, in list InflRegnBorder.
    //       
    
    PTVertex VUp = NULL;
    // "upper" starting point to visit the convex hull
    PTVertex VDown = NULL;
    // "lower" starting point to visit the convex hull
    
    //
    // It PLEdge->ET[1] (right adjacent triangle) is NULL, PLEdge is
    // directed counterclockwise on the convex hull.
    // If PLEdge->ET[0] (left adjacent triangle) is NULL, PLEdge is 
    // directed clockwise on the convex hull.
    // The orientation of PLEdge is important to maintain the
    // counterclockwise sorting w.r.t. VertexToIns.
    //
    
    InflRegnBorder.AddHead( PLEdge ); 
    PLEdge->Mark( INFL_BORDER );
        
    if ( PLEdge->ET[1] == NULL )
    {
       
        // only one of the two triangles adjacent to PLEdge is NULL
        
	#ifdef ROBUST
	   check( (PLEdge->ET[0] == NULL),
                  "TBuildDelaunay::InitInflRegnExternal(), <1> this should not happen!" );
	#endif
	
	VUp   = PLEdge->EV[1]; 
	VDown = PLEdge->EV[0];	 
    
    }
    else if ( PLEdge->ET[0] == NULL )
    {
    
        // only one of the two triangles adjacent to PLEdge is NULL

        #ifdef ROBUST
  	   check( (PLEdge->ET[1] == NULL),
                  "TBuildDelaunay::InitInflRegnExternal(), <2> this should not happen!" );
	#endif
		       
	VUp   = PLEdge->EV[0];
	VDown = PLEdge->EV[1];
    
    }
    else // ...( PLEdge->ET[0] != NULL && PLEdge->ET[1] != NULL )
    {
        // 
	// none of the triangles adiacent to PLEdge is NULL, this 
        // contraddicts PLEdge being on the convex hull
	//	
	error( "TBuildDelaunya::InitInflRegnExternal(), <3> this should not happen!");
    }
    
    
    //
    // travel "upwards" on the convex hull, starting from VUp
    // (i.e., in clockwise order w.r.t. VertexToIns)
    //
	
    //
    // Find the next vertex V after VUp, check if ( VertexToIns, VUp, V )
    // define a right turn; in case they do, add edge (V, VUp) to list
    // InflRegnBorder, set VUp = V and go on...
    // otherwise stop upward travel.
    //
	 
    //
    // Start from edge after PLEdge. E is the edge adjacent to VUp (on the
    // convex hull) different from PLEdge. V is the endpoint of E different
    // from VUp. 
    // For vertices on the convex hull (as V and VUp) the two edges stored
    // in relation VE are the ones lying on the convex hull.
    //
	 
    PTEdge E = ( VUp->VE[0] != PLEdge ? VUp->VE[0] : VUp->VE[1] );
    PTVertex V = ( E->EV[0] != VUp ? E->EV[0] : E->EV[1] );
	 
    while ( Geom::Turnxy( VertexToIns, VUp, V ) == TURN_RIGHT )
    {

        InflRegnBorder.AddHead( E );
	E->Mark( INFL_BORDER );
	    
        //
        // next edge and vertex
        //
	    
	VUp = V;
	E = ( VUp->VE[0] != E ? VUp->VE[0] : VUp->VE[1] );
	V = ( E->EV[0] != VUp ? E->EV[0] : E->EV[1] );

    } // end ...while
	 
   
    //
    // now travel "downward", starting from VDown
    //
    
    //
    // Everything is symmetric w.r.t. previous while loop.
    // Find the next vertex V after VDown, check if (VertexToIns, VDown, V)
    // define a left turn; in case they do, add edge (V, VDown) to list
    // InflRegnBorder, set VDown = V and go on...
    // otherwise stop upward travel.
    // During the upward travel edges are added to the end of
    // InflRegnBorder, here they are added to the beginning.
    // This maintains the list sorted counterclockwise w.r.t. VertexToIns.
    //

    E = ( VDown->VE[0] != PLEdge ? VDown->VE[0] : VDown->VE[1] );
    V = ( E->EV[0] != VDown ? E->EV[0] : E->EV[1] );
	  
    while ( Geom::Turnxy( VertexToIns, VDown, V ) == TURN_LEFT )
    {
 	   
        InflRegnBorder.AddTail( E );
        E->Mark( INFL_BORDER );
	      
        VDown = V;
        E = ( VDown->VE[0] != E ? VDown->VE[0] : VDown->VE[1] );
        V = ( E->EV[0] != VDown ? E->EV[0] : E->EV[1] );

    } // end ...while
   

}


// --------------------------------------------------------------------------------
//  
//  void TBuildDelaunay::CalcInflRegnMain()
//
//  Function called by CalcInfluenceRegion after initializing list
//  InflRegnBorder by calling one of InitInflRegnInternal() or 
//  InitInflRegnExternal() which inserts the first elements in list
//  InflRegnBorder. CalcInflRegnMain() performs the main loop to
//  enlarge the region of influence.
//
//  INPUT/OUTPUT: see description of TBuildDelaunay::CalcInfluenceRegion(),
//  an additional input is lista InflRegnBorder initialized by
//  InitInflRegnInternal/External().
//

void TBuildDelaunay::CalcInflRegnMain()
{

   #ifdef DEBUG
    DEBUG << "\nTBuildDelaunay::CalcInflRegnMain()" << endl;
   #endif // DEBUG
 
   PTEdge CurrEdg, Ea, Eb;
   PTTriangle NextTrg;
   PTVertex v0, v1;

   TDoubleListIterator<PTEdge> IRegnIter( &InflRegnBorder );
   IRegnIter.Restart();

   while( !IRegnIter.EndOfList() )
   {
   
      CurrEdg = IRegnIter.Current()->object;
      
      v0 = CurrEdg->EV[0];
      v1 = CurrEdg->EV[1];
      
      // ...piccolo controllo
      #ifdef ROBUST
         check( (v0 == NULL || v1 == NULL),
                "TBuildDelaunay::CalcInflRengMain(), <1> inconsistency detected");
      #endif
           
      //
      // The next triangle is the one adjacent to CurrEdg from the
      // "opposite" side w.r.t. VertexToIns.
      // To decide which side is the opposite one, we look at the
      // orientation of CurrEdg w.r.t. VertexToIns: "left" or "right".
      //
      
      if ( Geom::Turnxy( VertexToIns, CurrEdg->EV[0], CurrEdg->EV[1] ) == TURN_LEFT )
           NextTrg = CurrEdg->ET[1]; // next triangle on the right of CurrEdg
      else
           NextTrg = CurrEdg->ET[0]; // next triangle on the left of CurrEdg      
            
      if ( (NextTrg != NULL) && ( NextTrg->InCircle( VertexToIns ) ) )
      {
      
        int e;

        #ifdef DEBUG
         DEBUG << "marco E" << CurrEdg->EID << " e T" << NextTrg->TID << endl;
	#endif // DEBUG
	 
         CurrEdg->Mark( TO_DELETE );
	 NextTrg->Mark( TO_DELETE );
	 FirstTrgToDel = NextTrg;
	 
	 //
	 // add the edges of NextTrg different from CurrEdg to list
         // InflRegnBorder, and remove CurrEdg from such list
	 //
	 
	 for ( e=0; e<3; e++ )
	    if ( NextTrg->TE[e] == CurrEdg )
	       break;
	 
	 #ifdef ROBUST 
	    check( (e>=3), "TBuildDelaunay::CalcInflRengMain(), <2> inconsistency detected");
         #endif 
	   
	 Ea = NextTrg->TE[(e+1)%3];
	 Eb = NextTrg->TE[(e+2)%3];
	 
	 #ifdef ROBUST
  	    check( (Ea == NULL || Eb == NULL),
	           "TBuildDelaunay::CalcInflRengMain(), <2> inconsistency detected");
         #endif
	 
	 
	 // 
	 // ATTENTION: is it possible that Ea or Eb already belong to
	 //    InflRegnBorder? It seems NO. If Yes, we should check...
	 //
	 
	 #ifdef ROBUST
  	    check( (Ea->Marked( INFL_BORDER) || Eb->Marked( INFL_BORDER )),
	           "TBuildDelaunay::CalcInflRegnMain(), inconsistency detected, edge already marked" );
         #endif
	 
	 
         //
	 // Some care in deleting CurrEdg from InflRegnBorder, to avoid
	 // that the iterator points to the deleted node; and some case
	 // also in inserting Ea and Eb in order to maintain the edge
	 // sequence sorted counterclockwise w.r.t. VertexToIns.
	 //
	
  	 InflRegnBorder.AddBefore( IRegnIter.Current(), Ea );
	 Ea->Mark( INFL_BORDER );
	 
	 InflRegnBorder.AddBefore( IRegnIter.Current(), Eb );
	 Eb->Mark( INFL_BORDER );
	 
	 IRegnIter.GoPrev();
	 
	 InflRegnBorder.RemoveAfter( IRegnIter.Current() );
	 CurrEdg->UnMark( INFL_BORDER );
	 
	 IRegnIter.GoPrev();
	 // now IRegnIter points to Ea
	 
	 	 
      } 
      else // ...if (( NextTrg == NULL ) || ( !InCircle ))
      {
         IRegnIter.GoNext();
      }
         
   
   } // ...end while( !IRegnIter.EndOfList() )
 
   
   #ifdef DEBUG // print InflRegnBorder
     TDoubleListIterator<PTEdge> Iter( &InflRegnBorder );
     Iter.Restart();
     while ( !Iter.EndOfList() )
     {  
        PTEdge EIter = (PTEdge)(Iter.Current()->object);
        DEBUG << "-> E"<< EIter->EID << "(";
        if ( EIter->ET[0] != NULL ) DEBUG << "T" << EIter->ET[0]->TID << "/"; else DEBUG << "NULL/";
        if ( EIter->ET[1] != NULL ) DEBUG << "T" << EIter->ET[1]->TID << ")"; else DEBUG << "NULL)";
        Iter.GoNext();
     }
     DEBUG << endl;
   #endif // DEBUG
 
}


// --------------------------------------------------------------------------------
//  
//  void TBuildDelaunay::RetriangulateInfluenceRegion()
//
//  Re-triangulate the region of influence. Foe each edge E in list
//  InflRegnBorder, built by CalcInfluenceRegion(), create a new triangle
//  whose vertices are VertexToIns and the two endpoints of E.
//  The triangles of the region of influence and the edges inside such 
//  region are deleted.
//  
//  INPUT: list InflRegnBorder
//
//  NOTE (1): list InflRegnBorder contains the sequence of boundary edges
//  of the region of influence of VertexToIns, sorted counterclockwise
//  w.r.t. VertexToIns. Such sequence can be closed or open.
//  It is open if InflRegnBorder has been initialized by 
//  InitInflRegnExternal(), or by  in case InitInflRegnInternal()
//  VertexToIns falls on an edge of the convex hull
//
//  NOTE (2): The triangle pointer FirstTrgToDel may still be NULL.
//  This happens in case VertexToIns is outside the convex hull, and
//  its region of influence contains no triangles.
//
//  OUTPUT: a consistent status of the triangulation representing 
//  the updated triangulation after the insertion of the new point
//

void TBuildDelaunay::RetriangulateInfluenceRegion()
{
    

   #ifdef DEBUG
    DEBUG << "\nTBuildDelaunay::RetriangulateInfluenceRegion()" << endl;
   #endif // DEBUG
   
   
   // --- init variables used in this procedure ---
   
   //
   // first check InflRegnBorder is an open or a closed edge sequence
   //
   
   // first and last edge of the sequence
   PTEdge EFirst = InflRegnBorder.GetHead();
   PTEdge ELast  = InflRegnBorder.GetLast();
   
   PTVertex vf0, vf1, vl0, vl1;

   // vertices of the first edge of the sequence
   vf0 = EFirst->EV[0];
   vf1 = EFirst->EV[1];
   
   // vertices of the last edge of the sequence
   vl0 = ELast->EV[0];
   vl1 = ELast->EV[1];
   
   // check
   
   #ifdef ROBUST
      check( (vf0==NULL || vf1==NULL || vl0==NULL || vl1==NULL ),
   	  "TBuildDelaunay::RetriangulateInfluenceRegion(), inconsistency detected");
   #endif
   	
   //
   // check if InflRegnBorder is open or closed, by testing if
   // EFirst and Elast have a common vertex
   //
   
   //
   // If they have a common vertex, this may happen for two reasons:
   // because InflRegnBorder is a closed sequence,
   // or because it is an open sequence formed by just two edges;
   // variable closed is set to TRUE only in the former case.
   //        
   
   boolean closed;
   
   if (( vf0 == vl0 || vf0 == vl1 ) || ( vf1 == vl0 || vf1 == vl1 ) ) //...un vertice in comune
       closed = ( InflRegnBorder.Lenght() > 2 );
   else
       closed = FALSE;


   // first and last vertex of the sequence, if not closed, 
   // otherwise unique vertex common to the two edges
   
   PTVertex VFirst, VLast;
      
   if ( ! closed )
   {
      VFirst = ( ( Geom::Turnxy( VertexToIns, vf0, vf1 ) == TURN_RIGHT )
                 ? vf1 : vf0 );
		 
      VLast  = ( ( Geom::Turnxy( VertexToIns, vl0, vl1 ) == TURN_RIGHT )
                 ? vl0 : vl1 );
   }
   else
   {
      if ( vf0 == vl0 || vf0 == vl1 ) VFirst = VLast = vf0;
        else
      if ( vf1 == vl0 || vf1 == vl1 ) VFirst = VLast = vf1;   
   }
   
   //
   // Some variables
   //

   // iterator for InflRegnBorder
   TDoubleListIterator<PTEdge> IRegnIter( &InflRegnBorder );
   
   // current edge in list InflRegnBorder,
   // and edge preceding it in such list
   PTEdge   E, EPrev;
   
   // endpoints of current edge
   PTVertex v0, v1;
   
   // edges connecting VertexToIns to the two endpoints of current edge
   PTEdge E1, E2;
   
   // first created edge
   PTEdge E1Start;
   
   
   // --- re-triangulate influence region ---
    
   
   //
   // start traversal of InflRegnBorder to retriangulate the region
   //

   IRegnIter.Restart();

   //
   // at each iteration, E is the current edge on InflRegnBorder,
   // EPrev is the edge of the previous iteration
   //
   
   E = IRegnIter.Current()->object;
   EPrev = NULL;
   
   while( ! IRegnIter.EndOfList() )
   {
   
       //
       // find the two vertices of the current edge, take them in
       // counterclockwise order (w.r.t. VertexToIns).
       //
            
       v0 = E->EV[0];
       v1 = E->EV[1];
      
       if ( Geom::Turnxy( VertexToIns, v0, v1 ) != TURN_LEFT )
       {
          // ...Swap v0 e v1 
          PTVertex vTmp = v0;
	  v0 = v1;
	  v1 = vTmp;
       }
      
       //
       // create twoi edges from VertexToIns to v0, v1
       //      
      
       if ( E == EFirst ) // first iteration of while loop
       {
            // 
  	    // save E1Start, to be used in case the region of influence is
	    // a closed chain (see "else" branch...)
	    //

            #ifdef _GC_ON
               E1Start = E1 = GC::NewEdge( VertexToIns, v0 );
               E2 = GC::NewEdge( VertexToIns, v1 );
	    #else
               E1Start = E1 = new TEdge( VertexToIns, v0 );
               E2 = new TEdge( VertexToIns, v1 );
	       
	       check( (E1 == NULL || E2 == NULL), 
	          "TBuildDelaunay::RetriangulateInfluenceRegion(), insufficient memory" );
	    #endif // _GC_ON
	    
	   	    
       }
       else // next iterations, an edge has already been created
       {
            E1 = E2;
	  
	    if (( closed ) && ( E == ELast )) // close the chain
                E2 = E1Start;
	    else
	    { 
	        #ifdef _GC_ON
		   E2 = GC::NewEdge( VertexToIns, v1 );
		#else
   		   E2 = new TEdge( VertexToIns, v1 ); 
		   check( (E2 == NULL), 
		      "TBuildDelaunay::RetriangulateInfluenceRegion(), insufficient memory" );
		#endif // _GC_ON
		
            }
		
       }
       
       //
       // create triangle ( VertexToIns v0, v1 ), its edges, in
       // counterclockwise order, are E1, E, E2
       //
       
       #ifdef _GC_ON
          PTTriangle NewTriangle = GC::NewTriangle( E1, E, E2 );
       #else
          PTTriangle NewTriangle = new TTriangle( E1, E, E2 );
          check( (NewTriangle==NULL),
             "TBuildDelaunay::RetriangulateInfluenceRegion(), insufficient memory" );
       #endif // _GC_ON
      
       
       AddTriangle( NewTriangle );
             
       // adjust Edge-Triangle relations
               
       if ( E->EV[0] == v0 )          // se E directed "left"
          E->ET[0] = NewTriangle;     //    NewTriangle to the left of E
       else                           // otherwise, if E = ( v1, v0 )
          E->ET[1] = NewTriangle;     //    NewTriangle to the right of E
	  
       //
       // Since VertexToIns is the first vertex of E1, E2,
       // NewTriangle is to the right of E1 and to the left of E2.
       // ET adjacency relations are set to NULL by the constructor of
       // TEdge, thus they remain NULL if not modified. Therefore,
       // for new edges created on the convex hull, ET adjacencies
       // are already correctly set to NULL.
       //
       
       E1->ET[0] = NewTriangle;
       E2->ET[1] = NewTriangle;
       
       // go to the next one
       
       EPrev = E;
       IRegnIter.GoNext();
       
       //
       // Adjust VE adjacencies of ve (common vertex to E and EPrev).
       // If at least one of v1->VE[0] and v1->VE[1] is marked as TO_DELETE,
       // then v1 cannot be on the convex hull.
       // The following instructions, during the execution of the whole
       // while loop, set the VE relation of all vertices of InflRegnBorder,
       // with the exception of VFirst and VLast (they will be considered
       // separately, later).
       //
       
       if ( !IRegnIter.EndOfList() )
       {
           E = IRegnIter.Current()->object;
	  
	   boolean E0Del = v1->VE[0]->Marked( TO_DELETE );
	   boolean E1Del = v1->VE[1]->Marked( TO_DELETE );
	   
	   if ( E0Del && E1Del )
	   {
	       v1->VE[0] = EPrev;
	       v1->VE[1] = E;
	   }
	   else if ( !E0Del && E1Del )
	       v1->VE[1] = ( v1->VE[0] != E ? E : EPrev );
	   else if ( E0Del && !E1Del )
	       v1->VE[0] = ( v1->VE[1] != E ? E : EPrev );
	      
	   // else VE adjacency relations VE od v1 remain unchanged

//DEBUG
//CheckIsolatedPoint(v1);
       
       }
       
   
   } // end ...while( ! IRegnIter.EndOfList() )


   //
   // E2Stop is the LAST created edge, in both cases (closed TRUE or FALSE).
   //
   
   PTEdge E2Stop = ( closed ? E1 : E2 );
      
   
   // --- adjuct VE adjacency relation od VertexToIns, VFirst ans VLast ---
      
   //
   // VE relation of VertexToIns: set to ( E1Start, E2Stop ). This is ok
   // both in the case VertexToIns is on the convex hull (closed=FALSE) and
   // in the other case.
   //
   
   VertexToIns->VE[0] = E1Start;
   VertexToIns->VE[1] = E2Stop;

//DEBUG
//CheckIsolatedPoint(VertexToIns);

   //
   // VE relation of VFirst and VLast. 
   //
   
   if ( !closed )
   {
      //
      // NOTE: if !closed VFirst (as VLast) is CERTAINLY on the convex
      // hull: one of its adjacent edges is still on the hull, and is
      // not affected by the current update step. The other one is
      // certainly marked as TO_DELETE or INFL_BORDER. 
      // Since an edge can be INFL_BORDER, on the convex hull, and not
      // TO_DELETE, we first test if one of the two VE is TO_DELETE,
      // and only later we test if it is INFL_BORDER.
      //
      
      // first VFirst
      if ( VFirst->VE[0]->Marked( TO_DELETE ) )
          VFirst->VE[0] = E1Start;
      else
      if ( VFirst->VE[1]->Marked( TO_DELETE ) )
          VFirst->VE[1] = E1Start;
      else
      if ( VFirst->VE[0]->Marked( INFL_BORDER ) )
          VFirst->VE[0] = E1Start;
      else
      if ( VFirst->VE[1]->Marked( INFL_BORDER ) )
          VFirst->VE[1] = E1Start;
      else
          error( "TBuildDelaunay::RetriangulateInfluenceRegion(), <3> this should not happen");
	  
      // then VLast
      if ( VLast->VE[0]->Marked( TO_DELETE ) )
          VLast->VE[0] = E2Stop;
      else
      if ( VLast->VE[1]->Marked( TO_DELETE ) )
          VLast->VE[1] = E2Stop;
      else
      if ( VLast->VE[0]->Marked( INFL_BORDER ) )
          VLast->VE[0] = E2Stop;
      else
      if ( VLast->VE[1]->Marked( INFL_BORDER ) )
          VLast->VE[1] = E2Stop;
      else
          error( "TBuildDelaunay::RetriangulateInfluenceRegion(), <4> this should not happen");

   }
   else // closed == TRUE
   {
      // VFirst and VLast are the same! Deal with VFirst=VLast as with
      // any other vertex v1 on InflRegnBorder
      
      boolean E0Del = VFirst->VE[0]->Marked( TO_DELETE );
      boolean E1Del = VFirst->VE[1]->Marked( TO_DELETE );
      
      if ( E0Del && E1Del )
      {
         VFirst->VE[0] = EFirst;
	 VFirst->VE[1] = ELast;
      }
      else if ( !E0Del && E1Del )
         VFirst->VE[1] = ( VFirst->VE[0] != EFirst ? EFirst : ELast );
      else if ( E0Del && !E1Del )
         VFirst->VE[0] = ( VFirst->VE[1] != EFirst ? EFirst : ELast );
	 
      // else leave VE relation of VFirst unchanged
      
   } // end ... if (closed)
   
//DEBUG
//CheckIsolatedPoint(VFirst);
//CheckIsolatedPoint(VLast);

   // --- last operations ---

   //
   // delete internal triangles and edges of the influence region
   // 
   
   DeleteInfluenceRegion();
      
   //
   // remove INFL_BORDER mark from edges in InflRegnBorder
   //
   
   IRegnIter.Restart();
   while ( ! IRegnIter.EndOfList() )
   {
       IRegnIter.Current()->object->UnMark( INFL_BORDER );
       IRegnIter.GoNext();
   }
    
   // empty InflRegnBorder, ...it is no longer needed
   
   InflRegnBorder.ClearList();
      
}
