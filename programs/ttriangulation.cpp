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
//   file   : ttriangulation.cpp
//   author : Christian Melchiorre
//
//   Implementation of abstract class TTriangulation.
//


#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>

#include "defs.h"
#include "geom.h"
#include "utils.h" // PAOLA for CheckIsolatedPoint
#include "tdoublelist.h"
#include "ttriang.h"
#include "ttriangulation.h"


#define OUTTIME( tm ) ((tm)/60) << " min. / " << ((tm)%60) << " sec. "

// ---------------------------------------------------------------------------------
//
//   Constructor of class TTriangulation
//
//

TTriangulation::TTriangulation()  
{   
    #ifdef DEBUG
     DEBUG << "Constructor TTriangulation" << endl;
    #endif // DEBUG 

    Points = NULL;
    nPts = 0;

    VertexToIns = NULL;
    iNextPoint = 0;

    VertexToRemove = NULL;
    
    FirstTriangle = NULL;
    nTrg = 0;
    
    PLLocation = PL_UNDEFINED;
    PLTriangle = NULL;
    PLVertex = NULL;
    PLEdge = NULL;
}


// ---------------------------------------------------------------------------------
//
//   void TTriangulation::BuildTriangulation( char *fname )
//
//   The only function implemented by abstract class TTriangulation.
//   Call the other procedures for building the triangulation.
//

void TTriangulation::BuildTriangulation( const char *infname, const char *outfname )
{
   ReadData( infname );

   long starttime = time(NULL);
   long prevtime  = starttime;
   long currtime  = starttime;

   InitialTriangulation();
   
   int np = iNextPoint;
   
   cerr << endl;

   while ( !NoMoreUpdates() )
   {
         // ...give a signal each second
         currtime = time(NULL);
         if ( currtime != prevtime ) 
         {
	      prevtime = currtime;
              cerr << "processed " << np << " points   \r" << flush;
	 }
	 
         UpdateStep();
	 
	 np++;
   }
	 
   cerr << "\rprocessed " << np << " points                   " << endl;

   long stoptime = time(NULL);

   currtime = (stoptime-starttime);
   cerr << "running time for the triangulation: " << OUTTIME( currtime) << endl;

   WriteData( outfname );

   EndTriangulation();
  
}


// ---------------------------------------------------------------------------------
//
//   void TTriangulation::PointLocation ( PTPoint PointToLoc )
//
//   OUTPUT: set the value of status variables PLLocation, PLVertex, PLEdge,
//           and PLTriangle, as better explained in the following.
//
//   Locate point PointToLoc in the current triangulation. The result of
//   calling this function consists in setting the value of 4 status
//   variables of class TTriangulation:
//   
//     int PLLocation = PL_EXTERNAL if point is outside the triangulation,
//                      PL_TRIANGLE if point is inside a triangle,
//                      PL_EDGE if point is inside an edge,
//                      PL_VERTEX if it is coincident with a vertex
//
//     PTVertex PLVertex = if PLLocation=PL_VERTEX, point to coincident vertex
//	otherwise NULL. (this is needed for RepositionPoint in dechoppe.cpp).
//
//     PTTriangle PLTriangle = if PLLocation=PL_TRIANGLE, point to
//      containing triangle, otherwise NULL.
//
//     PTEdge PLEdge = if PLLocation=PL_EDGE, point to containing edge,
//      if PLLocatio=PL_EXTERNAL point to an edge of the convex hull which 
//      "visible" from point VertexToIns,
//      if PLLocation = PL_TRIANGLE, then NULL.
//

void TTriangulation::PointLocation ( PTPoint PointToLoc )
{
   
   #ifdef DEBUG
    SKIPPO
    DEBUG << "TTriangulation::PointLocation(" << (*PointToLoc) << ")" << endl;
   #endif // DEBUG

   PTVertex V[3];
   int i;
   PLVertex = NULL; // for default PLVertex=NULL; only in one case it is not

   //
   // start searching from first triangle (convention)
   //
   
   PLTriangle = FirstTriangle;
   int NrPLTriangle = 1;   
   
   while( TRUE )
   {
   
      //
      // find vertices of PLTriangle, and check if PointToLoc is equal
      // to one of the them (in such case, exit). Recall that after
      // calling GetTV, for each i=0..3, 
      // PLTriangle->TE[i] is the edge of endpoints V[i], V[(i+1)%3]
      //
      
      PLTriangle->GetTV( V[0], V[1], V[2] );
      
      for ( i=0; i<3; i++ )
      {
          //
          // if xy coordinates or PointToLoc are the same as V[i]...
	  //
          if ( PointToLoc->Equalsxy( V[i] ) )
	  {
	     PLLocation = PL_VERTEX;
	     PLVertex = V[i];
	     PLTriangle = NULL;
             PLEdge = NULL;
	     
	     return;
	  }
      }	
	
      //
      // find (by using turns) if PointToLoc is inside or outside
      // the boundary of the triangle
      //
            
      int external = -1;
      int aligned  = -1;
      
      for( i=0; i<3; i++ )
      {
         switch( Geom::Turnxy( V[i], V[(i+1)%3], PointToLoc ) )
	 {

             case ALIGNED: // allineato
	     
	          aligned = i;
		  break;
		  
	     case TURN_LEFT: // interno
	       
	          break;
	     
	     case TURN_RIGHT: // esterno
	     		  
		  external = i;
		  break;
		  
         } // end ...switch
		  
       } // end ...for
	 
       
       //
       // check results of point location w.r.t. current triangle
       //

       if ( external == -1 && aligned == -1 )
       {
       
          // 
	  // it means point inside triangle
	  //
	  
	  PLLocation = PL_TRIANGLE;
	  PLEdge = NULL;
	  // PLTriangle remains set to the current triangle
	  return;
        
	}
	else if ( external == -1 && aligned != -1 )
	{
	
	   //
	   // it means point inside edge V[aligned], V[(aligned+1)%3],
           // i.e., edge PLTriangle->TE[aligned].
	   // 
	   
	   PLLocation = PL_EDGE;
	   PLEdge = PLTriangle->TE[aligned];
	   PLTriangle = NULL;
	   return;
	   
	}
	else // external != -1
	{
	
	    //
	    // PointToLoc has been found outside edge
	    // PLTriangle->TE[external], continue searching in 
	    // adjacent triangle to PLTriangle in that direction
	    //
	    	    
	    PTTriangle PLNext = PLTriangle->GetTT( external );
	    	    
	    if (PLNext == NULL)
	    {
	       //
	       // no triangle in this direction! PointToLoc is
	       // outside convex hull
	       //
	       	       
	       PLLocation = PL_EXTERNAL;
	       PLEdge = PLTriangle->TE[external];
               PLTriangle = NULL;
	       return;
	    }
	    else 
	    {   PLTriangle = PLNext;       
		NrPLTriangle++;

		// if NrPLTriangle>nTrg then we have traversed more
                // triangles than the existing ones: we are in an infinite
		// loop, there must be an inconsistency somewhere!
   		check ((NrPLTriangle>nTrg), "TTriangulation::PointLocation: inconsistency !");
	    }
	       
	} // end ...(external !=-1)
	 
	
   } // end ...while(TRUE)
   
}



// --------------------------------------------------------------------------------
//  
//  void TTriangulation::DetachXXX( PTXXX )
//
//  Detach an entity (vertex, edge, triangle) from the triangulation.
//  To be called before deleting such entity. 
//  Remove pointers pointing to the detached entity from other 
//  entities adjacent to it.
//
//  ATTENTION: this function leaves the triangulation in inconsistent 
//  status. It must be called from other functions that guarantee the
//  consistency of the final status.
//

void TTriangulation::DetachVertex( PTVertex V )
{ 
   int e, v;
   for ( e=0; e<2; e++ )
     for ( v=0; v<2; v++ )
      if ( V->VE[e] != NULL && V->VE[e]->EV[v] == V ) V->VE[e]->EV[v] = NULL;
}


void TTriangulation::DetachEdge( PTEdge E )
{
  int v, e, t; 
 
  #ifdef DEBUG
   DEBUG << endl << "DETACHED E" << E->EID 
         << ": ( V" << E->EV[0]->VID << ", V" << E->EV[1]->VID << " )" << endl;
  #endif // DEBUG 
  
  #ifdef DEBUG
   for ( v=0; v<2; v++ )
      for ( e=0; e<2; e++ ) 
        if ( E->EV[v] != NULL && E->EV[v]->VE[e] == E )
	{
           DEBUG << "Warning, detaching E" << E->EID
	         << ", e V" << E->EV[v]->VID << "is VE-Adjacent to it" << endl;
           error("");
	}
  #endif // DEBUG

   for ( v=0; v<2; v++ )
     if ( E->EV[v] != NULL )
       if ( E->EV[v]->VE[0] == E ) E->EV[v]->VE[0] = NULL;
          else if ( E->EV[v]->VE[1] == E ) E->EV[v]->VE[1] = NULL;
	
   for ( t=0; t<2; t++ )
     if ( E->ET[t] != NULL ) 
       for ( e=0; e<3; e++ )
         if ( E->ET[t]->TE[e] == E ) E->ET[t]->TE[e] = NULL;

//   if (CheckIsolatedPoint(E->EV[0]) || CheckIsolatedPoint(E->EV[1]))
//      cout << "DetachEdge" << endl;
// it is normal that one of the points is isolated (the detached vertex)
}


void TTriangulation::DetachTriangle( PTTriangle T )
{
  int e;

  nTrg--;

  for ( e=0; e<3; e++ )
    if ( T->TE[e] != NULL )
       if ( T->TE[e]->ET[0] == T ) T->TE[e]->ET[0] = NULL;
         else if ( T->TE[e]->ET[1] == T ) T->TE[e]->ET[1] = NULL;

  #ifdef DEBUG
    DEBUG << endl << "DETACHED T" << T->TID << "(ntrg = " << nTrg << ")" << endl;
  #endif // DEBUG
  
}


// --------------------------------------------------------------------------------
//  
//  void TTriangulation::WriteData( const char * )
//
//  Save final triangulation to file.
//

void TTriangulation::WriteData( const char *outfname )
{

    int iv, it, i;

     #ifdef DEBUG
      DEBUG << "\nTTriangulation::WriteData(" << outfname << ")" << endl;
     #endif // DEBUG

    //
    // create arrays to contain output vertices and triangles 
    //
    
    // PAOLA: number of really existing vertices
    int nVrt = 0;

    PTVertex *VtxArray;   
    PTTriangle *TrgArray;

    VtxArray = new PTVertex[ nPts ];
    TrgArray = new PTTriangle[ nTrg ];
        
    check( ( VtxArray == NULL || TrgArray == NULL ),
           "TTriangulation::WriteData(), insufficient memory" );
	  
    for( iv=0; iv<nPts; iv++ ) VtxArray[iv] = NULL;
    for( it=0; it<nTrg; it++ ) TrgArray[it] = NULL;
    
    //
    // open output file
    //
           
    ofstream outFile;
    
    outFile.open( outfname );
    
    outFile.precision( 12 );
   
    TDoubleList<PTTriangle> Triangles;
   
    check( (FirstTriangle==NULL), "TTriangulation::WriteData(), No triangles?");
   
    Triangles.AddHead( FirstTriangle );
    FirstTriangle->Mark( VISITED );
   
    PTVertex v[3];
   
    PTTriangle NextTrg[3];
    PTTriangle CurTrg;
   
    it = 0;
    iv = 0;

    while( !Triangles.IsEmpty() )
    {
       
       CurTrg = Triangles.RemoveHead();
       
       TrgArray[it++] = CurTrg;

       // ...vertices of CurTrg
      
       CurTrg->GetTV( v[0], v[1], v[2] );
      
       #ifdef ROBUST
          check( (v[0] == NULL || v[1] == NULL || v[2] == NULL),
              "TTriangulation::WriteData(), inconsistency detected");
       #endif
      
       // ...adjacent triangles
       
       CurTrg->GetTT( NextTrg[0], NextTrg[1], NextTrg[2] );
       
       for( i=0; i<3; i++ )
       {
         
	 //
	 // put in VtxArray the pointers to the vertices, sorted by VID
	 //
	 
         if ( VtxArray[v[i]->VID] == NULL )
         {
              iv++;
	      VtxArray[v[i]->VID] = v[i];
         }
	 
	 //
	 // enqueue adjacent triangles of CurTrg
	 //
	      
         if ( NextTrg[i] != NULL )
            if (!( NextTrg[i]->Marked( VISITED ) ))
	    {
                 Triangles.AddTail( NextTrg[i] );
		 NextTrg[i]->Mark(VISITED);
	    }
       
       } // end ...for(i=0..2)
            
    } // end ...while( !Triangles.IsEmpty() )
     
    
    //
    // output ...first vertices...
    //
    
    // PAOLA: the vertex array may have holes, i.e., positions that are
    // NULL since the vertex with that VID has been deleted from the
    // triangulation. Scan the array and re-number VIDs of existing 
    // vertices with consecutive numbers, jumping over holes.
    // The number of non-NULL vertices must be equal to iv.
    // 
    nVrt = iv;
    {  int new_vid = 0;
       for ( iv=0; iv<nPts; iv++ )
       {
         if (VtxArray[iv]) VtxArray[iv]->VID = new_vid++;
       }
       check( (new_vid!=nVrt),
              "TTriangulation::WriteData(), error in vertex number");
    }

    // PAOLA: now we can simply print the non-NULL vertices.
    // Renumbering of VIDs guarantees that, when we print triangles,
    // vertex indices (which are VIDs) refer existing vertices.
    outFile << nVrt << endl;

    for( iv=0; iv<nPts; iv++ )
    {
       PTVertex CV = VtxArray[iv];
       
       if ( CV != NULL )
          outFile << CV->x << '\t' << CV->y << '\t' << CV->z << endl;
      
       if ( iv % 1000 == 0 ) cerr << "output " << iv << " vertices\r";
    }
    
    cerr << "output " << nVrt << " vertices from " <<
             nPts << " original points" << endl;
    
    //
    // ...then triangles
    //
    
    outFile << nTrg << endl;

    for( it=0; it<nTrg; it++ )
    {
       PTTriangle CT = TrgArray[it];

       check( (CT == NULL), "TTriangulation::WriteData(), <2> inconsistency detected");

       CT->GetTV( v[0], v[1], v[2] );
       
       outFile << v[0]->VID << '\t' << v[1]->VID << '\t' << v[2]->VID << endl;
       
       if ( it % 1000 == 0 ) cerr << "output " << it << " triangles\r";

       CT->UnMark( VISITED ); // unmark triangles
    }

    cerr << "output " << nTrg << " triangles" << endl;
    
    outFile.close();
    
}


// --------------------------------------------------------------------------------
//  
//  void TTriangulation::ConvertData(int *vNum, int *tNum, int *eNum,
//                                   float **vData, int **tData, int **eData)
//
//  Return vertices, triangles, and constraint edges (if present)
//  in indexed format. Put them in three arrays: vData (vertices),
//  tData (triangles), eData (constraints) which are passed as parameters.
//  Also return the number of vertices, triangles, constraints
//  in vNum, tNum, eNum.
//  If a number is zero, the corresponding array can be null.
//  Arrays are freed (if not null) and re-allocated inside this function.

void TTriangulation::ConvertData(int *vNum, int *tNum, int *eNum,
                                 float **vData, int **tData, int **eData)
{
    int iv, it, i;

     #ifdef DEBUG
      DEBUG << "\nTTriangulation::ConvertData(" << outfname << ")" << endl;
     #endif // DEBUG

    //
    // create arrays to contain output vertices and triangles 
    //
    
    // PAOLA: number or really existing vertices
    int nVrt = 0;

    PTVertex *VtxArray;   
    PTTriangle *TrgArray;

    if (*vData) free(*vData);
    *vData = NULL;
    if (*tData) free(*tData);
    *tData = NULL;
    if (*eData) free(*eData);
    *eData = NULL;
    *vNum = *tNum = *eNum = 0;
    if (nPts==0) return;

    VtxArray = new PTVertex[ nPts ];
    TrgArray = new PTTriangle[ nTrg ];
        
    check( ( VtxArray == NULL || TrgArray == NULL ),
           "TTriangulation::ConvertData(), insufficient memory" );
	  
    for( iv=0; iv<nPts; iv++ ) VtxArray[iv] = NULL;
    for( it=0; it<nTrg; it++ ) TrgArray[it] = NULL;
    
    TDoubleList<PTTriangle> Triangles;
   
    check( (FirstTriangle==NULL), "TTriangulation::ConvertData(), No triangles?");
   
    Triangles.AddHead( FirstTriangle );
    FirstTriangle->Mark( VISITED );
   
    PTVertex v[3];
   
    PTTriangle NextTrg[3];
    PTTriangle CurTrg;
   
    it = 0;
    iv = 0;

    while( !Triangles.IsEmpty() )
    {
       
       CurTrg = Triangles.RemoveHead();
       
       TrgArray[it++] = CurTrg;

       // ...vertices of CurTrg
      
       CurTrg->GetTV( v[0], v[1], v[2] );
      
       #ifdef ROBUST
          check( (v[0] == NULL || v[1] == NULL || v[2] == NULL),
              "TTriangulation::WriteData(), inconsistency detected");
       #endif
      
       // ...adjacent triangles
       
       CurTrg->GetTT( NextTrg[0], NextTrg[1], NextTrg[2] );
       
       for( i=0; i<3; i++ )
       {
         
	 //
	 // put in VtxArray the pointers to the vertices, sorted by VID
	 //
	 
         if ( VtxArray[v[i]->VID] == NULL )
         {
              iv++;
	      VtxArray[v[i]->VID] = v[i];
         }
	 
	 //
         // enqueue adjacent triangles of CurTrg
	 //
	      
         if ( NextTrg[i] != NULL )
            if (!( NextTrg[i]->Marked( VISITED ) ))
	    {
                 Triangles.AddTail( NextTrg[i] );
		 NextTrg[i]->Mark(VISITED);
	    }
       
       } // end ...for(i=0..2)
            
    } // end ...while( !Triangles.IsEmpty() )
     

    // PAOLA: the vertex array may have holes, i.e., positions that are
    // NULL since the vertex with that VID has been deleted from the
    // triangulation. Scan the array and re-number VIDs of existing 
    // vertices with consecutive numbers, jumping over holes.
    // The number of non-NULL vertices must be equal to iv.
    // 
    nVrt = iv;
    {  int new_vid = 0;
       for ( iv=0; iv<nPts; iv++ )
       {
         if (VtxArray[iv]) VtxArray[iv]->VID = new_vid++;
       }
       check( (new_vid!=nVrt),
              "TTriangulation::WriteData(), error in vertex number");
    }
    
    //
    // output data ...first of all allocate memory
    //

    *vNum = nVrt;
    *tNum = nTrg;
    *eNum = 0;
    *vData = (float*)malloc(3*nVrt*sizeof(float));
    *tData = (int*)malloc(3*nTrg*sizeof(int));

    //
    // ...first vertices...
    //

    // PAOLA: now we can simply print the non-NULL vertices.
    // Renumbering of VIDs guarantees that, when we print triangles,
    // vertex indices (which are VIDs) refer existing vertices.

    i = 0;
    for( iv=0; iv<nPts; iv++ )
    {
       PTVertex CV = VtxArray[iv];
      
       if ( CV != NULL )
       {  (*vData)[3*i] = CV->x;
          (*vData)[3*i+1] = CV->y;
          (*vData)[3*i+2] = CV->z;
          i++;
       }
       if ( iv % 1000 == 0 ) cerr << "converted " << iv << " vertices\r";
    }
    
    //
    // ...then triangles
    //
    
    for( it=0; it<nTrg; it++ )
    {
       PTTriangle CT = TrgArray[it];

       check( (CT == NULL), "TTriangulation::ConvertData(), <2> inconsistency detected");

       CT->GetTV( v[0], v[1], v[2] );
       
       (*tData)[3*it] = v[0]->VID;
       (*tData)[3*it+1] = v[1]->VID;
       (*tData)[3*it+2] = v[2]->VID;
       
       if ( it % 1000 == 0 ) cerr << "converted " << it << " triangles\r";

       CT->UnMark( VISITED ); // unmark triangles
    }

}
