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
//   file : refCDT.cpp
//   author : Alessio Calcagno
//

#include <stdlib.h>
#include "refCDT.h"
#include "geom.h"


#ifdef OUTPUT
  #include <time.h>
#endif


/**********************                   *********************************/



/**********************    RefCDTDelaunayTRA    ****************************/

//
// Insert in the triangulation all constraint segments, without assuming 
// that they are admissible (a constraint is admissible if no vertex 
// lies on it and if it does not intersect other constraints).
// Therefore, check for admissibility of a constraint segment and, if it
// is admissible, insert it.
//

void TRefCDT :: Add_Constrain()
{
   #ifdef DEBUG14
      DEBUG14 << "\n\nTRefCDT::Add_Constrain()\n";
   #endif

   int i;

   nConstr *= 2; 

   for (i=0; i<nConstr; i += 2)
   {
      #ifdef DEBUG14
         cerr << "TRefCDT::Add_Constrain() : constraint insertion loop:" << endl
              << "insert constraint C" << i/2 << "( " <<  constr_array[i] << ", "<< constr_array[i+1]
              << " ) = ( " << * (PTVertex) Points[constr_array[i]] << ", "
              <<  * (PTVertex) Points[constr_array[i+1]] << " )\n";
         Pause("");
      #endif
                 
      if (constr_array[i] == constr_array[i+1])
      {
         cerr << "Warning: constraint C" << i/2 << "( " << constr_array[i] << ", "<< constr_array[i+1] 
              << " ) is degenerate: not inserted\n";
         continue;
      }

      if ( constr_array[i] < 0   || constr_array[i] >= nPts   ||
           constr_array[i+1] < 0 || constr_array[i+1] >= nPts    )
      {
         cerr << "Warning: constraint C" << i/2 << "( " << constr_array[i] << ", "<< constr_array[i+1] 
              << " ) has non-existing vertices " << "( nPts = " << nPts << " ) : not inserted\n";
         continue;
      }

      #ifdef PAOLO
         cerr << "\n\n\n\n INSERT CONSTRAINT " << i/2 << " :  " << constr_array[i] << " "<< *((PTVertex)Points[constr_array[i]]);
         cerr << "    " << constr_array[i+1] <<" "<< *((PTVertex)Points[constr_array[i+1]]);
      #endif

      // check for admissibility of current constraint

      if( IsEdgeAdmissible( (PTVertex) Points[constr_array[i]], (PTVertex) Points[constr_array[i+1]] ) )
      {
         TTra_Add_Constraint( (PTVertex) Points[constr_array[i]], (PTVertex) Points[constr_array[i+1]] );

         #ifdef DEBUG
            DEBUG << "constraint inserted: ( " << Input_Constraints[i] << ", "
                  << Input_Constraints[i+1] << " )\n";
         #endif
      }
      else
      {
         cerr << "\nConstraint ( " << *((PTVertex)Points[constr_array[i]])
              << "cannot be inserted, " << *((PTVertex)Points[constr_array[i+1]]) << " ) because not admisible\n";
      }
//      if ( i % 1000 == 0 ) cerr << "\rread " << i << " constraints" << flush;  
      cerr << "\rread " << i << " constraints" << flush;  
   }  // endfor
   
   nConstr /= 2;
}



void TRefCDT::Delete_Constr_InfluenceRegion( PTEdge newConstr )
{

   #ifdef DEBUG
     DEBUG << "TRefCDT::Delete_Constr_InfluenceRegion()" << endl;
   #endif
    
   //
   // delete old triangles
   //
   
   if ( ! TDecimDelaunay::InitialPhase ) TDecimDelaunay::MT_KillInterference();
    
   //
   // Remove newConstr from Left_Border and Right_Border, and append the
   // two lists in InflRegnBorder in such a way that they form one
   // chain of edges bounding the region of influence of newConstr
   //

   check( newConstr != Right_Border.RemoveHead(), "TRefCDT::Delete_Constr_InfluenceRegion(..) : <1> inconsistency detected\n" );
   check( newConstr != Left_Border.RemoveHead(), "TRefCDT::Delete_Constr_InfluenceRegion(..) : <2> inconsistency detected\n" );
   #ifdef ROBUST
      check( ! InflRegnBorder.IsEmpty(), "TRefCDT::Delete_Constr_InfluenceRegion(..) : <3> inconsistency detected\n" );
      int nEdges = Right_Border.Lenght() + Left_Border.Lenght();
   #endif
   InflRegnBorder.Append( Right_Border );
   InflRegnBorder.Append( Left_Border );
   #ifdef ROBUST
      check( nEdges != InflRegnBorder.Lenght(), "TRefCDT::Delete_Constr_InfluenceRegion(..) : <4> inconsistency detected\n" );
      check( 0 != Right_Border.Lenght(), "TRefCDT::Delete_Constr_InfluenceRegion(..) : <5> inconsistency detected\n" );
      check( 0 != Left_Border.Lenght(), "TRefCDT::Delete_Constr_InfluenceRegion(..) : <6> inconsistency detected\n" );
   #endif


   //
   // call old version of DeleteInfluenceRegion() to delete the triangles
   //
   
//   TDestroyDelaunay::DeleteInfluenceRegion();
   TDelaunayBase::DeleteInfluenceRegion();

/*
   // 
   // TDestroyDelaunay::DeleteInfluenceRegion() scans the list
   // InflRegnBorder and re-considers the removability of each vertex V
   // (VertexToRemove can be among such vertices), for each vertex:
   //   - delete it from ElimVtxTree,
   //   - perform RecheckVertex to check its removability again
   //   - if Recheck returns true, insert it in ElimVtxTree.
   //  Obviously VertexToRemove must not be inserted in ElimVtxTree:
   //  thus, if it has been inserted, we remove it.
   //
   if( ElimVtxTree.IsIn( VertexToRemove ) )
      ElimVtxTree.Remove( VertexToRemove );
*/

   //
   // Reposition the "Detached" points from deleted triangles into new
   // triangles, in such a way that we can compute the error 
   // correctly later in MT_AddComponent().
   //
//   TDecimDelaunay::RepositionPoint( VertexToRemove ); 
      // here we have not retriangulated the region of influence of
      // VertexToRemove: it must not be repositioned

   while( ! TDecimDelaunay::DetachedPoints.IsEmpty() )
   {
       PTPoint p = TDecimDelaunay::DetachedPoints.RemoveHead();
       TDecimDelaunay::RepositionPoint( p );       
   } 

   //
   // adjust new triangles (those marked as NEW_TRIANGLE) and unmark them
   //
   
   if ( ! TDecimDelaunay::InitialPhase ) MT_AddComponent();


   //
   // remove INFL_BORDER mark from edges in InflRegnBorder
   //
      
   TDoubleListIterator<PTEdge> IRegnIter( &InflRegnBorder );
   IRegnIter.Restart();
   while ( ! IRegnIter.EndOfList() )
   {
       IRegnIter.Current()->object->UnMark( INFL_BORDER );
       IRegnIter.GoNext();
   }
      
   // empty InflRegnBorder, is is no longer needed
   
   InflRegnBorder.ClearList();


   //
   // remove INFL_BORDER mark from newConstr
   //
   newConstr->UnMark( INFL_BORDER );

}


/**********************    RefCDTDelaunay    ****************************/

//
// Call the version of InitialTriangulation() provided by TRefineDelaunay.
// Such function finds which points belong to the convex hull of the set
// of input points, and compute an initial triangulation.
// After executing TRefineDelaunay::InitialTriangulation(), the pointers
// to the points of the convex hull are in positions [0...nChPts-1] of 
// array Points, and other points (not on the convex hull) are in the 
// remaining positions [nChPts...nPts-1].
//

void TRefCDT :: InitialTriangulation()
{
    
    register int i;

    #ifdef PAOLO
       cerr << "TRefCDT::InitialTriangulation()" << endl;
    #endif 
   
    TRefineDelaunay::InitialPhase = TRUE;


    //
    // compute convex hull of set of input points
    //

    TRefineDelaunay::CalcConvexHull();

    #ifdef DEBUG
       DEBUG << " TRefCDT::convex hull computed" << endl;
    #endif
    
    //
    // After calling CalcConvexHull(), the pointers to the vertices of
    // the convex hull are found in  positions  [0...nChPts-1] of array
    // Points. We insert them in the initial triangulation.
    // By calling method InitialTriangulation() of base class TBuildDelaunay,
    // we build an initial triangulation formed by the first three
    // non-aligned points in Points.
    // Due to the position of the points in array Points (first the ones
    // on the convex hull, then the others), the first three non-aligned
    // points certainly belong to the convex hull.
    //


   //
   // Search for the first three non-aligned points in Points[]. 
   // Then, store them in Points[0], Points[1] e Points[2], in 
   // counterclockwise order.
   //

   int ip2 = 2;

   while( Geom::Alignedxy( Points[0], Points[1], Points[ip2] ) )
   {

      ip2++;
      check( (ip2 >= nPts), " TRefCDT::InitialTriangulation(), all points are aligned" );
   }

   // move Points[ip2] in 3rd position in the array (if not already there)
    
   if ( ip2 != 2 )
   {  
      PTPoint pTmp = Points[ip2];
      Points[ip2] = Points[2];
      Points[2] = pTmp;
   }
   

   //
   // make Points[0], Points[1], Points[2] a counterclockwise sequence
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
   // create three vertices corresponding to the three points, and then
   // delete the three points
   //
      
   PTVertex v0 = new TVertex( Points[0] ); 
   PTVertex v1 = new TVertex( Points[1] ); 
   PTVertex v2 = new TVertex( Points[2] ); 
   
   check( ( v0 == NULL || v1 == NULL || v2 == NULL ),
	   " TRefCDT::InitialTriangulation(), insufficient memory" );
    
   
   OrderInitial[Points[0]->PID] = 0;
   delete Points[0]; Points[0] = v0;

   OrderInitial[Points[1]->PID] = 1;
   delete Points[1]; Points[1] = v1;

   OrderInitial[Points[2]->PID] = 2;
   delete Points[2]; Points[2] = v2;
   
   iNextPoint = 3;
   
   //
   // create the three edges of triangle formed by vertices v0,v1,v2,
   // in counterclockwise order
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
         " TRefCDT::InitialTriangulation(), insufficient memory" );
   #endif // _GC_ON
   
   //
   // now create the triangle
   //
   
   #ifdef _GC_ON
      PTTriangle t0 = GC::NewTriangle( e0, e1, e2 );
   #else
      PTTriangle t0 = new TTriangle( e0, e1, e2 );
      check( (t0 == NULL),  " TRefCDT::InitialTriangulation(), insufficient memory" );
   #endif // _GC_ON   
   
   
   //
   // we still have to set "inverse" relations, i.e., adjacency relations
   // from one entity to another entity of larger dimension. Such relations
   // are not initialized by the constructors.
   //

   v0->VE[0] = e2;    v0->VE[1] = e0;
   v1->VE[0] = e0;    v1->VE[1] = e1;
   v2->VE[0] = e1;    v2->VE[1] = e2;

   e0->ET[0] = t0;    e0->ET[1] = NULL;   
   e1->ET[0] = t0;    e1->ET[1] = NULL;
   e2->ET[0] = t0;    e2->ET[1] = NULL;
   
   
   AddTriangle( t0 );
   
       CheckIsolatedPoint(v0);
       CheckIsolatedPoint(v1);
       CheckIsolatedPoint(v2);

  //
  //  the initial triangulation of the first three points has been built,
  //  now add the remaining points of the convex hull
  //


    for ( i=3; i<nChPts; i++ )
    { 
       VertexToIns = new TVertex( Points[i] );
       check( (VertexToIns == NULL), "TRefineDelaunay::InitialTriangulation(), insufficient memory");
       

       OrderInitial[Points[i]->PID] = VertexToIns->VID;
       delete (Points[i]);
       Points[i] = VertexToIns;

       
       TRefineDelaunay::InsertVertex();
       #ifdef OUTPUT
         if ( i % 1000 == 0 ) cerr << "inserted " << i << " points in the convex hull   \r" << flush;
       #endif
    }

    #ifdef OUTPUT
      cerr << "inserted " << nChPts << " points in the convex hull" << endl;
    #endif
    

    //
    // the remaining points after building the triangulation of the convex
    // hull, i.e., Points[nChPts]...Points[nPts-1], are inserted in the
    // pointlists of edges/triangles inside which they fall
    //
    
    //for ( i=nChPts; i<nPts; i++ )
    //     RepositionPoint( Points[i] );
       
    //
    // the initial triangulation has been built, the next point to be
    // inserted is the one in position nChPts in array Points
    //

    iNextPoint = nChPts;
    #ifdef OUTPUT
      cout << "\ntriangulation of convex hull completed" << endl;
    #endif    

    MT_Initial();
    
    TRefineDelaunay::InitialPhase = FALSE;
 
}



//
//   void TRefCDT :: NextPoint()  derived from  TBuildDelaunay:: NextPoint()
//
//   Put into variable VertexToIns the pointer to the next point to be
//   inserted in the triangulation.
//
//   INPUT: iNextPoint, index in array Points of the first uninserted point
//
//   OUTPUT: VertexToIns, the next point to be inserted
//
//   SIDE EFFECT: increment the value of iNextPoint by one,
//   update array OrderInitial, which maintains the initial order 
//   of the points
//

void TRefCDT :: NextPoint()
{
   PTPoint PointToIns = Points[iNextPoint];
   
   #ifdef ROBUST
     check( (PointToIns == NULL), "TBuildDelaunay::NextPoint(), inconsistency detected" );
   #endif

   #ifdef DEBUG
       cout << "TRefCDT:: NextPoint()" << endl;
   #endif


   VertexToIns = new TVertex( PointToIns );
   OrderInitial[PointToIns->PID] = VertexToIns->VID;
   Points[iNextPoint] = VertexToIns;
   check( (VertexToIns == NULL), "TBuildDelaunay::NextPoint(), insufficient memory");

   //
   // Maintaining the vertices in array Points, in order of their creation,
   // is useful for function WriteData(), when we write them in order
   // of creation (in such a way that references from a triangle to the 
   // VIDs of its vertices are coincident with the order in which
   // vertices are written).
   //

   delete PointToIns;

   iNextPoint++;
}



void TRefCDT :: ReadData( const char *infname )
{
   register int i;

   #ifdef DEBUG
    DEBUG << "\nTRefCDT::ReadData()" << endl;
   #endif // DEBUG
    
   ifstream inFile;
   
   inFile.open( infname );
   
   check( !inFile, "TRefCDT::ReadData(), cannot open input file");

   //
   // read number of points
   //
   
   inFile >> nPts;  
   
   check( (nPts < 3), "input with less than 3 points" );
   
   //
   // read input points and store them in array Points
   //
   
   OrderInitial = new int[nPts];
   Points = new PTPoint[nPts];
   check( (Points == NULL), "TRefCDT::ReadData(), insufficient memory for Points");
   check( (OrderInitial == NULL), "TRefCDT::ReadData(), insufficient memory for OrderInitial");
   
   int ip = 0;
   PTPoint p = NULL;
   
   
   for( i=0; i<nPts; i++ )
   {
       p = new TPoint();
       p->PID = i;
      
       check( (p == NULL), "TBuildDelaunay::ReadData(), insufficient memory for TPoint");
       check( (inFile.eof()), "TBuildDelaunay::ReadData(), unexpected End Of File");
      
       //
       // read point (its three coordinates)
       //
     
       inFile >> (*p); 
       
       #ifdef DEBUG
        cerr << "letto punto: " << (*p) << endl;
       #endif // DEBUG
       
       //
       // add point to the end of list Points
       //
       //OrderInitial[ip] = p;
       Points[ip++] = p; 
       #ifdef OUTPUT
         if ( ip % 1000 == 0 ) cerr << "\rread " << ip << " points" << flush;  
       #endif
   }
   
   #ifdef OUTPUT
     cerr << "\rread " << ip << " points" << endl;
   #endif

   //
   // read constraints (as pairs of indices), store each pair
   // in constr_array[i] and constr_array[i+1]
   //

   inFile >> nConstr;

   if (nConstr>0)
   {
     int count=0;
     constr_array = new int[nConstr*2];
     check( (constr_array == NULL), "TBuildDelaunay::ReadData(), insufficient memory for constr_array");

     nConstr *= 2;
     for( i=0; i<nConstr; /* && !inFile.eof();*/ )
     {
       inFile >> ip;
       constr_array[i++] = ip;
       inFile >> ip;
       constr_array[i++] = ip;
       count++;
     }

     nConstr = nConstr/2; 
     #ifdef OUTPUT
       cerr << "\n read " << nConstr << " constraints  \n";
     #endif
   }
      
   inFile.close();
   
}



void TRefCDT :: BuildTriangulation( const char *infname, const char *outfname )
{
   ReadData( infname );

#ifdef OUTPUT
   long starttime = time(NULL);
   long prevtime  = starttime;
   long currtime  = starttime;
#endif

   InitialTriangulation();
   
   #ifdef OUTPUT  
     cerr << endl;
     long prevtime = time(NULL);
     while ( !NoMoreUpdates() )
     {  
         // ...give a signal each second
         long currtime = time(NULL);
         if ( currtime != prevtime ) 
         {
	      prevtime = currtime;
              cerr << "processed " << iNextPoint << " points   \r" << flush;
	 }

         UpdateStep();
     }
	 
     cerr << "\rprocessed " << iNextPoint << " points " << endl;

   #else
     while ( !NoMoreUpdates() ) UpdateStep();
   #endif


   //
   // update constraints to be inserted with the new order of points
   // and insert constraints
   //

   PrepareToEnd();
  
   WriteData( outfname );
   EndTriangulation();
}

void TRefCDT :: PrepareToEnd()
{
   //
   // update constraints to be inserted with the new order of points
   //

   int num_pairs = 2*nConstr;
   int * aux_array = (int *) malloc (num_pairs*sizeof(int));
   int i;

   for (i=0; i<num_pairs; i++)
      aux_array[i] = OrderInitial[constr_array[i]];
   for (i=0; i<num_pairs; i++)
      constr_array[i] = aux_array[i];
   delete(aux_array);
   delete(OrderInitial);

   //
   // insert constraints
   //
  
   Add_Constrain();
}


void TRefCDT::WriteData( const char *outfname )
{
   TDecCDT::WriteData(outfname);
}

void TRefCDT ::ConvertData(int *vNum, int *tNum, int *eNum,
                           float **vData, int **tData, int **eData)
{
   TDecCDT::ConvertData(vNum,tNum,eNum,vData,tData,eData);
}

