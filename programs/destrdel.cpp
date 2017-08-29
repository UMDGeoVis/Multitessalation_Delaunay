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
//   file   : destrdel.cpp
//   author : Christian Melchiorre
//
//   Implementation  of class TDestroyDelaunay, for the decimation of an
//   initial triangulation that has been read from an input file.
//

//
//  ATTENTIONE: array points is declared as an array of  PTPointin base 
//  class TTriangulation, in such a way that its elements can point 
//  to TPoint as well as to TVertex (sub-class of TPoint). 
//  In TDestroyDelaunay, the elements of such array always point to TVertex.
//  In spite of that, I though it was not necessary to redefine Points
//  as an array of PTVertex. The only umpleasant thing is that we have 
//  to make a cast ...(PTVertex)Points[v]...
//


#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "defs.h"
#include "error.h"
#include "geom.h"
#include "utils.h" // PAOLA for CheckIsolatedPoint
#include "tdoublelist.h"
#include "tbtree.h"
#include "ttriang.h"
#include "basedel.h"
#include "destrdel.h"

// -------------------------------------------------------------------------
//
//  int compare( PTVertex, PTVertex )
//
//  Compare two points based on the maximum error. Return -1, 0, +1
//  according to the result of the comparison.
//  This function is required for using type TBTree<PTPoint>.
//
//  Since in class TDestroyDelaunay no value is assigned to field Error
//  of a vertex (thus, it is always 0.0), the comparison is always done
//  based on vertex coordinates.
//  


int compare( PTVertex P0, PTVertex P1 )
{

   double err0 = P0->Error;
   double err1 = P1->Error;
   
   if ( err0 < err1 ) 		// T0 < T1
      return -1;
   else if ( err0 > err1 ) 	// T0 > T1
      return +1;
   else
   {
   
       //
       // NB: it is necessary to sort well two points. Value 0 must
       // be returned only if the two points are perfectly equal.
       // In case error are equal, therefore, sort the two points based
       // on their coordinates.
       //

       int sx = Geom::Sign( P0->x - P1->x );
       int sy = Geom::Sign( P0->y - P1->y );
       int sz = Geom::Sign( P0->z - P1->z );
       
       if ( sx != 0 )
          return sx;
       else if ( sy != 0 )
          return sy;
       else 
          return sz;
       
   }
 
}   
    
      
// ---------------------------------------------------------------------------------
//
//   Constructor of class TDecimDealunay.
//

TDestroyDelaunay::TDestroyDelaunay( int iK ) : TDelaunayBase(), SwapEdgeQueue(), ElimVtxTree()
{   
    #ifdef DEBUG
     DEBUG << "TDestroyDelaunay Constructor" << endl;
    #endif // DEBUG 

    //
    // upper bound on degree of removable vertices, if KDegree == 0,
    // no upper bound is set
    //
    
    KDegree = iK;
}

   

// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::ReadData( const char * )
//
//   Read initial triangulation from input file. File format is the
//   following:
//
//        n
//        v1x  v1y  v1z        ...vertex coordinates
//        .
//        .
//        .
//        vnx  vny  vnz
//        t
//        v10  v11  v12        ...vertex indices of avery triangle
//        .
//        .
//        .
//        vt0  vt1  v12
//
//   OUTPUT: put the pointers to triangulation vertices into array Points.
//

void TDestroyDelaunay::ReadData( const char *infname )
{

   #ifdef DEBUG
    DEBUG << "\nTDestroyDelaunay::ReadData()" << endl;
   #endif
    
   ifstream inFile;
   
   inFile.open( infname );
   
   check( !inFile, "TBuildDelaunay::ReadData(), cannot open input file");

   ReadVertices( inFile );
   
   ReadTriangles( inFile );

}


// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::ReadVertices( ifstream &inFile )
//
//   Called from TDestroyDelaunay::ReadData(), read vertices from
//   input file.
//


void TDestroyDelaunay::ReadVertices( ifstream &inFile )
{

   int i;

   #ifdef DEBUG
       DEBUG << "\nTDestroyDelaunay::ReadVertices()" << endl;
   #endif
   
   
   //
   // read vertex number
   //
   
   inFile >> nPts;
   
   check( (nPts < 3), "input with less than 3 points" );

   //
   // read vertex coordinates, pointers to such vertices are put 
   // in array Points
   //
   
   Points = new PTPoint[nPts];
   
   check( (Points == NULL), "TDestroyDelaunay::ReadData(), insufficient memory" );
   
   
   PTPoint p = new TPoint;
   check( (p == NULL), "TDestroyDelaunay::ReadData(), insufficient memory");
   
   for( i=0; i<nPts; i++ )
   {
      
       check( (inFile.eof()), "TDestroyDelaunay::ReadData(), unexpected End Of File");
  
       //
       // read punto (its three coordinate) 
       //
       
       inFile >> (*p); 
       
       #ifdef DEBUG
        DEBUG << "read vertex: " << (*p) << endl;
       #endif
              
       if ( i % 1000 == 0 ) cerr << "\rread " << i << " points" << flush;
       
       Points[i] = new TVertex( p );
       check( (Points[i] == NULL), "TDestroyDelaunay::ReadData(), insufficient memory");
       
   }
   
    cerr << "\rread " << nPts << " points" << endl;
}



// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::ReadTriangles( ifstream &inFile )
//
//   Called from TDestroyDelaunay::ReadData(), read triangles from input
//   file, recover adjacency relations of the initial triangulation.
//

/*
New versione by Paola Magillo (february 2001).
For each triangle define three elements made up of edge+triangle.
Sort all these elements in lessicographic order w.r.t. to the pair
of vertex indices defining the edge.
Then, pairs of adjacent triangles are found as consecutive elements
in the sorted list, in which the edge is the same.
Some auxiliary functions are needed, which are defined below as 
"static" functions (visibility limited to this file).
*/

/* Auxiliary type for adjacency reconstruction. It represents one edge
   plus one of its adjacent triangles. */
typedef struct
{
  int ind1, ind2; /* vertex indexes of edge endpoints */
  PTTriangle t;   /* one triangle adjacent to edge */
  PTEdge e;       /* the edge (it will be created later) */
} aux;
   
/* Function that compares two such structures according to a lexicographic
   order of the two vertex indexes */
static int cmp_aux( const void *p, const void *q )
{
  aux *a = (aux *) p ;
  aux *b = (aux *) q ;
  int min_a, min_b, max_a, max_b;
  if ( a->ind1 < a->ind2 ) {  min_a = a->ind1; max_a = a->ind2;  }
  else {  min_a = a->ind2; max_a = a->ind1;  }
  if ( b->ind1 < b->ind2 ) {  min_b = b->ind1; max_b = b->ind2;  }
  else {  min_b = b->ind2; max_b = b->ind1;  }
  if ( min_a > min_b ) return 1;
  if ( min_a < min_b ) return -1;
  if ( max_a > max_b ) return 1;
  if ( max_a < max_b ) return -1;
  return 0;
}

/* Function that creates a point located at the middle point of an edge. */
static PTPoint MiddlePoint ( PTEdge e )
{
   PTPoint p;
   double xm = ( (e->EV[0]->x + e->EV[1]->x) / 2 );
   double ym = ( (e->EV[0]->y + e->EV[1]->y) / 2 );
   
   p = new TPoint (xm,ym,0.0);
   return p;
}         

/* Fix the orientation of a triangle: it must be counterclockwise.
   Use the middle points of edges to calculate orientation. */
static void OrientTriangle ( PTTriangle t )
{
  PTPoint p[3];
  PTEdge e;
  int i;
  for (i=0; i<3; i++)
    p[i] = MiddlePoint (t->TE[i]);
  if ( Geom::Turnxy( p[0], p[1], p[2] ) != TURN_LEFT )
  {  // Change orientation by swapping two edges in TE
     e = t->TE[0];
     t->TE[0] = t->TE[1];
     t->TE[1] = e;
  }
  for (i=0; i<3; i++)
    delete p[i];
}

void TDestroyDelaunay::ReadTriangles( ifstream &inFile )
{

   int e, t, v, i, k;
   aux *triedge_vec ; /* auxiliary vector for adjacency reconstuction */
   int triedge_ind = 0;

   PTTriangle NewT;
   PTEdge NewE, TmpE;
   int vidx[3];
   int nt = 0;
   PTVertex TmpV[3];

   #ifdef DEBUG
       DEBUG << "\nTDestroyDelaunay::ReadTriangles()" << endl;
   #endif 
  
   //
   // read triangles: for each triangle read three vertex indices
   //
   
   inFile >> nt;
   
   #ifdef DEBUG
      DEBUG << nt << " input triangles" << endl;
   #endif

   check( (nt < 1), "TDestroyDelaunay::ReadTriangles(), input with less than 1 triangle" );

   triedge_vec = (aux*) calloc( nt*3 , sizeof(aux) ) ;
   check( (!triedge_vec), "TDestroyDelaunay::ReadTriangles(), insufficient memory");

   for ( t=0; t<nt; t++ )
   {
       check( (inFile.eof()), "TDestroyDelaunay::ReadTriangles(),  unexpected end of file" );
       
       inFile >> vidx[0] >> vidx[1] >> vidx[2];
       
       if ( t % 1000 == 0 ) cerr << "read " << t << " triangles\r" << flush;
       
       //
       // vertices must be non-aligned, and sorted counterclockwise
       // 
              
       switch( (Geom::Turnxy( Points[vidx[0]], Points[vidx[1]], Points[vidx[2]] ) ) )
       {
          case ALIGNED:
	          error( "Triangle with three aligned vertex detected" );
	       
	  case TURN_RIGHT:
	       {
	          // ...swap vidx[1] e vidx[2]
		  int vtmp = vidx[1];
		  vidx[1] = vidx[2];
		  vidx[2] = vtmp;       
	       }
	       break;
	  case TURN_LEFT:
	       // ...Ok
	       break;
	       
       } // end ...switch
       
       #ifdef DEBUG
          DEBUG << "read triangle ( V" << vidx[0] << ", V" << vidx[1] << ", V" << vidx[2] <<  " )" << endl;
       #endif

       NewT = new TTriangle( NULL, NULL, NULL );
       check( (NewT==NULL), "TDestroyDelaunay::ReadTriangles(),  insufficient memory" );
       AddTriangle( NewT );

       for ( e=0; e<3; e++ )
       {
           triedge_vec[triedge_ind].t = NewT;
           triedge_vec[triedge_ind].ind1 = vidx[e];
           triedge_vec[triedge_ind].ind2 = vidx[(e+1)%3];
           triedge_ind++;
       }

   } // end ...for(t) (reading triangles)

   cerr << "read " << nt << " triangles" << endl;

   // Now sort all (edge,triangle) pairs. Pairs referring to the two adjacent 
   // triangles of the same edge will be consecutive in the sorted array.
   qsort(triedge_vec,3*nt,sizeof(aux),cmp_aux) ;
   
   // First scan of edge-triangle vector. I find all pairs of
   // adjacent triangles and I create all edges.
   k = 0;
   while ( k < 3*nt )
   {  if ( (k<(3*nt-1)) && ( triedge_vec[k].ind1 == triedge_vec[k+1].ind2 ) &&
           ( triedge_vec[k].ind2 == triedge_vec[k+1].ind1 ) )
      {
         // Same edge: the two triangles of the pair have a common edge,
         // create edge and set TE, ET, EV.
         NewE = new TEdge( (PTVertex)(Points[triedge_vec[k].ind1]),
                           (PTVertex)(Points[triedge_vec[k].ind2]) );
         for (t=0; t<2; t++)
         {
           for (e=0; e<3; e++)
           { // This edge is put in the 1st free place, edges will be
             // sorted counterclockwise when the last one is set.
             if (!triedge_vec[k+t].t->TE[e]) 
             {  triedge_vec[k+t].t->TE[e] = NewE; 
                if (e==2) // now all three edges are defined
                {   OrientTriangle(triedge_vec[k+t].t);
                    triedge_vec[k+t].t->CalcCircle();
                }
                break;
             }
           }
           check( (e==3), "TDestroyDelaunay::ReadTriangles(), inconsistency");
           NewE->ET[t] = triedge_vec[k+t].t;
           triedge_vec[k+t].e = NewE;
         }
         k+= 2; // go to next pair 
      }
      else
      {
         // Different edge: the first triangle of the pair is a boundary
         // triangle, create edge and set TE, ET, EV.
         NewE = new TEdge( (PTVertex)(Points[triedge_vec[k].ind1]),
                           (PTVertex)(Points[triedge_vec[k].ind2]) ); 
         for (e=0; e<3; e++)
         { // This edge is put in the 1st free place, edges will be 
           // sorted counterclockwise when the last one is set.
           if (!triedge_vec[k].t->TE[e])
           {  triedge_vec[k].t->TE[e] = NewE; 
              if (e==2) // now all three edges are defined
              {   OrientTriangle(triedge_vec[k].t);
                  triedge_vec[k].t->CalcCircle();
              }
              break;
           }
         }
         check( (e==3), "TDestroyDelaunay::ReadTriangles(), inconsistency");
         NewE->ET[0] = triedge_vec[k].t;
         triedge_vec[k].e = NewE;
         k++; // go to next element
      }

      // Set VE for the two endpoints of NewE (if still null).
      for (v=0; v<2; v++)
      {
        if ( NewE->EV[v]->VE[0] == NULL )
	   NewE->EV[v]->VE[0] = NewE;
        else if ( NewE->EV[v]->VE[1] == NULL ) 
           NewE->EV[v]->VE[1] = NewE;
      }
      CheckIsolatedPoint(NewE->EV[0]);
      CheckIsolatedPoint(NewE->EV[1]);

   } // while k

   /* Second scan of edge-triangle vector. Now I fix the 
      relation VE for boundary vertices */
   k = 0;
   while ( k < (3*nt) )
   {  if ( (k<(3*nt-1)) && ( triedge_vec[k].ind1 == triedge_vec[k+1].ind2 ) &&
           ( triedge_vec[k].ind2 == triedge_vec[k+1].ind1 ) )
      {
         // Same edge: it is not a boundary edge.
         k+= 2; /* go to next pair */
      }
      else
      {
         // Different edge: the edge in the first element is a
         // boundary edge.
         // Fix VE for vertices.
         TmpE = triedge_vec[k].e;
         for (v=0; v<2; v++)
         {  
           if ( (TmpE->EV[v]->VE[0]!=TmpE) && (TmpE->EV[v]->VE[1]!=TmpE) )
           {
             if (!TmpE->EV[v]->VE[0]->OnConvexHull())
                TmpE->EV[v]->VE[0] = TmpE;
             else
                TmpE->EV[v]->VE[1] = TmpE;
           }
         } // for v
         k++; /* go to next element */
      }
   } // while k

/* I DON'T THINK THIS IS NECCESSARY, EVEN IF IT WORKS
   // Ensure that, for all boundary vertices, going from VE[0] to VE[1] 
   // is counterclockwise
   
   for( i=0; i<nPts; i++ )
   {
      PTVertex vv = (PTVertex)(Points[i]);   

      if (vv->VE[0]->OnConvexHull())
      {
        check (!(vv->VE[1]->OnConvexHull()),
               "TDestroyDelaunay::ReadTriangles(), inconsistency detected in VE");
        if ( vv == vv->VE[0]->EV[0] ) k = 0; else k = 1;
        // if k == 0 null triangle of VE[0] must be on its left side
        // if k == 1 null triangle of VE[0] must be on its right side
        if ( vv->VE[0]->ET[k] != NULL )
        {
          TmpE = vv->VE[0];
          vv->VE[0] = vv->VE[1];
          vv->VE[1] = TmpE;
        }
      }
   }
END OF UNNECESSARY THING */

    #ifdef DEBUG
      cerr << "VE of vertices after adjusting: " << endl;
      for( i=0; i<nPts; i++ )
      {
	PTVertex vv = (PTVertex)(Points[i]);   
	PTEdge eee;
//        cerr << "V" << vv->VID << ":  E" << vv->VE[0]->EID << flush 
//	                        << ", E" << vv->VE[1]->EID << endl;
        vv->GetVE();
        cerr << "Edges of V" << vv->VID << "= " ;
        while (!(vv->EdgeList.IsEmpty()))
        {
          eee = vv->EdgeList.RemoveHead();
          cerr << " E" << eee->EID << "(V" << eee->EV[0]->VID << ",V"
               << eee->EV[1]->VID ;
        }
        cerr << endl;
      }
      cerr << endl;
   #endif

}

// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::InitialTriangulation()
//
//   Since function ReadData() reads the initial triangulation from input
//   file, function InitialTriangulation is not needed.
//   But this is the most suitable place to find out which vertices
//   of the  initial triangulationare removable.
//   A vertex is removable if either it not on the convex hull,
//   or the two boundary edges incident in it are aligned.
//

void TDestroyDelaunay::InitialTriangulation()
{

   int v;

   #ifdef DEBUG
       DEBUG << "\nTDestroyDelaunay::InitialTriangulation()" << endl;
   #endif
   
   
   for( v=0; v<nPts; v++ )
   {
     if( ((PTVertex) Points[v])->VE[0] != NULL ) 
     // in order to read also files that contain extra points, not
     // used as vertices of the triangulation, for such vertices only
     // coordinates are initialized
     // PAOLA: this trick should no longer be necessary after the new 
     // version of writing function!
     {
       //
       // A vertex is removable if either it not on the convex hull,
       // or the two boundary edges incident in it are aligned
       // (i.e., the vertex and the other two endpoints of such edges
       // are aligned). Removable vertices that satisfy the upper bound
       // on the degree (if present) are put in ElimVtxTree.
       //
       
       PTVertex V = (PTVertex)(Points[v]);
      
       if ( IsVtxElim( V ) && ( OkDegree( CalcDegree(V) ) ) )
           ElimVtxTree.Insert( V );
       
     } // end ...if(  Points[v]->VE[0] != NULL
   } // end ...for(v)

   #ifdef DEBUG
      DEBUG << "Removable vertices : " << endl;
      ElimVtxTree.Print();
   #endif

   
}



// -----------------------------------------------------------------------------
//  
//   boolean TDestroyDelaunay::OkDegree( int deg )
//
//   Return true if there is no upper bound on the degree of removable
//   vertices (KDegree == 0) o if deg is within the bound ( deg <= KDegree )
//

boolean TDestroyDelaunay::OkDegree( int deg )
{
   return( KDegree == 0 || deg <= KDegree );
}


// -----------------------------------------------------------------------------
//  
//   int TDestroyDelaunay::CalcDegree( PTVertex )
//
//   Compute the number of edges incident in a vertex V
//

int TDestroyDelaunay::CalcDegree( PTVertex V )
{
   int deg = 0;
   PTEdge   TmpE;

   #ifdef DEBUG
     DEBUG << "TDestroyDelaunay::CalcDegree( V" << V->VID << " )" ;
   #endif
   
   V->GetVE();
   while (!(V->EdgeList.IsEmpty()))
   {
      TmpE = V->EdgeList.RemoveHead();
      deg++;
   }
   V->EdgeList.ClearList();
   return( deg );
}


// -----------------------------------------------------------------------------
//  
//   boolean TDestroyDelaunay::IsVtxElim( PTVertex )
//
//   Return TRUE iff the given vertex is removable (as explained before).
//

boolean TDestroyDelaunay::IsVtxElim( PTVertex V )
{

    boolean elim = FALSE;
       
    if ( ! V->VE[0]->OnConvexHull() )
    {
          #ifdef ROBUST
	      // either VE[0] and VE[1] both on the convex hull,
              // or none of them is
	      check( (V->VE[1]->OnConvexHull()),
	             "TDestroyDelaunay::InitialTriangulation(), inconsistency detected" );
	  #endif
	  
          elim = TRUE;
    }
    else // VE[0]->OnConvexHull()
    {
          #ifdef ROBUST
	      // either VE[0] and VE[1] both on the convex hull,
              // or none of them is
	      check( (!V->VE[1]->OnConvexHull()),
	             "TDestroyDelaunay::InitialTriangulation(), inconsistency detected" );
	  #endif
	  
	  
	  // va and vb are the two vertices on the convex hull which are
          // "neighbors" of Points[v]
	  
	  PTVertex va = ( V->VE[0]->EV[0] != V ? 
	                  V->VE[0]->EV[0] :
			  V->VE[0]->EV[1] );
          
	  PTVertex vb = ( V->VE[1]->EV[0] != V ? 
	                  V->VE[1]->EV[0] :
			  V->VE[1]->EV[1] );
			 
          elim = ( Geom::Alignedxy( va, V, vb ) );
	  	  
    }
    
    return( elim );

}


// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::UpdateStep()
//
//   Perform one step of the update process on the current triangulation,
//   by removing a vertex from it.
//   First, find the next vertex to be removed (the min of tree ElimVtxTree)
//   by calling NextVertex().
//   Then, remove it by calling RemoveVertex().
//

void TDestroyDelaunay::UpdateStep()
{
   #ifdef DEBUG
      DEBUG << "\nTDestroyDelaunay::UpdateStep()" << endl;
     DEBUG << "ntrg = " << nTrg << endl;
   #endif
   
   //
   // TDestroyDelaunay::NextVertex(), put the pointer to the next vertex
   // to be removed in VertexToRemove.
   // 
   
   NextVertex();
   
   //
   // Now VertexToRemove is removed from the triangulation
   //
   
   RemoveVertex();
  
}


// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::NextVertex()
//
//   Find the next vertex to be removed from the current triangulation.
//   The pointer to such vertex is put in the status variable
//   VertexToRemove.
//

void TDestroyDelaunay::NextVertex()
{

    #ifdef DEBUG
       DEBUG << "\nTDestroyDelaunay::NextVertex()" << endl;
    #endif

    #ifdef ROBUST
       check( (ElimVtxTree.IsEmpty()), "TDestroyDelaunay::NextVertex(), no more points to remove" );
    #endif
 
    VertexToRemove = ElimVtxTree.RemoveMin();

    #ifdef DEBUG    
       DEBUG << "VertexToRemove: V" << VertexToRemove->VID  << endl;
    #endif
             
}


// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::RemoveVertex()
//
//   Remove the vertex pointed by the status variable VertexToRemove from
//   the triangulation. This function is symmetric w.r.t. InsertVertex in
//   sub-classes of TBuildDelaunay.
//
//   INPUT: VertexToRemove
//
//   OUTPUT: a consistent status of TDestroyDelaunay, representing the new
//   triangulation obtained after removing VertexToRemove.
//

void TDestroyDelaunay::RemoveVertex()
{
     boolean closed = CalcInfluenceRegion();
     RetriangulateInfluenceRegion( closed );
}


// -----------------------------------------------------------------------------
//  
//   boolean TDestroyDelaunay::CalcInfluenceRegion()
//
//   Compute the influence region of the vertex to be removed, i.e., the
//   list of boundary edges of the set of triangles adjacent to
//   VertexToRemove. Such edges are sorted counterclockwise around
//   VertexToRemove.
//
//   INPUT: VertexToRemove
//
//   OUTPUT: fill the doubly linked list InflRegnBorder, as explained above.
//   Pointer to TTriangle FirstTrgToDel is set to point one of the triangles
//   of the influence region by one of the functions InitInflRegnXXX().
//
//   RETURN VALUE: TRUE if InflRegnBorder is a closed sequence of edges
//   (i.e., if VertexToRemove is inside the triangulation domain),
//   FALSE otherwise (i.e., if VertexToRemove is a removable vertex 
//   on the convex hull).
//

boolean TDestroyDelaunay::CalcInfluenceRegion()
{
   
    #ifdef DEBUG
       DEBUG << "\nTDestroyDelaunay::CalcInfluenceRegion()" << endl;
    #endif
    
    #ifdef ROBUST
       check( (!InflRegnBorder.IsEmpty()),
            "TDestroyDelaunay::CalcInfluenceRegion(), <a> inconsistency detected" );
    #endif
    
    boolean closed;

    PTEdge NextEdge = VertexToRemove->VE[0]; // ...one ot the two...
    
    #ifdef ROBUST
       check( (NextEdge == NULL),
          "TDestroyDelaunay::CalcInfluenceRegion(), <b> inconsistency detected" );
    #endif
    
    
    if ( VertexToRemove->VE[0]->OnConvexHull() ) 
    // ...removable vertex on convex hull
    {      
       #ifdef DEBUG
          DEBUG << "removing V" << VertexToRemove->VID << endl;
       #endif
       
       #ifdef ROBUST
          
	  // some consistency check
	  
          check( (!VertexToRemove->VE[1]->OnConvexHull()),
                "TDestroyDelaunay::CalcInfluenceRegion(), <2> inconsistency detected" );
		
          PTVertex V0 = ( VertexToRemove->VE[0]->EV[0] != VertexToRemove ?
	                  VertexToRemove->VE[0]->EV[0] :
			  VertexToRemove->VE[0]->EV[1] );
          PTVertex V1 = ( VertexToRemove->VE[1]->EV[0] != VertexToRemove ?
	                  VertexToRemove->VE[1]->EV[0] :
			  VertexToRemove->VE[1]->EV[1] );
			  
	  check( (!Geom::Alignedxy( V0, VertexToRemove, V1)),
                "TDestroyDelaunay::CalcInfluenceRegion(), <3> inconsistency detected" );
			  
       #endif // ROBUST
       
       InitInflRegnHull();      
       closed = FALSE; 

    }
    else // ...interior vertex
    {
       #ifdef ROBUST
          check( (VertexToRemove->VE[1]->OnConvexHull()),
                "TDestroyDelaunay::CalcInfluenceRegion(), <4> inconsistency detected" );
       #endif
       
       InitInflRegnInternal();  
       closed = TRUE;
       
    }
    
    //
    // Now InflRegnBorder has been initialized to contain one of the
    // boundary edges of the influence region of VertexToRemove. The 
    // following edges are inserted in such list by calling function
    // CalcInflRegnMain(). Parameter closed specifies if the boundary
    // of the influence region is closed (the vertex was internal) or
    // open (the vertex was on the convex hull). 
    //
    #ifdef DEBUG
        TDoubleListIterator<PTEdge> Iter( &InflRegnBorder );
        Iter.Restart();
	while( !Iter.EndOfList() )
	{
	   DEBUG << " ->E" << (Iter.Current()->object)->EID;
	   Iter.GoNext();
	}
	DEBUG << endl;
    #endif	 
	
    CalcInflRegnMain();     

    #ifdef DEBUG
        Iter.Restart();
	while( !Iter.EndOfList() )
	{
	   DEBUG << " ->E" << (Iter.Current()->object)->EID;
	   Iter.GoNext();
	}
	DEBUG << endl;
        DEBUG << "FirstTrgToDel = ";
        if ( FirstTrgToDel == NULL ) DEBUG << "NULL"; else DEBUG << "T" << FirstTrgToDel->TID;
        DEBUG << endl;
    #endif	 
    
    return( closed );
    
}    
   

// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::InitInflRegnHull()
//
//   Insert the first edge in the sequence of boundary edges of the 
//   influence region in case VertexToRemove is on the convex hull.
//
//   INPUT: VertexToRemove
//
//   OUTPUT: List InflRegnBorder contains the first edge of the influence
//   region. Variable FirstTrgToDel is set to point to one of the triangles
//   of the influence region (it is used as an entry point to the triangles
//   of the influence region).
//

void TDestroyDelaunay::InitInflRegnHull()
{

    int i;

    #ifdef DEBUG
        DEBUG << "\nTDestroyDelaunay::InitInflRegnHull()" << endl;
    #endif

    //
    // E0 and E1 are the two edges adjacent to VertexToRemove. Since 
    // VertexToRemove is on the convex hull, the two edges must be on
    // the convex hull, too (we don't check, we have already checked
    // in CalcInfluenceRegion).
    //

    PTEdge E0 = VertexToRemove->VE[0];
    PTEdge E1 = VertexToRemove->VE[1];
    
    //
    // V0 and V1 the two endpoints of E0 and E1 (respectively) different 
    // from VertexToRemove. Since we have imposed that E0 and E1 are
    // aligned, V0 and V1 are distinct.
    //
    
    PTVertex V0 = ( E0->EV[0] != VertexToRemove ? E0->EV[0] : E0->EV[1] );
    PTVertex V1 = ( E1->EV[0] != VertexToRemove ? E1->EV[0] : E1->EV[1] );    

    //
    // T0 are T1 the two triangles having one edge on the convex hull
    // and equal to E0 and E1 (respectively). Since we have imposed that E0
    // and E1 are aligned, T0 and T1 are distinct.
    //
        
    PTTriangle T0 = ( E0->ET[0] != NULL ? E0->ET[0] : E0->ET[1] );
    PTTriangle T1 = ( E1->ET[0] != NULL ? E1->ET[0] : E1->ET[1] ); 
    
    #ifdef ROBUST
       check( (T0 == NULL || T1 == NULL), "TDestroyDelaunay::InitInflRegnHull(), <0> inconsistency detected");
    #endif
    
    FirstTrgToDel = T0;

    //
    // V0b and V1b are the vertices of T0 and T1, different from
    // VertexToRemove/V0 and VertexToRemove/V1 (respectively).
    // It may happen that V0b == V1b. E0b and E1 are the adjacent edges
    // to V0b/V1b different from E0/E1 (i.e., the edges of T0/T1 which 
    // are NOT adjacent to VertexToRemove).
    //
    
    PTVertex tv[3];
    
    // compute V0b/E0b

    PTVertex V0b = NULL;
    PTEdge E0b = NULL;
    
    T0->GetTV( tv[0], tv[1], tv[2] );
        
    for( i=0; i<3; i++ ) 
    {
       if ( tv[i] != VertexToRemove && tv[i] != V0 )
          V0b = tv[i];
	  
       if ( tv[(i+2)%3] == VertexToRemove )
          E0b = T0->TE[i]; // e' lo spigolo tv[i]..tv[(i+1)%3]
    }
    
    #ifdef ROBUST
       check( (V0b==NULL || E0b==NULL), "TDestroyDelaunay::InitInflRegnHull(), <1> inconsistency detected" );
    #endif
    
    // compute V1b/E1b 

    PTVertex V1b = NULL;
    PTEdge E1b = NULL;
    
    T1->GetTV( tv[0], tv[1], tv[2] );
        
    for( i=0; i<3; i++ ) 
    {
       if ( tv[i] != VertexToRemove && tv[i] != V1 )
          V1b = tv[i];
	  
       if ( tv[(i+2)%3] == VertexToRemove )
          E1b = T1->TE[i]; // it is edge tv[i]..tv[(i+1)%3]
    }
    
    #ifdef ROBUST
       check( (V1b==NULL || E1b==NULL), "TDestroyDelaunay::InitInflRegnHull(), <2> inconsistency detected" );
    #endif

    if ( Geom::Turnxy( VertexToRemove, V0, V0b ) == TURN_LEFT )
    {
       FirstTrgToDel = T0;
       InflRegnBorder.AddHead( E0b );  E0b->Mark( INFL_BORDER );
    }
    else if ( Geom::Turnxy( VertexToRemove, V1, V1b ) == TURN_LEFT )
    {
       FirstTrgToDel = T1;
       InflRegnBorder.AddHead( E1b );  E1b->Mark( INFL_BORDER );
    } 
    else
    {
      error( "TDestroyDelaunay::InitInflRegnHull(), <3> inconsistency detected" );
    }

}


// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::InitInflRegnInternal()
//
//   Insert the first
//   Insert the first edge in the sequence of boundary edges of the
//   influence region in case VertexToRemove is an internal vertex to
//   the domain of the triangulation.
//
//   INPUT: VertexToRemove
//
//   OUTPUT: List InflRegnBorder contains one (arbitrary) edge of the
//   influence region. Such edge is marked as INFL_BORDER.
//

void TDestroyDelaunay::InitInflRegnInternal()
{

    int i;

    #ifdef DEBUG
        DEBUG << "TDestroyDelaunay::InitInflRegnInternal()" << endl;
    #endif
    
    #ifdef ROBUST
        check( (VertexToRemove->VE[0] == NULL),
	    "TDestroyDelaunay::InitInflRegnInternal(), <1> inconsistency detected" );
	
	// edge must be internal (not on convex hull)
	
        check( (VertexToRemove->VE[0]->OnConvexHull()),
	    "TDestroyDelaunay::InitInflRegnInternal(), <2> inconsistency detected" );
    #endif // ROBUST
    
    
    PTTriangle TFirst = VertexToRemove->VE[0]->ET[0];
    
    #ifdef ROBUST
       check( (TFirst == NULL), "TDestroyDelaunay::InitInflRegnInternal(), <3> inconsistency detected" );
    #endif
    
    FirstTrgToDel = TFirst;
    
    //
    // EFirst is the only edge of TFirst not adjacent to VertexToRemove
    //
    
    for ( i=0; i<3; i++ )
       if ( TFirst->TE[i]->EV[0] != VertexToRemove && TFirst->TE[i]->EV[1] != VertexToRemove )
          break;
	  
    #ifdef ROBUST
        check( (i>=3), "TDestroyDelaunay::InitInflRegnInternal(), <4> inconsistency detected" );
    #endif

    PTEdge EFirst = TFirst->TE[i];
    InflRegnBorder.AddHead( EFirst );
    EFirst->Mark( INFL_BORDER );    

}


// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::CalcInflRegnMain()
//
//   Complete the determinatione of the influence region of vertex
//   VertexToRemove.
//
//   INPUT: VertexToRemove, plus the content of list InflRegnBorder,
//   this has been initialized by InitInflRegnXXX() with one (if
//   internal vertex) of three (if convex hull) edges.
//
//   OUTPUT: The above list InflRegnBorder, that will contain the set
//   of edges belonging to the boundary of the influence region.
//   Edges in InflRegnBorder will be marked as INFL_BORDER, while 
//   edges and triangles lying inside the influence region will be marked
//   as TO_DELETE.
//

void TDestroyDelaunay::CalcInflRegnMain()
{

   int i, e;

   #ifdef DEBUG
       DEBUG << "TDestroyDelaunay::CalcInflRegnMain()" << endl;
   #endif
   
   PTEdge EFirst = InflRegnBorder.GetHead();
   PTEdge E = EFirst;

   PTTriangle T;

   for( i=0; i<2; i++ ) // i=0/1
   {
      // T is the triangolo adjacent to EFirst and to VertexToRemove
      
      T = E->ET[i];
      
      if ( T != NULL )
      {
         PTVertex v0, v1, v2;
         T->GetTV( v0, v1, v2 );

         if ( v0 == VertexToRemove || v1==VertexToRemove || v2==VertexToRemove ) break;
      }
   }
   
   #ifdef ROBUST
      check( (i>=2), "TDestroyDelaunay::CalcInflRegnMain(), <1> inconsistency detected" );
   #endif
      
   PTEdge ENext;
   PTTriangle TNext;   

   //
   // main loop: turn around VertexToRemove in counterclockwise order
   //
   
   while( ! T->Marked( TO_DELETE ) && T != NULL )
   {

      #ifdef DEBUG 
         DEBUG << "esamino T" << T->TID << endl;
      #endif
            
      // mark T and all its adjacent edges different from E as TO_DELETE
      
      T->Mark( TO_DELETE );
      for ( e=0; e<3; e++ )
        if ( T->TE[e] != E )
	   T->TE[e]->Mark( TO_DELETE );

      // next triangle of influence region
      
      for( e=0; e<3; e++ ) if ( T->TE[e] == E ) break;
   
      #ifdef ROBUST
          check( (e>=3), "TDestroyDelaunay::CalcInflRegnMain(), <1> inconsistency detected" );
      #endif
            
      ENext = T->TE[(e+1)%3];
      TNext = ( ENext->ET[0] != T ? ENext->ET[0] : ENext->ET[1] );
            
      if ( TNext == NULL || TNext->Marked( TO_DELETE ) ) break;
           // exit while
      
      // next edge on boundary of influence region
      
      for ( e=0; e<3; e++ ) if ( TNext->TE[e] == ENext ) break;

      #ifdef ROBUST
         check( (e>=3), "TDestroyDelaunay::CalcInflRegnMain(), <2> inconsistency detected" );
      #endif
      
      E = TNext->TE[(e+1)%3];
      InflRegnBorder.AddTail( E ); 
      // ...enqueue so that edges are inserted in counterclockwise order
      E->Mark( INFL_BORDER );
      T = TNext;
      
   } // end ...while( !T->Marked( TO_DELETE ) )
    
}


// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::RetriangulateInfluenceRegion( boolean closed )
//
//   Retriangulate the region of influence whose boundary is formed by
//   the edges in the doubly linked list InflRegnBorder.
//
//   INPUT: list InflRegnBorder and vertex VertexToRemove. Parameter
//   closed (boolean) specifies if InflRegnBorder is a closed sequence
//   of edges (VertexToRemove is inside the domain) or not (VertexToRemove
//   is a removable vertex on the convex hull). e' un vertice
//   Edges on the boundary of the influence region are marked as
//   INFL_BORDER, and edges and triangles inside the influence region are
//   marked as TO_DELETE.
//
//   Note: list InflRegnBorder represents a sequence of edges bounding the
//   regione of influence of vertex VertexToRemove, to be removed from
//   the triangulation, sorted in counterclockwise order around it.
//   It may be a closed or open sequence (as explained before), and this
//   is encoded in the boolean parameter closed. In case it is open, the
//   two edges of the convex hull incident in VertexToRemove are marked as
//   TO_DELETE.
//
//   OUTPUT: a new consistent status of TDestroyDelaunay representing the
//   triangulation after the removal of VertexToRemove. List InflRegnBorder
//   is emptied and all edges contained in it loose the INFL_BORDER mark.
//   Because of call to DeleteInfluenceRegion(), all edges and triangles
//   marked as TO_DELETE are deleted.
//

void TDestroyDelaunay::RetriangulateInfluenceRegion( boolean closed )
{

   #ifdef DEBUG
       DEBUG << "TDestroyDelaunay::RetriangulateInfluenceRegion()" << endl;
  #endif
   
   //
   // Step 1: If InflRegnBorder is open, create a new edge connecting the
   // first and last vertex of the boundary of the influence region.
   // Relations ET of such vertex are NULL by now. Do not use them before
   // they have been set!
   //
   
   PTEdge NewEdge;
   
   //
   // find VFirst and VLast, first and last vertex of the sequence of
   // edges in InflRegnBorder.
   //
   
   PTEdge EFirst = InflRegnBorder.GetHead();
   PTEdge ELast  = InflRegnBorder.GetLast();
   
   PTVertex VFirst = ( (Geom::Turnxy( VertexToRemove, EFirst->EV[0], EFirst->EV[1] ) == TURN_LEFT ) ?
                       EFirst->EV[0] :
		       EFirst->EV[1] );
		   
   PTVertex VLast  = ( (Geom::Turnxy( VertexToRemove, ELast->EV[0], ELast->EV[1] ) == TURN_RIGHT ) ?
                       ELast->EV[0] :
	               ELast->EV[1] );
		       
   #ifdef ROBUST
      check (  ( (closed && VFirst!=VLast) || (!closed && VFirst==VLast) ),
               "TDestroyDelaunay::RetriangulateInfluenceRegion(), <1> inconsistency detected" );
   #endif // ROBUST   

   if ( !closed )
   {
       // create new edge and put it in list
       
       NewEdge = new TEdge( VFirst, VLast );
       check( (NewEdge == NULL), "TDestroyDelaunay::RetriangulateInfluenceRegion(), insufficient memory" );
       
       // temporarily adjust VE relations of VFirst and VLast (now 
       // they point to edges of the convex hull to be deleted
       
       if ( VFirst->VE[0]->Marked( TO_DELETE ) ) VFirst->VE[0] = NewEdge;
       else
       if ( VFirst->VE[1]->Marked( TO_DELETE ) ) VFirst->VE[1] = NewEdge;
       else
       error( "TDestroyDelaunay::RetriangulateInfluenceRegion(),  <a> inconsistency detected" );
       
       if ( VLast->VE[0]->Marked( TO_DELETE ) ) VLast->VE[0] = NewEdge;
       else
       if ( VLast->VE[1]->Marked( TO_DELETE ) ) VLast->VE[1] = NewEdge;
       else
       error( "TDestroyDelaunay::RetriangulateInfluenceRegion(),  <b> inconsistency detected" );
       
       CheckIsolatedPoint(VFirst);
       CheckIsolatedPoint(VLast);
       
       #ifdef DEBUG
          DEBUG << "Create edge E" << NewEdge->EID 
	        << "( V" << VFirst->VID << ", V" << VLast->VID << " )" << endl;
       #endif
       
       InflRegnBorder.AddTail( NewEdge );
       NewEdge->Mark( INFL_BORDER );
       
       #ifdef ROBUST
         // now we can set closed = TRUE and VLast = VFirst (probably
         // we will not use such variables any more, but for consistency
         // it is better to set them).
	 closed = TRUE;
	 VLast  = VFirst;
       #endif // ROBUST       
       
   }
   
   //
   // Step 2: Adjust VE relations for all endpoint vertices of edges in
   // InflRegnBorder. Set such relations in such a way that every vertex
   // always points to the two edges in InflRegnBorder adjacent to it.
   //

   TDoubleListIterator<PTEdge> IRegnIter( &InflRegnBorder );
   IRegnIter.Restart();
   
   PTEdge ECurr;
   PTEdge EPrev = InflRegnBorder.GetLast();
   
   PTVertex VCurr = NULL;
   
   while ( ! IRegnIter.EndOfList() )
   {
      ECurr = IRegnIter.Current()->object;
      
      // VCurr is the common vertex of EPrev and ECurr, i.e., the second
      // vertex, in counterclockwise order w.r.t. VertexToRemove, of EPrev
      
      if ( EPrev->EV[0] == ECurr->EV[0] || EPrev->EV[0] == ECurr->EV[1] )
        VCurr = EPrev->EV[0];
      else
        VCurr = EPrev->EV[1];
      
      #ifdef ROBUST
         check( (ECurr->EV[0]!=VCurr && ECurr->EV[1]!=VCurr),
	      "TDestroyDelaunay::RetriangulateInfluenceRegion(), <2> inconsistency detected" ); 
      #endif
      
      // adjacency relations VE of VCurr are set to point ECurr/EPrev, if
      // they were pointing to internal edges of the influence region (to 
      // be deleted).

      boolean E0Del = VCurr->VE[0]->Marked( TO_DELETE );
      boolean E1Del = VCurr->VE[1]->Marked( TO_DELETE );
      
      /* #ifdef DEBUG
         DEBUG << "EPrev : E" << EPrev->EID << "( V"
	       << EPrev->EV[0]->VID << ", V" << EPrev->EV[1]->VID << " )" << endl;
         DEBUG << "ECurr : E" << ECurr->EID << "( V"
	       << ECurr->EV[0]->VID << ", V" << ECurr->EV[1]->VID << " )" << endl;
	 DEBUG << "VCurr : V" << VCurr->VID 
	       << "( E" << VCurr->VE[0]->EID << ", E" << VCurr->VE[1]->EID << " )"
	       << endl;
	 DEBUG << "E0Del : " << E0Del << endl;
	 DEBUG << "E1Del : " << E1Del << endl;
      #endif */

      if ( E0Del && E1Del )
      {
          VCurr->VE[0] = EPrev;
          VCurr->VE[1] = ECurr;      
      }
      else if ( !E0Del && E1Del )
      {
          VCurr->VE[1] = ( VCurr->VE[0] != EPrev ? EPrev : ECurr );
      }
      else if ( E0Del && !E1Del )
      {
          VCurr->VE[0] = ( VCurr->VE[1] != EPrev ? EPrev : ECurr );
      }
    
      CheckIsolatedPoint(VCurr);

      // else ... !E0Del && !E1Del... leave VE as they are

     
      // next edge
      
      EPrev = ECurr;
      IRegnIter.GoNext();
      
      
   } // end ... while ( ! EndOfList() )
  
    
   //
   // Main steps of RetriangulateInfluenceRegion():
   // RIRInitialTrg() crete an initial triangulation of the influence region
   // (trying to join a vertex with all the others),
   // then RIRDelOptTrg() optimize such triangulation according to the
   // Delaunay criterion.
   //
     
   {
       RIRInitialTrg( InflRegnBorder, SwapEdgeQueue );
       RIRDelOptTrg();
   }


   //
   // last operations to empty InflRegnBorder and remove INFL_BORDER marks
   // from all edges in it
   //
   
   // in order that DetachEdge in DeleteInfluenceRegion gives no error
   VertexToRemove->VE[0] = VertexToRemove->VE[1] = NULL;

   DeleteInfluenceRegion();   

   //
   // remove INFL_BORDER marks from edges in InflRegnBorder
   //
   
   IRegnIter.Restart();
   while ( ! IRegnIter.EndOfList() )
   {
       IRegnIter.Current()->object->UnMark( INFL_BORDER );
       IRegnIter.GoNext();
   }
      
   // empty InflRegnBorder, it is no longer needed
   
   InflRegnBorder.ClearList();
   
}


// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::RIRInitialTrg( TDoubleList<PTEdge> & InflRegnBorder,
//                                         TDoubleList<PTEdge> & SwapEdgeQueue   );
//
//   The name stands for (R)etriangulate (I)nfluence (R)egion (Initial)
//   (Tr)iangulation. One of the
//   auxiliary functions used by RetriangulateInfluenceRegion
//   (the names of all such functions start with RIR).
//
//   This function creates an initial triangulation of the star-shaped
//   polygon whose boundary is formed by the edges in InflRegnBorder.
//   If possible, join one of its vertices with all the others.
//   Retriangulate the polygon by using the 'ears' algorithm.
//   Scan the boundary edges in counterclockwise order and, for each pair
//   ( E1, E2 ) of adjacent edges, if they define a left turn and the 
//   potential new edge E (making a new triangle with E1 and E2) intersects
//   no boundary edge, then add the new triangle and replace E1 and E2
//   with E in list InflRegnBorder.
//
//   INPUT: InflRegnBorder and VertexToRemove.
//   Edges of the boundary of the influence region are marked as
//   INFL_BORDER, and edges and triangles lying inside the influence region 
//   are marked as TO_DELETE.
//
//   Note (1): InflRegnBorder, after step 1 in RetriangulateInfluenceRegion,
//   is always a closed chain of edges.
//
//   Note (2): One of the edges in InflRegnBorder may be adjacency
//   relations ET not yet set (both equal to NULL, if it has been created
//   to close an open InflRegnBorder). Therefore, we do not rely on the
//   consistency of such relations, while we will set them inside this
//   function.
//
//   OUTPUT: a consistent status corresponding to a triangulation
//   obtained after retriangulating the influence region of the removed
//   vertex VertexToRemove. Such triangulation may not be a Delaunay one.
//   Variable FirstTriangle is modified by the calls to AddTriangle,
//   in such a way to point the to last created triangle.
//   Newly created triangles are marked as NEW_TRIANGLE.
//
//   Pointers to edges created inside this function are put in a doubly
//   linked list (SwapEdgeQueue) that will be used fro optimizing the
//   triangulation (function RIRDelOptTrg); such edges are also
//   marked as SWAP_EDGE_QUEUE.
//

void TDestroyDelaunay::RIRInitialTrg( TDoubleList<PTEdge> & InflRegnBorder,
                                      TDoubleList<PTEdge> & SwapEdgeQueue   )
{

    int e;

    #ifdef DEBUG
       DEBUG << "TDestroyDelaunay::RIRInitialTrg()" << endl;
    #endif
    
    PTTriangle NewTriangle;
  
    //
    // Make a copy of list InflRegnBorder. We will modify this second
    // list, while we need to maintain the content of InflRegnBorder 
    // when we call function RIRDelOptTrg(). Similarly, we use a different
    // mark value, INFL_BORDER_AUX, to mark edges belonging to this 
    // second list.
    //
    
    TDoubleList<PTEdge> InflRegnAux;
    
    TDoubleListIterator<PTEdge> IRegnIter( &InflRegnBorder );
    
    IRegnIter.Restart();
    while( !IRegnIter.EndOfList() )
    {
       PTEdge E = IRegnIter.Current()->object;
       
       E->Mark( INFL_BORDER_AUX );
       InflRegnAux.AddTail( E );
       
       IRegnIter.GoNext();
    }
    
    //
    // start the main loop of retriangulation of the star-shaped region
    // whose boundary edges are in InflRegnBorder/InflRegnAux
    //

    TDoubleListIterator<PTEdge> IAuxIter( &InflRegnAux );
    
    IAuxIter.Restart();
    
    // current edge and next edge on InflRegnAux
    PTEdge ECurr = IAuxIter.Current()->object;
    PTEdge ENext;

    // next three vertices (in counterclockwise order) on InflRegnAux
    PTVertex v0, v1, v2;

    while( InflRegnAux.Lenght() > 3 )
    {
 
       IAuxIter.GoNext();
       if ( IAuxIter.EndOfList() ) IAuxIter.Restart(); 
       // manage list as if it were circular
 
       // next edge (the one after ECurr)
       ENext = IAuxIter.Current()->object;
       
       // v0, v1 and v2 are the endpoints of ECurr and ENext (v1 is
       // common) in counterclockwise order.
              
       v0 = ECurr->EV[0];
       v1 = ECurr->EV[1];

       if ( ENext->EV[0] == v0 || ENext->EV[1] == v0 )
       {
          // the common vertex is v0! swap v0 and v1
	  PTVertex vTmp = v0;
	  v0 = v1;
	  v1 = vTmp;
	  // now the common vertex is v1
       }

       #ifdef ROBUST
           check( (ENext->EV[0] != v1 && ENext->EV[1] != v1 ),
	      "TDestroyDelaunay::RIRInitialTrg(), <1> inconsistency detected" );
       #endif
	
       v2 = ( ENext->EV[0] != v1 ? ENext->EV[0] : ENext->EV[1] );
       
	   
       /* #ifdef DEBUG
          TDoubleListIterator<PTEdge> Iter( &InflRegnAux );
	  DEBUG << endl;
          DEBUG << "ECurr : E" << ECurr->EID << endl;
	  DEBUG << "ENext : E" << ENext->EID << endl;	  
	  DEBUG << "v0    : V" << v0->VID << endl;
	  DEBUG << "v1    : V" << v1->VID << endl;
	  DEBUG << "v2    : V" << v2->VID << endl;
	  while(!Iter.EndOfList())
	  {
	      DEBUG << " -> E" << Iter.Current()->object->EID << flush;      
	      Iter.GoNext();
	  }
	  DEBUG << endl;
       #endif */
       
    
       if ( Geom::Turnxy( v0, v1, v2 ) == TURN_LEFT && OkTriangle( v0, v1, v2, InflRegnAux ) )
       {
           
   	   // create triangle of vertices v0, v1, v2, i.e., whose
           // boundary edges are ECurr, ENext, and a new edge v2-v0
	   
	   PTEdge NewEdge = new TEdge( v2, v0 );
	   check( (NewEdge == NULL), "TDestroyDelaunay::RIRInitialTrg(), insufficient memory" );
	   
	   NewTriangle = new TTriangle( ECurr, ENext, NewEdge );	   
	   	   
	   // since NewEdge = (v2, v0), NewTriangle is on the left of NewEdge
	   
	   NewEdge->ET[0] = NewTriangle;
	   NewEdge->ET[1] = NULL;

	   // put new edge in queue for RIRDelOptTrg()
	   
	   SwapEdgeQueue.AddHead( NewEdge );
	   NewEdge->Mark( SWAP_EDGE_QUEUE );
	   	   
	   // adjacency ET relations of ECurr and ENext from the "interior"
           // side (the side of the not-yet-triangulated part of the region of
           // influence) must point to NewTriangle. Recall that 
           // ECurr = (v0, v1) and ENext = (v1, v2)

	   if ( ECurr->EV[0] == v0 && ECurr->EV[1] == v1 ) 
              // ...NewTriangle to the left of ECurr
	      ECurr->ET[0] = NewTriangle;
	   else 
	   if ( ECurr->EV[0] == v1 && ECurr->EV[1] == v0 )
              // ...NewTriangle to the right of ECurr
	      ECurr->ET[1] = NewTriangle;
	   else
	      error( "TDestroyDelaunay::RIRInitialTrg(), <> inconsistency detected" );
	      
	   // same thing for ENext
	   
	   if ( ENext->EV[0] == v1 && ENext->EV[1] == v2 )
              // ...NewTriangle to the left ENext
	      ENext->ET[0] = NewTriangle;
	   else 
	   if ( ENext->EV[0] == v2 && ENext->EV[1] == v1 )
              // ...NewTriangle to the right of ENext
	      ENext->ET[1] = NewTriangle;
	   else
	      error( "TDestroyDelaunay::RIRInitialTrg(), <> inconsistency detected" );
	      
           AddTriangle( NewTriangle );
  
           #ifdef DEBUG
	      DEBUG << endl << "creato T" << NewTriangle->TID << "( E"
	            << ECurr->EID << ", E" << ENext->EID << ", E" << NewEdge->EID << " )" << endl;
	      DEBUG << "E" << ECurr->EID << " : ";
	       if ( ECurr->ET[0] == NULL ) DEBUG << "NULL, "; else DEBUG << "T" << ECurr->ET[0]->TID << " ";
	       if ( ECurr->ET[1] == NULL ) DEBUG << "NULL "; else DEBUG << "T" << ECurr->ET[1]->TID;
	       DEBUG << endl;
	      DEBUG << "E" << ENext->EID << " : ";
	       if ( ENext->ET[0] == NULL ) DEBUG << "NULL, "; else DEBUG << "T" << ENext->ET[0]->TID << " ";
	       if ( ENext->ET[1] == NULL ) DEBUG << "NULL "; else DEBUG << "T" << ENext->ET[1]->TID;
	       DEBUG << endl;
	      DEBUG << "E" << NewEdge->EID << " : ";
	       if ( NewEdge->ET[0] == NULL ) DEBUG << "NULL, "; else DEBUG << "T" << NewEdge->ET[0]->TID << " ";
	       if ( NewEdge->ET[1] == NULL ) DEBUG << "NULL "; else DEBUG << "T" << NewEdge->ET[1]->TID;
	       DEBUG << endl;
	  #endif // DEBUG    
	      
	   
	   // remove ECurr, ENext from list InflRegnAux (and remove their
	   // mark INFL_BORDER_AUX), and replace them with NewEdge (and
           // mark it as BORDER_AUX). Recall that
	   // IAuxIter points to the node containing ENext. 
           // In such operations, the list is managed as a circular list,
	   // so we must to pay attention to the case in which ECurr
	   // is the last element of the lista and ENext is the first one.
	   
	   ECurr->UnMark( INFL_BORDER_AUX );
	   ENext->UnMark( INFL_BORDER_AUX );
	   
	   if ( IAuxIter.StartOfList() ) // IAuxIter points to first of list
	   {
	      // ECurr is the last element, ENext is the first one
	      
	      InflRegnAux.RemoveHead(); // remove ENext
	      InflRegnAux.RemoveLast(); // remove ECurr
	      
	      NewEdge->Mark( INFL_BORDER_AUX ); // put NewEdge
	      InflRegnAux.AddHead( NewEdge );
	      
	      IAuxIter.Restart();
	         
	   }
	   else // ! IAuxIter.BeginOfList()
	   {
	      // remove ECurr
	      InflRegnAux.RemoveBefore( IAuxIter.Current() ); 

              // attach NewEdge after ENext, and mark it
	      NewEdge->Mark( INFL_BORDER_AUX );
	      InflRegnAux.AddAfter( IAuxIter.Current(), NewEdge );
	      
	      // go to new edge (NewEdge)	      
	      IAuxIter.GoNext();
	      
	      // remove ENext (now the one before NewCurrent())
	      InflRegnAux.RemoveBefore( IAuxIter.Current() );	      
	   }
	   	   
       }
       
       ECurr = IAuxIter.Current()->object;
	          
    } // end ...while( InflRegnAux.Lenght() > 3 )
    

    
    //
    // Now there are three edges left in list InflRegnAux, create the
    // last triangle with such edges
    //
    
    #ifdef ROBUST
       check( (InflRegnAux.Lenght() != 3), "TDestroyDelaunay::RIRInitialTrg(), <> inconsistency detected" );
    #endif
    
    PTEdge E[3];
    
    for( e=0; e<3; e++ )
    {
        E[e] = InflRegnAux.RemoveHead();
	E[e]->UnMark( INFL_BORDER_AUX );
    }
    
    NewTriangle = new TTriangle( E[0], E[1], E[2] );
    check( (NewTriangle == NULL), "TDestroyDelaunay::RIRInitialTrg(), insufficient memory" );
    
    
    // adjust relations ET of the three edges

    PTVertex v[3];
        
    NewTriangle->GetTV( v[0], v[1], v[2] );
    
    
    // from GetTV, V[i] and V[i+1 % 3] are the vertices of E[i]
    // sorted in counterclockwise order
    
    for( e=0; e<3; e++ )
    {
       if ( E[e]->EV[0] == v[e] && E[e]->EV[1] == v[(e+1)%3] )
          E[e]->ET[0] = NewTriangle; // to the left of E[e]
       else
       if ( E[e]->EV[1] == v[e] && E[e]->EV[0] == v[(e+1)%3] )
          E[e]->ET[1] = NewTriangle; // to the right of E[e]
       else
          error( "TDestroyDelaunay::RIRInitialTrg(), <> inconsistency detected" );
    }
    
    // then adjust the last triangle
            
    AddTriangle( NewTriangle );
    
    #ifdef DEBUG
	DEBUG << endl << "created T" << NewTriangle->TID << "( E"
	            << E[0]->EID << ", E" << E[1]->EID << ", E" << E[2]->EID << " )" << endl;
	 DEBUG << "E" << E[0]->EID << " : ";
	       if ( E[0]->ET[0] == NULL ) DEBUG << "NULL, "; else DEBUG << "T" << E[0]->ET[0]->TID << " ";
	       if ( E[0]->ET[1] == NULL ) DEBUG << "NULL "; else DEBUG << "T" << E[0]->ET[1]->TID;
	       DEBUG << endl;
	 DEBUG << "E" << E[1]->EID << " : ";
	       if ( E[1]->ET[0] == NULL ) DEBUG << "NULL, "; else DEBUG << "T" << E[1]->ET[0]->TID << " ";
	       if ( E[1]->ET[1] == NULL ) DEBUG << "NULL "; else DEBUG << "T" << E[1]->ET[1]->TID;
	       DEBUG << endl;
	 DEBUG << "E" << E[2]->EID << " : ";
	       if ( E[2]->ET[0] == NULL ) DEBUG << "NULL, "; else DEBUG << "T" << E[2]->ET[0]->TID << " ";
	       if ( E[2]->ET[1] == NULL ) DEBUG << "NULL "; else DEBUG << "T" << E[2]->ET[1]->TID;
	       DEBUG << endl;
    #endif // DEBUG    

}



// -----------------------------------------------------------------------------
//  
//   boolean TDestroyDelaunay::OkTriangle( PTVertex v0, v1, v2, 
//                                       TDoubleList<PTEdge> &InflRegnAux )
//
//   Boolean function, return true iff it is possible to create
//   a triangle of vertices v0, v1, v2, w.r.t. the region whose boundary
//   is in list InflRegnAux (the same as the InflRegnAux list
//   in the calling function RIRInitialTrg()).
//
//   Check that each vertex, different from v0, v1 and v2, which is
//   an endpoint of an edge in InflRegnAux, does not fall inside
//   the potential triangle v0, v1, v2.
//
//   This check makes the initial retriangulation of the influence region 
//   a quadratic time complexity in the dimension of the boundary
//   of the region. But influence regions are typically small...
//

int TDestroyDelaunay::OkTriangle( PTVertex v0, PTVertex v1, PTVertex v2,
                                    TDoubleList<PTEdge>& InflRegnAux )
{
    #ifdef DEBUG
       DEBUG << "TDestroyDelaunay::OkTriangle( V"
             << v0->VID << ", V" << v1->VID << ", V" << v2->VID
             << " )" << endl;
       
    #endif       

    #ifdef ROBUST
       check ( (InflRegnAux.Lenght() <= 3), 
           "TDestroyDelaunay::OkTriangle(), <0> inconsistency detected" );
    #endif
    
    PTEdge EPrev = InflRegnAux.GetLast(), ECurr = NULL;
    
    TDoubleListIterator<PTEdge> Iter( &InflRegnAux );
           
    Iter.Restart();   
    while ( !Iter.EndOfList() )
    {
  
       ECurr = Iter.Current()->object;
       
       // v is the vertex of ECurr which is not common to EPrev
       // (the "second" vertex of ECurr)
       
       PTVertex v = ( ECurr->EV[0] != EPrev->EV[0] && ECurr->EV[0] != EPrev->EV[1] ?
                ECurr->EV[0] :
	        ECurr->EV[1] );
		
       #ifdef ROBUST
          // check if the other veretex of ECurr is really common to EPrev
          PTVertex vx = ( v != ECurr->EV[0] ? ECurr->EV[0] : ECurr->EV[1] );
	  check( (vx != EPrev->EV[0] && vx != EPrev->EV[1]), 
	     "TDestroyDelaunay::OkTriangle(), <1> inconsistency detected" );
       #endif

       if ( v!=v0 && v!=v1 && v!=v2 )
          if ( Geom::InTrianglexy( v0, v1, v2, v ) )
	     return( FALSE );
	     
       // else go to next vertex
       
       EPrev = ECurr;
       Iter.GoNext();
         
    }
    
    return( TRUE );
}

// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::RIRDelOptTrg()
//
//   Another auxiliary function used by RetriangulateInfluenceRegion().
//   Its name stands for: (R)etriangulate (I)nfluence (R)egion
//   (Del)aunay (Opt)imize (Tr)iangulation.
//
//   Optimize the initial triangulation built by RIRInitialTrg() inside
//   the influence region, according to the delaunay criterion.
//   This is done through a sequence of edge swaps.
//
//   INPUT: list SwapEdgeQueue, initially containing the set of edges
//   created in RIRInitialTrg().
//
//   OUTPUT: a consistent status of TDestroyDelaunay.
//   FirstTriangle is modified by calls to AddTriangle in such a way that
//   it always points to the last created triangle.
//

void TDestroyDelaunay::RIRDelOptTrg()
{
    int e;

    #ifdef DEBUG
       DEBUG << "TDestroyDelaunay::RIRDelOptTrg()" << endl;
    #endif
    
    PTEdge ECurr, NewEdge, E[4];
    
    while( !SwapEdgeQueue.IsEmpty() )
    {
	ECurr = SwapEdgeQueue.RemoveHead();
	ECurr->UnMark( SWAP_EDGE_QUEUE );
	
	if ( EdgeToSwap( ECurr ) )
	{
	    // swap edge
	    
	    NewEdge = SwapEdge( ECurr );
	    
	    // enqueue the boundary edges of the quadrilateral
	    // formed by the two triangles adjacent to NewEdge
	       
	    GetQuadBorder( NewEdge, E[0], E[1], E[2], E[3] );
	    
	    for( e=0; e<4; e++ )
	    {
		if ( !E[e]->Marked( SWAP_EDGE_QUEUE ) && !E[e]->Marked( INFL_BORDER ) )
		{
		    E[e]->Mark( SWAP_EDGE_QUEUE );
		    SwapEdgeQueue.AddTail( E[e] );
		}
	    }
	    
	} // end ...if( EdgeToSwap )
	
    } // end ...while( !SwapEdgeQueue.IsEmpty() )
    
}


// -----------------------------------------------------------------------------
//  
//   boolean TDestroyDelaunay::EdgeToSwap( PTEdge E )
//
//   Return TRUE iff the given edge does not satisfy the Delaunay
//   criterion, i.e.:
//
//    - quadrilateral formed by the union of the two triangles
//      adjacent to edge E is strictly convex, and
//
//    - the circum-circle of one of the two triangles contains the
//      vertex opposite to E of the other triangle
//

boolean TDestroyDelaunay::EdgeToSwap( PTEdge E )
{
    
    int i;

    // the two vertices of the quadrilateral which are endpoints of E
    
    PTVertex ve0, ve1;
    
    ve0 = E->EV[0];
    ve1 = E->EV[1];
    
    // the other two vertices of the quadrilateral

    PTTriangle T0 = E->ET[0];
    PTTriangle T1 = E->ET[1];    
   
    PTVertex vo0 = NULL, vo1 = NULL;
    
    PTVertex V[3];
    
    // find vo0
    T0->GetTV( V[0], V[1], V[2] );
    for( i=0; i<3; i++ )
       if ( V[i] != ve0 && V[i] != ve1 )
          vo0 = V[i];
	  
    // find vo1
    T1->GetTV( V[0], V[1], V[2] );
    for( i=0; i<3; i++ )
       if ( V[i] != ve0 && V[i] != ve1 )
          vo1 = V[i];
    
    #ifdef ROBUST
       check( (vo0 == NULL || vo1 == NULL), 
         "TDestroyDelaunay::EdgeToSwap(), inconsistency detected" );
    #endif
    
    
    // check turns
    
    int Turn0 = Geom::Turnxy( vo0, ve0, vo1 );
    int Turn1 = Geom::Turnxy( vo0, ve1, vo1 );
    
    
    // cases in which quadrilateral is not strictly convex
    
    if ( Turn0 == ALIGNED || Turn1 == ALIGNED || Turn0 == Turn1 )
       return( FALSE );
   
   
    // otherwise, quadrilateral is strictly convex
    
    return( T0->InCircle( vo1 ) || T1->InCircle( vo0 ) );
         
}


// -----------------------------------------------------------------------------
//  
//   PTEdge TDestroyDelaunay::SwapEdge( PTEdge &E )
//
//   Swap edge E with the opposite edge of the quadrilateral formed by the
//   two triangles adjacent to E.
//   
//   Note (1): This function is called only during the optimization of the
//   influence region, for internal edges of the region. 
//   Since, before that, VE adjacencies of all vertices of the region
//   boundary have been set to point edges of the same boundary, there is
//   no danger that one of the two endpoints of E may point to E through
//   its VE relation.
//
//   Note (2): Pointer E is passed by reference in order to set it to NULL
//   after deleting the edge it points. This makes the code more robust.
//
//   The pointer to the edge created to replace E is returned by the
//   function.
//
 
PTEdge TDestroyDelaunay::SwapEdge( PTEdge &E )
{

    #ifdef DEBUG
       DEBUG << "TDestroyDelaunay::SwapEdge(E" << E->EID << ")" << endl;
    #endif
    
    PTTriangle OldT0 = E->ET[0];
    PTTriangle OldT1 = E->ET[1];
    
    #ifdef ROBUST
       check( (OldT0 == NULL || OldT1 == NULL), "TDestroyDelaunay::SwapEdge() <1> inconsistency detected" );
    #endif

    PTEdge E0a, E0b, E1a, E1b;
    GetQuadBorder( E, E0a, E0b, E1a, E1b );    
    
    // vo0 and vo1 are the vertices of OldT0 and OldT1 (respectively),
    // not endpoints of E
    // (vo0 is the common vertex of E0a and E0b, and vo1 of E1a and E1b)
        
    PTVertex vo0 = NULL, vo1 = NULL;
    
    vo0 = ( (E0a->EV[0] == E0b->EV[0] || E0a->EV[0] == E0b->EV[1]) ?
            E0a->EV[0] :
	    E0a->EV[1] );
	    
    vo1 = ( (E1a->EV[0] == E1b->EV[0] || E1a->EV[0] == E1b->EV[1]) ?
            E1a->EV[0] :
	    E1a->EV[1] );
	    
    #ifdef ROBUST
       check( ( vo0 != E0b->EV[0] && vo0 != E0b->EV[1] ),
            "TDestroyDelaunay::SwapEdge(), <2> inconsistency detected");
       check( ( vo1 != E1b->EV[0] && vo1 != E1b->EV[1] ),
            "TDestroyDelaunay::SwapEdge(), <3> inconsistency detected");
       check( (vo0 == NULL || vo1 == NULL),
            "TDestroyDelaunay::SwapEdge(), <4> inconsistency detected");
    #endif // ROBUST
               
    PTEdge NewE = new TEdge( vo0, vo1 ); // directed from vo0 to vo1
    check( (NewE == NULL), "TDestroyDelaunay::SwapEdge(), insufficient memory" );
    
    PTTriangle NewT0 = new TTriangle( NewE, E1b, E0a );
    PTTriangle NewT1 = new TTriangle( NewE, E0b, E1a );
    check( (NewT0 == NULL || NewT1 == NULL), "TDestroyDelaunay::SwapEdge(), insufficient memory" );
    
    // adjust ET relations of NewE, E0a, E0b, E1a, E1b
    
    NewE->ET[0] = NewT0;
    NewE->ET[1] = NewT1;
    
    if ( E0a->ET[0] == OldT0 ) E0a->ET[0] = NewT0;
    else 
    if ( E0a->ET[1] == OldT0 ) E0a->ET[1] = NewT0;
    else
    error( "TDestroyDelaunay::SwapEdge(), <5> inconsistency detected");
    
    if ( E0b->ET[0] == OldT0 ) E0b->ET[0] = NewT1;
    else 
    if ( E0b->ET[1] == OldT0 ) E0b->ET[1] = NewT1;
    else
    error( "TDestroyDelaunay::SwapEdge(), <6> inconsistency detected");
    
    if ( E1a->ET[0] == OldT1 ) E1a->ET[0] = NewT1;
    else 
    if ( E1a->ET[1] == OldT1 ) E1a->ET[1] = NewT1;
    else
    error( "TDestroyDelaunay::SwapEdge(), <7> inconsistency detected");
    
    if ( E1b->ET[0] == OldT1 ) E1b->ET[0] = NewT0;
    else 
    if ( E1b->ET[1] == OldT1 ) E1b->ET[1] = NewT0;
    else
    error( "TDestroyDelaunay::SwapEdge(), <8> inconsistency detected");
    
    DetachEdge( E );
    
    delete(E); E = NULL;
    
    DetachTriangle( OldT0 );
    delete( OldT0 );
    DetachTriangle( OldT1 );
    delete( OldT1 );
    
    AddTriangle( NewT0 );
    AddTriangle( NewT1 );
    
    #ifdef DEBUG 
       DEBUG << "T" << NewT0->TID 
                    << ": E" << NewT0->TE[0]->EID 
		    << ", E" << NewT0->TE[1]->EID 
		    << ", E" << NewT0->TE[2]->EID << endl;
       DEBUG << "T" << NewT1->TID 
                    << ": E" << NewT1->TE[0]->EID 
		    << ", E" << NewT1->TE[1]->EID 
		    << ", E" << NewT1->TE[2]->EID << endl;
       DEBUG << "E" << NewE->EID 
             << ": T" << NewE->ET[0]->TID 
	     << ", T" << NewE->ET[1]->TID 
	     << endl;
    #endif
    
    return( NewE );

}


// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::GetQuadBorder()
//
//   Given an edge E, return the four boundary edges E0a, E1a, E0b, E1b 
//   of the quadrilateral formed by the two triangles adjacent to E to
//   the right and to the left.
//   If E is on the convex hull, two of such pointers are set to NULL,
//   but this function is always called on internal edges...
//   The four edges form the boundary of the quadrilateral in
//   counterclockwise order, starting from the "second" endpoint of E
//   (E->EV[1]).
//

void TDestroyDelaunay::GetQuadBorder( PTEdge E, PTEdge &E0a, PTEdge &E0b,
				    PTEdge &E1a, PTEdge &E1b )
{
   
   int e;

   E0a = E0b = E1a = E1b = NULL;
   
   #ifdef ROBUST
     check( (E==NULL), "TDestroyDelaunay::GetQuadBorder(), <1> inconsistency detected" );
   #endif
   
   PTTriangle T0 = E->ET[0];
   PTTriangle T1 = E->ET[1];
      
   if ( T0 != NULL )
   {
        for( e=0; e<3; e++ ) if ( T0->TE[e] == E ) break;

        #ifdef ROBUST
          check( (e>=3), "TDestroyDelaunay::GetQuadBorder(), <2> inconsistency detected" );
        #endif
    
        E0a = T0->TE[(e+1)%3];
        E0b = T0->TE[(e+2)%3];
    }
    
    if ( T1 != NULL )
    {
    
        for( e=0; e<3; e++ ) if ( T1->TE[e] == E ) break;

        #ifdef ROBUST
           check( (e>=3), "TDestroyDelaunay::SwapEdge(), <2> inconsistency detected" );
        #endif
    
        E1a = T1->TE[(e+1)%3];
        E1b = T1->TE[(e+2)%3];
    }
 
    
}
   
     

// -----------------------------------------------------------------------------
//  
//   void TDestroyDelaunay::DeleteInfluenceRegion()
//
//   Redefined to obtain that, at the end of each update step, before
//   deleting the content of list InflRegnBorder, we recompute the
//   degree of vertices belonging to such sequence of edges.
//   Such vertices are taken away from tree ElimVtxTree, and inserted 
//   again after recomputing their degree.
//

void TDestroyDelaunay::DeleteInfluenceRegion()
{
   #ifdef DEBUG
      DEBUG << "TDestroyDelaunay::DeleteInfluenceRegion()" << endl;
   #endif

   //
   // Call the version of DeleteInfluenceRegion() provided in the base class
   // to delete the region of influence
   //

   // it was at the end: moved here in such a way that ReCheckVertex(V)
   // can see the correct number of nIncConstr that has been decremented
   // by DetachEdge().
   TDelaunayBase::DeleteInfluenceRegion(); 

   TDoubleListIterator<PTEdge> IRegnIter( &InflRegnBorder );
   IRegnIter.Restart();

   PTEdge   ECurr;
   PTEdge   EPrev = InflRegnBorder.GetLast();
   
   //
   // Recall that InflRegnBorder, at this step, always represent a
   // closed sequence of edges (a new edge has been created in
   // TDestroyDelaunay::RetriangulateInfluenceRegion(), if it was open)
   //

   PTVertex VCurr = NULL;

   while( !IRegnIter.EndOfList() ) 
   {   
      
      ECurr = IRegnIter.Current()->object;
      
      //
      // VCurr is the common vertex of EPrev and ECurr, i.e., the second
      // vertex, in counterclockwise order w.r.t. VertexTORemove, of EPrev
      //
      
      if ( EPrev->EV[0] == ECurr->EV[0] || EPrev->EV[0] == ECurr->EV[1] )
        VCurr = EPrev->EV[0];
      else
        VCurr = EPrev->EV[1];
	
      #ifdef ROBUST
        check( (ECurr->EV[0] != VCurr && ECurr->EV[1] != VCurr),
	       "TDestroyDelaunay::DeleteInfluenceRegion() <0> inconsistency detected" );	
      #endif
		
      //
      // remove vertex VCurr from the tree, call ReCheckVertex, and put
      // it in the tree again;
      // VCurr may not be in the tree, since until now it was not removable,
      // but changing its degree it may have become removable
      //
       
      if ( ElimVtxTree.IsIn( VCurr ) )
      {
         PTVertex VN = ElimVtxTree.Remove( VCurr );
            
         #ifdef ROBUST
            check( (VN != VCurr), "TDestroyDelaunay::DeleteInfluenceRegion(), <1> inconsistency detected" );
         #endif
      }

      if ( ReCheckVertex( VCurr ) )
      {
         ElimVtxTree.Insert( VCurr );
      }
      
      //
      // the next one
      //
      
      EPrev = ECurr;
      
      IRegnIter.GoNext();
      
   } // end ...while( !IRegnIter.EndOfList() )
   
}



// -----------------------------------------------------------------------------
//  
//   boolean TDestroyDelaunay::ReCheckVertex( PTVertex V )
//
//   Check if a vertex is removable after the retriangulation of the 
//   influence region on the boundary of which the vertex is.
//

boolean TDestroyDelaunay::ReCheckVertex( PTVertex V )
{
   return( IsVtxElim(V) && OkDegree( CalcDegree(V) ) );
}

