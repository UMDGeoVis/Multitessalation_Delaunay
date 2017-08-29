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
//   file   : refrnddel.cpp
//   author : Christian Melchiorre
//
//   Implementation of class TRefRndDelaunay, sub-class of TRefineDelaunay,
//   for building a Delaunay triangulation through refinement, with a
//   random choice of the point to be inserted.
//   The only difference With respect to the base class, is that points
//   not belonging to the convex hull are sorted based on a random key.
//   Thus, method InitialTriangulation() ius redefined: it performs
//   such random sorting after the operation that it performed in the
//   base class.
//

#include <iostream>

#ifdef CC_SILICON /* PAOLA */
#include <rand48.h>
#else
#include <stdlib.h>
#endif /* CC_SILICON */
#include <time.h>

#include "defs.h"
#include "error.h"

#include "refdel.h"
#include "refrnddel.h"

// aggiunta del tra
#ifdef CC_VISUAL5
#define srand48 srand
#define lrand48() ((long)rand())
#endif
// fine aggiunta

   
// -----------------------------------------------------------------------------------------
//
//  void TRefRndDelaunay::InitialTriangulation()
//
//  Call the version of InitialTriangulation() profided by the base class,
//  and, after that, re-sort remaining points (after making the initial
//  triangulation of the convex hull) according to a random key,
//  by using function TRefRndDelaunay::RandomSortPoints().
//

void TRefRndDelaunay::InitialTriangulation()
{
    
    #ifdef DEBUG
       DEBUG << "Chiamata TRefRndDelaunay::InitialTriangulation()" << endl;
    #endif // DEBUG
           
    //
    // The following call to TRefineDelaunay::InitialTriangulation() 
    // finds the convex hull of the input points and builds a triangulation
    // of the convex hull. At the end, the pointers to the vertices of
    // the convex hull are found in  positions  [0...nChPts-1] of array
    // Points, while other points are in positions [nChPts...nPts-1].
    //
    
    TRefineDelaunay::InitialTriangulation();
    
    //
    // RandomSortPoints() sorts the points contained in positions
    // [nChPts...nPts-1] according to a random key
    //
    
    RandomSortPoints();
    
}
    

// -----------------------------------------------------------------------------------------
//
//  void TRefRndDelaunay::RandomSortPoints()
//
//  Sort, according to a random key, the points of array Points which do
//  not belong to the convex hull, i.e., those points pointed by 
//  Points[nChPts..nPts-1] after calling CalcConvexHull.
//
//  The random sorting is obtained by simply making a number 
//  n = nPts - nChPts of random swaps between elements of the relevant
//  part of the array.
//  The alternative is a call to qsort with a random key.
//  This method has a lower time complexity, O( n ) w.r.t. O( n log n )
//  in the average case of qsort.
//  Moreover, it is not necessary to store the random key.
//  Results are comparable.
//

void TRefRndDelaunay::RandomSortPoints()
{

   #ifdef DEBUG
     DEBUG << "Chiamata TRefRndDelaunay::RandomSortPoints()" << endl;
   #endif
   
   int n = nPts - nChPts;
   int j;
   
   PTPoint pTmp;
   
   srand48( (int)time(NULL) );
   
   for( int i=nChPts; i<nPts; i++ )
   {

      j = nChPts + (int)( lrand48() % (long)n );
      // j is another random number in [nChPts...nPts-1]

      #ifdef ROBUST  
         check( (j<nChPts || j>=nPts),
            "TRefRndDelaunay::RandomSortPoints(), index out of bounds" );
      #endif  

      //
      // swap Points[i] and Points[j]
      //
      
      pTmp = Points[i];
      Points[i] = Points[j];
      Points[j] = pTmp;
   }
   
}

