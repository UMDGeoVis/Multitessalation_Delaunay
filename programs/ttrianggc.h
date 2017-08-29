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

// --------------------------------------------------------------------------
//
//  file   : ttrianggc.h
//  author : Christian Melchiorre
//
//  Definition of class GC, which provides a number of static functions 
//  to allocate / deallocate entities defined in ttriang.h
//  (TTriangle/TEdge...). For such entities, several consecutive new/delete
//  operations are performed during the process of computing a
//  Delaunay triangulation.
//
//  Only the version of New/Delete for TTriangle/TEdge has been implemented.
//  For other entities (TPoint/TVertex), Garbage Collection is not needed,
//  because no consecutive new/delete operations are performed.
//
//  The GC management is a bit different from that performed for template
//  container classes (TList/TDoubleList/TBTree). This is due to the fact
//  that entities involved here have not a field of type pointer to 
//  entity of the same type, which can be used as next fiels in the GCList.
//  Therefore, we use an array of pointers to the involved entities.
//


#ifndef _GC_H
#define _GC_H

#include "defs.h"
#include "error.h"
#include "ttriang.h"


#ifdef _GC_ON


//
// Costants for the maximum dimension of GCxxxLists. It must be set
// depending on the application issues.
//

#define GC_TRG_CAPACITY 255
#define GC_EDG_CAPACITY 255


class GC
{
   private:

      // Array of pointers to the Garbage Collected structures
      static PTTriangle TrgArray[GC_TRG_CAPACITY];            
      
      // Next free entry in the array
      static int NextFreeTrg;                                 	     
	     
      // Array of pointers to the Garbage Collected structures
      static PTEdge EdgArray[GC_EDG_CAPACITY];
      
      // Next free entry in the array
      static int NextFreeEdg;
            
   public:
    
      static PTTriangle NewTriangle( PTEdge, PTEdge, PTEdge );
      static void DeleteTriangle( PTTriangle );
      
      static PTEdge NewEdge( PTVertex, PTVertex );
      static void DeleteEdge( PTEdge );

};


#endif // _GC_ON

#endif // _GC_H

