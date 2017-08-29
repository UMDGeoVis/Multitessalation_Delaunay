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
//   author : Christian Melchiorre
//
//  Implementation of class, which provides a number of static functions 
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


#include "defs.h"
#include "error.h"
#include "ttriang.h"
#include "ttrianggc.h"

#ifdef _GC_ON

//
// Init (static) dimension of GCList
//

int GC::NextFreeTrg = 0;
int GC::NextFreeEdg = 0;

PTTriangle GC::TrgArray[];
PTEdge GC::EdgArray[];


// --------------------------------------------------------------------------
//
//  new TTriangle
//

PTTriangle GC::NewTriangle( PTEdge E0, PTEdge E1, PTEdge E2 )
{
    PTTriangle NewTriangle = NULL;
    
    if ( GC::NextFreeTrg > 0 )
    {
       //
       // there is at least one element in GCList TrgArray[]
       //
       
       NextFreeTrg--;
       
       NewTriangle = TrgArray[ NextFreeTrg ];
       
       #ifdef ROBUST
          check( (NewTriangle == NULL), "GC::NewTriangle(), <1> this should not happen" );
       #endif
	  
       TrgArray[ NextFreeTrg ] = NULL;
       
       NewTriangle->TE[0] = E0;
       NewTriangle->TE[1] = E1;
       NewTriangle->TE[2] = E2;
         
       NewTriangle->TID = TTriangle::NextTID++;
    
       NewTriangle->CalcCircle();
   
       NewTriangle->Error = 0.0;
       
       NewTriangle->MarkReset();
       
       check( (!NewTriangle->PointList.IsEmpty()),
           "GC::NewTriangle(), found a Garbage Collected TTriangle with non void PointList");
	   
    }
    else // NextFreeTrg == 0
    {
       NewTriangle = new TTriangle( E0, E1, E2 );
    }
    
    check( (NewTriangle == NULL), "GC::NewTriangle(), insufficient memory" );
    
    return( NewTriangle );

}



// --------------------------------------------------------------------------
//
//  delete TTriangle
//

void GC::DeleteTriangle( PTTriangle T )
{

   if ( NextFreeTrg >= GC_TRG_CAPACITY )
   {
       delete(T);
   }
   else
   {
       TrgArray[ NextFreeTrg ] = T;
       
       NextFreeTrg++;
   }
}   
  
  
// --------------------------------------------------------------------------
//
//  new TEdge
//
 
PTEdge GC::NewEdge( PTVertex V0, PTVertex V1 )
{
    PTEdge NewEdge = NULL;
    
    if ( GC::NextFreeEdg > 0 )
    {
       //
       // there is at least one element in GCList EdgArray[]
       //
       
       NextFreeEdg--;
       
       NewEdge = EdgArray[ NextFreeEdg ];
       
       #ifdef ROBUST
          check( (NewEdge == NULL), "GC::NewEdge(), <1> this should not happen" );
       #endif
	  
       EdgArray[ NextFreeEdg ] = NULL;

       NewEdge->EID = TEdge::NextEID++;

       NewEdge->EV[0] = V0;
       NewEdge->EV[1] = V1;

       NewEdge->ET[0] = NULL;
       NewEdge->ET[1] = NULL;

       NewEdge->Error = 0.0;          
       
       NewEdge->MarkReset();
       
       check( (!NewEdge->PointList.IsEmpty()),
           "GC::NewTriangle(), found a Garbage Collected TTriangle with non void PointList");
	   
    }
    else // NextFreeEdg == 0
    {
       NewEdge = new TEdge( V0, V1 );
    }
    
    check( (NewEdge == NULL), "GC::NewEdge(), insufficient memory" );
    
    return( NewEdge );
    
}

 
// --------------------------------------------------------------------------
//
//  delete TEdge
//

void GC::DeleteEdge( PTEdge E )
{

   if ( NextFreeEdg >= GC_EDG_CAPACITY )
   {
       delete(E);
   }
   else
   {
       EdgArray[ NextFreeEdg ] = E;
       
       NextFreeEdg++;
   }
}   
   
     
#endif // _GC_ON
