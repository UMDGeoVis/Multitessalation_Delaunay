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

// --------------------------------------------------------------------------------
//
//  file   : basedel.cpp
//  author : Christian Melchiorre
//
//  Implementation of base class TDelaunayBase.
//

#include <iostream>

#include "defs.h"
#include "error.h"
#include "utils.h" // PAOLA for CheckIsolatedPoint
#include "tdoublelist.h"

#include "ttriang.h"
#include "ttrianggc.h"
#include "basedel.h"

// --------------------------------------------------------------------------------
//  
//  Constructor of TDelaunayBase.
//

TDelaunayBase::TDelaunayBase() : TTriangulation(), InflRegnBorder() 
{
   #ifdef DEBUG
      DEBUG << "TDelaunayBase Constructor" << endl;
   #endif
 
   FirstTrgToDel = NULL; 
}
       
// --------------------------------------------------------------------------------
//  
//  void TDelaunayBase::DeleteInfluenceRegion()
//
//  Delete triangles and edges internal to the influence region.
//

void TDelaunayBase::DeleteInfluenceRegion()
{

     #ifdef DEBUG
      DEBUG << "\nTDelaunayBase::DeleteInfluenceRegion( "
            << "FirstTrgToDel = ";
      if ( FirstTrgToDel == NULL ) DEBUG << "NULL"; else 
      DEBUG << "T" << FirstTrgToDel->TID;
      DEBUG << " )" << endl;
     #endif // DEBUG

     // ...if there are no triangles to be deleted
     if ( FirstTrgToDel == NULL ) return; 
     
     TDoubleList<PTTriangle> TrgToDel;
     PTTriangle CurrTrg, NextTrg;
     PTEdge NextEdg;
          
     TrgToDel.AddHead( FirstTrgToDel );
     
     while ( ! TrgToDel.IsEmpty() )
     {
 	    
         CurrTrg = TrgToDel.RemoveHead();
	 
	 #ifdef ROBUST
  	    check( (CurrTrg == NULL), "TDelaunayBase::DeleteInfluenceRegion(), <1> this should not happen");
	 #endif
	 	  
	 for( int e=0; e<3; e++ )
	 {	    
	     NextEdg = CurrTrg->TE[e];
	     
	     //
	     // edge may have been deleted in previous iterations of the
	     // while loop, must check
	     //
	     
	     if ( NextEdg != NULL ) 
	     {	     
	     
		if ( NextEdg->Marked(TO_DELETE) )
		{
                    NextTrg = ( NextEdg->ET[0] != CurrTrg ? NextEdg->ET[0] : NextEdg->ET[1] ); 
		    
                    DetachEdge( NextEdg );
	            
		    #ifdef _GC_ON
		       GC::DeleteEdge( NextEdg );
		    #else
   		       delete( NextEdg );
		    #endif // _GC_ON
		    
		    NextEdg = NULL;
		}
		else
		    NextTrg = NULL;
		
		
	        if ( NextTrg != NULL )
		{
		   if ( NextTrg->Marked( TO_DELETE ) )
		   {
		       NextTrg->UnMark( TO_DELETE );
		       TrgToDel.AddTail( NextTrg );
		   }
		}
						     
	     } // end ...if( NextEdg != NULL )...
	     
	     
	 } // end...for
	 
	 //
	 // detach CurrTrg from triangulation. The destructor of TTriangle
	 // set to NULL the pointers to CurrTrg from edges adjacent to it.
	 //
	 	 
	 DetachTriangle( CurrTrg );
	 
	 #ifdef _GC_ON
	    GC::DeleteTriangle( CurrTrg );
	 #else
   	    delete( CurrTrg );
	 #endif // _GC_ON
	 
	 CurrTrg = NULL;

	    
     } // end ...while( !TrgToDel.IsEmpty() )
     
      #ifdef DEBUG
         DEBUG << "Exit DeleteInfluenceRegion" << endl;
      #endif 
}
