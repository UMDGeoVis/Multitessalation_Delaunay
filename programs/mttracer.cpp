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

// -----------------------------------------------------------------------------------------
//
//  file   : mttracer.cpp
//   author : Christian Melchiorre
//
//  Implementation of static class MTTracer. It acts as an interface between
//  the Delaunay triangulation programs and the MT History Tracer for
//  building a multiresolution model (MT).
//  It implements the termination criteria for the process of triangulation
//  update, calls the functions of the MT library, and takes track of
//  the error of the current triangulation and of the mumber of updates
//  performed so far.
//


#ifdef MT_TRACER
#include "fieldmt.h"
#endif

#include "tbtree.h"
#include "ttriang.h"
#include "mttracer.h"
#include "decdel.h"    // for TDecimDelaunay::InitialPhase


// for TBTree Errors (used only if norm = NORM_MAX)

int compare( double a , double b )
{
   if ( a<b ) return -1;
   else
   if ( a>b ) return +1;
   else
   return( 0 );
}



// -----------------------------------------------------------------------------------------
//
//  Constructor of class MTTracer.
//

#ifdef MT_TRACER
MTTracer::MTTracer(void) : FieldWithErrorBuildingInterfaceClass()
#else
MTTracer::MTTracer(void)
#endif
{
   TotError = 0.0;
 
   nTrgs = 0;
   nUpd = -1;
   // trick to ensure that the first time we call MTTracer::MeshOK()
   // (after creation of the initial triangulation), nUpd is correctly
   // set to zero before starting the update process.

   // ...by default
   Term = NO_TERM;
    
}



// -----------------------------------------------------------------------------------------
//
//  Functions to interface with MT History Tracer library
//

void MTTracer::StartHistory ( int HistoryType )
{
   #ifdef MT_TRACER

      MT_StartHistory( 3, 2, HistoryType );
      StartTileErrorHistory();

   #endif
   
   Type = HistoryType;
}


void MTTracer::MakeTriangle( PTTriangle T )
{
   double Error = T->GetError();
   
   PTVertex V[3];

   T->GetTV( V[0], V[1], V[2] );

   #ifdef MT_DEBUG
   {
      int i;
      for(i=0; i<3; i++ )
         MT_DEBUG << "V" << V[i]->VID << " : ( " << (*(PTPoint)(V[i])) << " ) " << endl;
      MT_DEBUG << "MTTracer::MakeTriangle( T" << T->TID << ", Error = " << Error << " )" << endl;
   }
   #endif

   #ifdef MT_TRACER
   {
      float coord[3];
      int i;
      for( i=0; i<3; i++ )
      {  coord[0] = V[i]->x;
         coord[1] = V[i]->y;
         coord[2] = V[i]->z;
         V[i]->MTIdx = MT_UseVertex( V[i]->MTIdx, coord );
      }
      T->MTIdx = MT_MakeTile();
      MakeTileError( T->MTIdx, Error );
   }
   #endif
   
   if ( Term == TERM_ERR )
      UpdateError( ADDTRG, Error );

}



void MTTracer::KillTriangle( PTTriangle T )
{
  
   double Error = T->GetError();
   
   #ifdef MT_DEBUG
      PTVertex TV[3];
      T->GetTV( TV[0], TV[1], TV[2] );

      MT_DEBUG << "MTTracer::KillTriangle( T" << T->TID << " (V" << TV[0]->VID << " ,V" << TV[1]->VID
               << " ,V" << TV[2]->VID << "), Error = " << Error << " )" << endl;
   #endif
   
   #ifdef MT_TRACER
      MT_KillTile( T->MTIdx );
   #endif
   
   if ( Term == TERM_ERR )
      UpdateError( DELTRG, Error );
   
}


void MTTracer::MeshOk()
{
   if( Term == TERM_NUPD  /* &&  TDecimDelaunay::InitialPhase == FALSE */ )
      nUpd++;
   
   #ifdef MT_DEBUG
      MT_DEBUG << "MTTracer::MeshOk(";
      if ( Term == TERM_NUPD ) MT_DEBUG << "nUpd = " << nUpd << ", UpdLev = " << UpdLev;
      else 
      if ( Term == TERM_ERR ) MT_DEBUG << "TotError" << TotError << ", ErrLev = " << ErrLev; 
      MT_DEBUG << ")" << endl;
   #endif

   #ifdef MT_TRACER
      MT_EndUpdate();
   #endif
}


void MTTracer::EndHistory()
{
   #ifdef MT_TRACER
      FieldWithError mt;
      FILE * fd;
   #endif

   #ifdef MT_DEBUG
      MT_DEBUG << "MT::EndHistory()" << endl;
   #endif
   
   #ifdef MT_TRACER
      MT_EndHistory(); 
      EndTileErrorHistory();

      mt = new FieldWithErrorClass(3,2);
      MT_SetTarget(mt);
      SetTargetTileErrorTable(mt);

      MT_Convert();
      ConvertTileErrors();

      fd = fopen("output.mtf","w");
      mt->MT_SetDescription(MT_TheDescription());
      mt->SetTileErrorDescription(MT_TheDescription());
      mt->MT_Write(fd);
      fclose(fd);
      fd = fopen("output.err","w");
      mt->WriteTileErrors(fd);
      fclose(fd);
      delete mt;
   #endif
}


// -----------------------------------------------------------------------------------------
//
//  Functions for maintaining the error of the triangulation, and
//  for termination.
//  
//


void MTTracer::SetTerminateCondition( int iTerm, int iUpdLev )
{
   check( (iTerm != TERM_NUPD), "MTTracer::SetTerminateCondition(), invalid arguments" );
   Term = iTerm;   
   UpdLev = iUpdLev;
      
}

void MTTracer::SetTerminateCondition( int iTerm, int iNorm, double iErrLev )
{
   check( (iTerm != TERM_ERR), "MTTracer::SetTerminateCondition(), invalid arguments" );
   Term   = iTerm;
   Norm   = iNorm;
   ErrLev = iErrLev;
}

boolean MTTracer::TerminateCondition()
{
   
   switch( Term )
   {
   
      case NO_TERM: return(FALSE);
      
      case TERM_NUPD:
           {
	      return( nUpd >= UpdLev );
	   }
	   
      case TERM_ERR:
           {
	      switch( Type )
	      {
	          case MT_REFINING: return( TotError <= ErrLev );
		  case MT_COARSENING: return( TotError >= ErrLev );
		  default:
                    error( "MTTracer::TerminateCondition(), <1> inconsistency detected" );
   	      }
	   }
	   
      default: 
         error( "MTTracer::TerminateCondition(), <2> inconsistency detected" );
	 
   }

   return( FALSE ); // ...to avoid warning

}


void MTTracer::UpdateError( int op, double Error )
{

   switch( op )
   {
      case ADDTRG:
           {
	      nTrgs++;
	      
	      switch( Norm )
	      {
	         case NORM_MAX:
		      {
		         Errors.Insert( Error );
			 TotError = Errors.GetMax();
		      }
		      break;
		      
		 case NORM_MED:
		      {
		         if ( nTrgs == 1 )
			   TotError = Error;
			 else
                           TotError = ( ( TotError * (nTrgs-1) ) + Error ) / nTrgs;
	              }
		      break;
		      
		 case NORM_SQM:
		      {
		         if ( nTrgs == 1 )
			   TotError = Error * Error;
			 else
			 {
		           TotError *= (nTrgs-1);
			   TotError += ( Error * Error );
			   TotError /= nTrgs;
			 }
	              }
		      break;
		      
		 default:
		    error( "MTTracer::UpdateError(), norm not implemented" );
		    		  
	      }
	   }
           break;
      
      case DELTRG:
           {
	      nTrgs--;
	      
	      switch( Norm )	      
	      {
	         case NORM_MAX:
		      {
		         Errors.Remove( Error );
			 if ( Errors.IsEmpty() )
			    TotError = 0.0;
			 else
			    TotError = Errors.GetMax();
		      }
		      break;
		      
		 case NORM_MED:
		      {
		         if ( nTrgs == 1 )
			    TotError = Error;
			 else			    
		            TotError = ( ( TotError * (nTrgs+1) ) - Error ) / nTrgs;
	              }
		      break;
		      
		 case NORM_SQM:
		      {
		         if ( nTrgs == 1 )
			    TotError = Error * Error;
			 else
			 {
		            TotError *= (nTrgs+1);
			    TotError -= ( Error * Error );
			    TotError  = TotError / nTrgs;
			 }
	              }
		      break;
		      
		 default:
		    error( "MTTracer::UpdateError(), norm not implemented" );
		    		  
	      }
	   }
           break;
	   
      default:
        error( "MTTracer::UpdateError(), inconsistency detected" );
   }
   
}
