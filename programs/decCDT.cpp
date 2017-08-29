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
//   file   : decCDT.cpp
//   author : Alessio Calcagno
//
//   Implementation of Class TDecCDT, sub-class of TDecimDelaunay, for the
//   decimation of a constrained Delaunay triangulation (CDT) read from
//   file.
//   It re-implements methods ReadData, WriteData, InitialTriangulation...
//

#include "defs.h"
#include "error.h"
#include "geom.h"
#include <stdlib.h>  // for exit in FindConstraintInVE
/*
#ifdef CC_VISUAL5        // for WriteConstraints
	#include <fstream.h>
#else
	#include <fstream>
#endif
*/
#include "ttriang.h"

#include "utils.h"
#include "decCDT.h"

// ----------------------------------------------------------------------------
// 
//  Constructor of class TDecCDT
//

// Constructor without special options.
TDecCDT::TDecCDT( int iK, PMTTracer iMT )
		: TDecimDelaunay( iK, iMT ), TDelaunayBase()
{

   #ifdef DEBUG
      DEBUG << "TDecCDT Constructor" << endl;
   #endif 
   
   nIncConstr = 0;
   ExtActive = TRUE;
   AllowFeaturesDel = FALSE;
   AllowChainBrk    = FALSE;

   nRemovedVertex = 0;
	Constraints = NULL;
	nConstrInFile = 0;
   nConstrInTRI = 0;
   ReadingConstr = FALSE;

}

// Constructor that allows specifying options for the extension of
// optimization, and for removing vertices incident in one or three
// constraints.
TDecCDT::TDecCDT( int iK, PMTTracer iMT, boolean EXTActive, boolean ALLOWFeaturesDel, boolean ALLOWChainBrk )
		: TDecimDelaunay( iK, iMT ), TDelaunayBase()
{

   #ifdef DEBUG
      DEBUG << "TDecCDT(Options) Constructor" << endl;
   #endif 
   
   nIncConstr = 0;
   ExtActive = EXTActive;
   AllowFeaturesDel = ALLOWFeaturesDel;
   AllowChainBrk    = ALLOWChainBrk;

   nRemovedVertex = 0;
	Constraints = NULL;
	nConstrInFile = 0;
   nConstrInTRI = 0;
   ReadingConstr = FALSE;

}


//
//  boolean TDecCDT::IsNewConstrOK( PTVertex V )
//
//  Decide whether a vertex having two incident constraints is removable.
//
//  INPUT: a vertex V having two incident constraints C1 and C2.
//
//  OUTPUT: true if the new constraint newConstr that should replace the
//  two existing constraints is admissible, i.e., it does not intersect
//  any other constraint, and it does not contain any vertex inside it.
//
//  SPECIAL CASES: if newConstr is an edge of the triangulation that is
//  not a constraint, then return true.
//  If newConstr is the same as an existing constraint C, the sequence of
//  edges < C C1 C2 > is closed, thus it is a polygon and we don't want to
//  eliminate it: thus return false (this choice is to be revised
//  according to the option that allows breaking closed chains).
//

boolean TDecCDT::IsNewConstrOK( PTVertex V ){

   #ifdef DEBUG3
      DEBUG3 << "\nTDecCDT::IsNewConstrOK( " << *V << " )\n";
   #endif 

   // find constraints incident in V

   TDoubleList<PTEdge> IncConstr;

   int numIncConstr = FindConstraintsInVE( V, IncConstr );
	if( numIncConstr < 2 )
	{  cerr << "BUG: TDecCDT::IsNewConstrOK( V ), " << *V << " has less than 2 incident edges Marked( CONSTRAINED )\n"; exit( -1); }
	if( numIncConstr > 2 )
	{  cerr << "BUG: TDecCDT::IsNewConstrOK( V ), " << *V << " has more than 2 incident edges Marked( CONSTRAINED )\n"; exit( -1); }

	PTEdge newConstr = new TEdge( Other_v_in_e( IncConstr.GetHead(), V ),
	                              Other_v_in_e( IncConstr.GetLast(), V )   );

   #ifdef CHIAMATE_DI_PROVA
      TDoubleList<PTEdge> Right_Border, Left_Border;
      if( CalcEdgInflRegn( newConstr, Right_Border, Left_Border ) )
      {
         CheckEdgeInflBorder( newConstr, Right_Border );
      //   Pause("Appena eseguito check del bordo dx");
         CheckEdgeInflBorder( newConstr, Left_Border );
      //   Pause("Appena eseguito check del bordo sx");
      }
   #endif

   PTVertex v0 = newConstr->EV[0], v1 = newConstr->EV[1];

   int intersection;
   // result of intersection of newConstr with edges of examined triangles

   //
   //  examine v0
   //

   // first and last triangle of VT( v0 )
   PTTriangle firstTrg = Next_CCW_t_of_e_by_v( FirstEdgeInVE( v0 ) , v0 );
   PTTriangle lastTrg = LastTrgInVE( v0 );

   PTTriangle currTrg;  // current triangle of VT( v0 )
   PTEdge oppositeEdg;  // edge of currTrg opposite to v0


   PTTriangle adj_v0_Trg, adj_v1_Trg;
   // adjacent triangles to v0/v1 and properly intersected by newConstr

   //  scan VT of v0 in counterclockwise order searching for a triangle
   //  Trg such that the edge opposite to v0 intersects E
   for( currTrg = firstTrg; /*currTrg != lastTrg*/; currTrg = Next_CCW_t_of_t_by_v( currTrg, v0 ) )
   {
      oppositeEdg = Opposite_e_of_t_by_v( currTrg, v0 );

      intersection = Tra_intersect( newConstr, oppositeEdg );

      if( intersection == PROPER_INTER )
      {
         if( oppositeEdg->Marked( CONSTRAINED ) )
            return FALSE;

         adj_v0_Trg = currTrg;
         break;
      }
      else
      if( intersection == ONLY_1CV )
      {
         // currTrg can be firstTrg: in this case the candidate edge
         // to be the same as newConstr can be the one preceding 
         // it in CCW order w.r.t. oppositeEdg.
         PTEdge nextEdg = Next_CCW_e_of_t_by_e( currTrg, oppositeEdg ),
                prevEdg = Prev_CCW_e_of_t_by_e( currTrg, oppositeEdg );

         // if edge newConstr already exists and is not a constraint,
         // then return true
         if( (  EqualsEV( newConstr, nextEdg ) && ! nextEdg->Marked( CONSTRAINED )  )
             ||
             (  EqualsEV( newConstr, prevEdg ) && ! prevEdg->Marked( CONSTRAINED )  )
                                                                                      )
         {
            return TRUE;
         }
         else
            return FALSE;
      }

      check( currTrg == lastTrg, "TDecCDT::IsNewConstrOK, triangle intersecting V0-V1 not found in VT( V0 )" );
   }

   //
   //  examine v1
   //
   firstTrg = Next_CCW_t_of_e_by_v( FirstEdgeInVE( v1 ) , v1 );
   lastTrg = LastTrgInVE( v1 );;

   for( currTrg = firstTrg; /*currTrg != lastTrg*/; currTrg = Next_CCW_t_of_t_by_v( currTrg, v1 ) )
   {
      oppositeEdg = Opposite_e_of_t_by_v( currTrg, v1 );

      intersection = Tra_intersect( newConstr, oppositeEdg );

      if( intersection == PROPER_INTER )
      {
         if( oppositeEdg->Marked( CONSTRAINED ) )
            return FALSE;

         adj_v1_Trg = currTrg;
         break;
      }
      else
      if( intersection == ONLY_1CV )  
         // if newConstr existed, we would have found it while examining v0
         return FALSE;

      check( currTrg == lastTrg, "TDecCDT::IsNewConstrOK, triangle intersecting V0-V1 not found in VT( V1 )" );
   }

   //
   // examine triangles intersected by newCostr between v0 and v1:
   // start from adj_v1_Trg and move through adjacency links until we
   // arrive at adj_v0_Trg
   //

   PTTriangle prevTrg;
   // adjacent triangle to the current triangle currTrg along edge
   // intersectingEdg p[roperly intersected by newConstr

   prevTrg = adj_v1_Trg;

   PTEdge intersectingEdg = oppositeEdg;

   currTrg = Opposite_t_of_t_by_e( currTrg, intersectingEdg );

 
   while( currTrg != adj_v0_Trg )
   {
      // search among edges of currTrg for the edge 
      // properly intersecting newConstr
      for( int i = 0; i < 3; i++ )
      {
         if( currTrg->TE[i] != intersectingEdg )
         {
            if( (intersection = Tra_intersect( newConstr , currTrg->TE[i] )) )
            {
               if( intersection == PROPER_INTER  &&  ! currTrg->TE[i]->Marked( CONSTRAINED ) )
               {
                  intersectingEdg = currTrg->TE[i];
                  currTrg = Opposite_t_of_t_by_e( currTrg, intersectingEdg );
                  break;
               }
               else return FALSE;
            }
         }
      }
   }
   return TRUE;
}


//
// Check that the boundary list is a closed sequence of edges.
//
void TDecCDT::CheckEdgeInflBorder( PTEdge E, TDoubleList<PTEdge> &border)
{
   #ifdef DEBUG7
      DEBUG7 << "\nTDecCDT::CheckEdgeInflBorder( " << *E << " ...): \n"; 
   #endif

   // search for vertex common to the last and the first edge
   // of border list
   PTEdge currEdg = border.GetLast();
   PTEdge nextEdg = border.GetHead();

   PTVertex currVtx; // common vertex to currEdg and nextEdg

   if( Is_v_in_e( currEdg->EV[0], nextEdg ) )
   { 
      if( Is_v_in_e( currEdg->EV[1], nextEdg ) )
         Pause( "TDecCDT::CheckEdgeInflBorder, inconsistency <a>" );

      currVtx = currEdg->EV[0];
      #ifdef DEBUG7
         DEBUG7 << *currEdg << " OK\n";
      #endif
   }
   else if( Is_v_in_e( currEdg->EV[1], nextEdg ) )
   {
      currVtx = currEdg->EV[1];
      #ifdef DEBUG7
         DEBUG7 << *currEdg << " OK\n";
      #endif
   }
   else
      Pause( "TDecCDT::CheckEdgeInflBorder, inconsistency <b>" );

   TDoubleListIterator<PTEdge> bordersIter( &border );
   bordersIter.Restart();

   currEdg = nextEdg;
   currVtx = Other_v_in_e(  currEdg, currVtx);
   bordersIter.GoNext();

   while( ! bordersIter.EndOfList() )
   {
      nextEdg = bordersIter.Current()->object;

      if( Is_v_in_e( currVtx, nextEdg ) )
      {
         #ifdef DEBUG7
            DEBUG7 << *currEdg << " OK\n";
         #endif
      }
      else
      {
         cerr << "ERROR: TDecCDT::CheckEdgeInflBorder( " << *E << " ...): \n\t"
              << *currEdg << " is not followed by " << *nextEdg << "\n";
         Pause( " " );
      }
      currEdg = nextEdg;
      currVtx = Other_v_in_e(  currEdg, currVtx);
      bordersIter.GoNext();
   }
}


void TDecCDT::CalcEdgeWideInflRegn( PTEdge E,  TDoubleList<PTTriangle> &InternalTrgs,
                                TDoubleList<PTEdge> &Right_border, TDoubleList<PTEdge> &Left_border )
{   
   if( InitCalcEdgeWideInflRegn( E, InternalTrgs, Right_border, Left_border ) == EDGE_ALREADY_EXISTS )
      return;


   PTTriangle currTrg, oppTrg;
   // current triangle and triangle opposite to the given edge of currTrg

   int i;

   while( ! InternalTrgs.IsEmpty() )
   {
      currTrg = InternalTrgs.RemoveHead();

      for( i = 0 ; i < 3; i++ )  // scan edges of currTrg
      {
         if( ! Tra_intersect( E, currTrg->TE[i] ) )
         {
            currTrg->TE[i]->Mark( INFL_BORDER );
            if( Geom::Turnxy( E->EV[0], E->EV[1], currTrg->TE[i]->EV[0] ) == TURN_RIGHT )
               Right_border.AddTail( currTrg->TE[i] );
            else
               Left_border.AddHead( currTrg->TE[i] );
         }
         else
         {
            oppTrg = Opposite_t_of_t_by_e( currTrg, currTrg->TE[i] );
            if( ! oppTrg->Marked( VISITED )  )
            {
               oppTrg->Mark( VISITED );
               InternalTrgs.AddTail( oppTrg );
            }
         }
      }
   }
   PTEdge tmpEdg;  
   // adjust order of edges inserted by InitCalcEdgeWideInflRegn 
   // during examination of E->EV[1].

   tmpEdg = Right_border.RemoveHead();
   Right_border.AddTail( tmpEdg );

   tmpEdg = Left_border.RemoveLast();
   Left_border.AddHead( tmpEdg );

   CheckEdgeInflBorder( E, Right_border );
   CheckEdgeInflBorder( E, Left_border );   
}


int TDecCDT::InitCalcEdgeWideInflRegn( PTEdge E, TDoubleList<PTTriangle> &InternalTrgs,
                    TDoubleList<PTEdge> &Right_border, TDoubleList<PTEdge> &Left_border )
{
   PTVertex v0 = E->EV[0], v1 = E->EV[1];

   //
   //  examine v0
   //

   // primo ed ultimo triangolo di VT( v0 )
   PTTriangle firstTrg = Next_CCW_t_of_e_by_v( FirstEdgeInVE( v0 ) , v0 );
   PTTriangle lastTrg = LastTrgInVE( v0 );

   PTTriangle currTrg;  // current trianglr of VT( v0 )
   PTEdge oppositeEdg;  // edge of currTrg opposite to v0

   PTEdge EdgR, EdgL;  
   // edge of Right_border and Left_border intersecting v0/v1 

   PTTriangle adjTrg;
   // triangle adjacent to the triangle currTrg currently scanned

   //  scan VT of v0 in counterclockwise order searching from a triangle
   //  Trg such that the edge opposite to v0 intersects E
   for( currTrg = firstTrg; currTrg != lastTrg; currTrg = Next_CCW_t_of_t_by_v( currTrg, v0 ) )
   {
      oppositeEdg = Opposite_e_of_t_by_v( currTrg, v0 );

      int intersection = Tra_intersect( E, oppositeEdg );

      if( intersection == PROPER_INTER )
      {
         EdgR = Prev_CCW_e_of_t_by_e( currTrg, oppositeEdg );
         EdgL = Next_CCW_e_of_t_by_e( currTrg, oppositeEdg );

         Right_border.AddTail( EdgR );   EdgR->Mark( INFL_BORDER );
         Left_border.AddHead( EdgL );    EdgL->Mark( INFL_BORDER );

         currTrg->Mark( VISITED );

         // Enqueue the triangle opposite to currTrg w.r.t. oppositeEdg.
         // If such Trg does not contain v1 then do not enqueue it because
         // it will be marked and processed during the scan of v1.
         adjTrg = Opposite_t_of_t_by_e( currTrg, oppositeEdg );
         if( ! Is_v_in_t( v1, adjTrg ) )
         {
            adjTrg->Mark( VISITED );
            InternalTrgs.AddTail( adjTrg );
         }
         break;
      }
      else
      if( intersection == ONLY_1CV )
      {
         if( Is_v_in_e( v1, oppositeEdg ) ) 
         // the edge E that we are going to insert already exists
         {
            return EDGE_ALREADY_EXISTS;
         }
         
         EdgR = Prev_CCW_e_of_t_by_e( currTrg, oppositeEdg );

         PTTriangle currTrg2 = Next_CCW_t_of_t_by_v( currTrg, v0 );

         EdgL = Left_e_of_t_by_v( currTrg2, v0 );

         Right_border.AddTail( EdgR );   EdgR->Mark( INFL_BORDER );
         Left_border.AddHead( EdgL );    EdgL->Mark( INFL_BORDER );

         currTrg->Mark( VISITED );
         currTrg2->Mark( VISITED );

         // insert opposite triangles as in case PROPER_INTER
         adjTrg = Opposite_t_of_t_by_e( currTrg, oppositeEdg );
         if( ! Is_v_in_t( v1, adjTrg ) )
         {
            adjTrg->Mark( VISITED );
            InternalTrgs.AddTail( adjTrg );
         }

         adjTrg = Opposite_t_of_t_by_e( currTrg2, Opposite_e_of_t_by_v( currTrg2, v0 ) );
         if( ! adjTrg->Marked( VISITED )  &&  ! Is_v_in_t( v1, adjTrg ) ) 
         {  // if triangle opposite to currTrg2 is NOT the same
            // as the triangle opposite to currTrg
            adjTrg->Mark( VISITED );
            InternalTrgs.AddTail( adjTrg );
         }
         break;
      }
   }

   //
   //  examine v1
   //
   firstTrg = Next_CCW_t_of_e_by_v( FirstEdgeInVE( v1 ) , v1 );
   lastTrg = LastTrgInVE( v1 );;

   //  Scan the VT relation of v1 in counterclockwise order searching for
   //  a triangle Trg such that the edge of Trg opposite to v1 intersects E.
   //  Unlike for v0, we put no edge in InternalTrgs:
   //  simply we mark the two boundary edges and the (one or two)
   //  triangle(s) incident in v1 and intersecting E.
   //  The idea is that, during the traversal of the influence region of E
   //  se move from a generic triangle t to those triangles adjacent to
   //  edges NOT Marked( INFL_BORDER ) and not marked VISITED; we mark as
   //  visited al, the edges of t that do not intersect E.
   //  This works for any triangle internal to the region of influence,
   //  with the exclusion of those incident in a vertex of E, because
   //  their boundary edges intersect E.
   //  Therefore, we visit such triangle in a preliminary step, we mark
   //  then as VISITED and mark their boundary edge(s) as INFL_BORDER.

   for( currTrg = firstTrg; currTrg != lastTrg; currTrg = Next_CCW_t_of_t_by_v( currTrg, v1 ) )
   {
      oppositeEdg = Opposite_e_of_t_by_v( currTrg, v1 );

      int intersection = Tra_intersect( E, oppositeEdg );

      if( intersection == PROPER_INTER )
      {
         EdgR = Prev_CCW_e_of_t_by_e( currTrg, oppositeEdg );
         EdgL = Next_CCW_e_of_t_by_e( currTrg, oppositeEdg );

         Right_border.AddHead( EdgR );
         // put it here by now, but at the end of CalcEdgeWideInflRegn
         // we must put it back to its place
         EdgR->Mark( INFL_BORDER );
         Left_border.AddTail( EdgL );  
         // put it here by now, but at the end of CalcEdgeWideInflRegn
         // we must put it back to its place
         EdgL->Mark( INFL_BORDER );

         currTrg->Mark( VISITED );

         break;
      }
      else
      if( intersection == ONLY_1CV )
      {         
         EdgR = Prev_CCW_e_of_t_by_e( currTrg, oppositeEdg );

         PTTriangle currTrg2 = Next_CCW_t_of_t_by_v( currTrg, v1 );

         EdgL = Left_e_of_t_by_v( currTrg2, v1 );

         Right_border.AddHead( EdgR );
         // put it here by now, but at the end of CalcEdgeWideInflRegn
         // we must put it back to its place
         EdgR->Mark( INFL_BORDER );
         Left_border.AddTail( EdgL );
         // put it here by now, but at the end of CalcEdgeWideInflRegn
         // we must put it back to its place
         EdgL->Mark( INFL_BORDER );

         currTrg->Mark( VISITED );
         currTrg2->Mark( VISITED );

         break;
      }
   }
   return 0;
}

//
//  Return values: see #define's at the beginning of decCDT.h
//  

int TDecCDT::Tra_intersect( PTEdge e1, PTEdge e2 )
{
   check( e1 == NULL || e2 == NULL || ! Is_EV_Defined( e1 ) || ! Is_EV_Defined( e2 ),
          "TDecCDT::intersect, called on NULL edge" );
   check( e1->EV[0] == e1->EV[1] || e2->EV[0] == e2->EV[1], "TDecCDT::intersect, called on degenerate edge" );

  // i due edge sono coincidenti ( hanno entrambi i vertici in comune )
   if (( e1->EV[0] == e2->EV[0] &&  e1->EV[1] == e2->EV[1] ) ||
       ( e1->EV[0] == e2->EV[1] &&  e1->EV[1] == e2->EV[0] )  )
         return UPON_AND_2CV;
          
  // if we arrive here, then e1 and e2  are NOT the same

  int turn_e1_first_e2  = Geom::Turnxy( e1->EV[0], e1->EV[1], e2->EV[0] );
  // turn from e1 to first vertex of e2
  int turn_e1_second_e2 = Geom::Turnxy( e1->EV[0], e1->EV[1], e2->EV[1] );
  // turn from e1 to second vertex of e2
  int turn_e2_first_e1  = Geom::Turnxy( e2->EV[0], e2->EV[1], e1->EV[0] );
  // turn from e2 to first vertex of e1
  int turn_e2_second_e1 = Geom::Turnxy( e2->EV[0], e2->EV[1], e1->EV[1] );
  // turn from e2 to second vertex of e1

  // e1 and e2 are collinear (but they may not intersect)
  if( turn_e1_first_e2 == turn_e1_second_e2  &&  turn_e1_first_e2 == ALIGNED  &&
      turn_e2_first_e1 == turn_e2_second_e1  &&  turn_e2_first_e1 == ALIGNED )
      // if all == ALIGNED
  {
   //   check( turn_e1_first_e2 != ALIGNED, "TDecCDT::Tra_intersect, BUG: logical error" );
      
      // the four vertices are aligned but they may not intersect
      //  
      if( e1->EV[0]->x == e1->EV[1]->x )
      // aligned along y axis, use y coordinate to compare them
      {
         double max_y_e1 = MAX( e1->EV[0]->y , e1->EV[1]->y );
         double min_y_e1 = MIN( e1->EV[0]->y , e1->EV[1]->y );
         double max_y_e2 = MAX( e2->EV[0]->y , e2->EV[1]->y );
         double min_y_e2 = MIN( e2->EV[0]->y , e2->EV[1]->y );

         if( max_y_e1 < min_y_e2 )  return NO_INTER;
         if( max_y_e1 == min_y_e2 || max_y_e2 == min_y_e1 )  return ONLY_1CV;
         if( min_y_e1 == min_y_e2 || max_y_e1 == max_y_e2 )  return UPON_AND_1CV;
         return  UPON_NO_CV;
      }
      else // NOT aligned along y axis, use x coordinate to compare them
      {
         double max_x_e1 = MAX( e1->EV[0]->x , e1->EV[1]->x );
         double min_x_e1 = MIN( e1->EV[0]->x , e1->EV[1]->x );
         double max_x_e2 = MAX( e2->EV[0]->x , e2->EV[1]->x );
         double min_x_e2 = MIN( e2->EV[0]->x , e2->EV[1]->x );

         if( max_x_e1 < min_x_e2 )  return NO_INTER;
         if( max_x_e1 == min_x_e2 || max_x_e2 == min_x_e1 )  return ONLY_1CV;
         if( min_x_e1 == min_x_e2 || max_x_e1 == max_x_e2 )  return UPON_AND_1CV;
         return  UPON_NO_CV;
      }

  }

  //
  // e1 and e2 are NOT collinear
  // 

  //  e1 and e2 have only one common vertex (and are not collinear)
   if( (e1->EV[0] == e2->EV[0]) ||   (e1->EV[1] == e2->EV[1]) ||
       (e1->EV[0] == e2->EV[1]) ||   (e1->EV[1] == e2->EV[0])    )  return ONLY_1CV;

  // vertex of one edge is inside the other edge
  // 
  if( ( turn_e1_first_e2  == ALIGNED  || turn_e1_second_e2 == ALIGNED ) && turn_e2_first_e1  != turn_e2_second_e1 )
     return NOT_PROPER;
  if( ( turn_e2_first_e1  == ALIGNED  || turn_e2_second_e1 == ALIGNED ) && turn_e1_first_e2  != turn_e1_second_e2 )
     return NOT_PROPER;

  // now either intersection is proper (interior of e1 intersects 
  // interior of e2), or there is no intersection
  if( turn_e1_first_e2  != turn_e1_second_e2 && turn_e2_first_e1  != turn_e2_second_e1 )
     return PROPER_INTER;
  else
     return NO_INTER;
}

// Used inside function Tra_intersect
void TDecCDT::Tra_intersect_test()
{
      TVertex a(0,0,0), b(2,0,0), c(2,2,0), d(0,2,0), e(1,1,0), f(3,-1,0), g(0,-1,0),
               h(0,-2,0), i(4,0,0), j(5,0,0);
      TEdge A(&a, &e), B(&b, &c), C(&a, &d),   D(&a, &c), E(&d, &e), F(&d, &b), G(&e, &f),
            H(&g,&h), I(&i,&j), J(&a,&b), K(&h,&a), L(&g,&d), M(&j,&b), N(&i,&a), O(&b,&f) ;
      cerr << "Tra_intersect( A, B ) = " << Tra_intersect( &A, &B ) << endl;
      cerr << "Tra_intersect( H, C ) = " << Tra_intersect( &H, &C ) << endl;
      cerr << "Tra_intersect( J, I ) = " << Tra_intersect( &J, &I ) << endl;
      cerr << "Tra_intersect( E, O ) = " << Tra_intersect( &E, &O ) << endl << endl;

      cerr << "Tra_intersect( A, C ) = " << Tra_intersect( &A, &C ) << endl;
      cerr << "Tra_intersect( E, G ) = " << Tra_intersect( &E, &G ) << endl;
      cerr << "Tra_intersect( L, H ) = " << Tra_intersect( &L, &H ) << endl;
      cerr << "Tra_intersect( J, M ) = " << Tra_intersect( &J, &M ) << endl << endl;

      cerr << "Tra_intersect( A, F ) = " << Tra_intersect( &A, &F ) << endl;
      cerr << "Tra_intersect( J, L ) = " << Tra_intersect( &J, &L ) << endl;
      cerr << "Tra_intersect( N, B ) = " << Tra_intersect( &N, &B ) << endl << endl;

      cerr << "Tra_intersect( D, F ) = " << Tra_intersect( &D, &F ) << endl;
      cerr << "Tra_intersect( G, N ) = " << Tra_intersect( &G, &N ) << endl << endl;

      cerr << "Tra_intersect( F, G ) = " << Tra_intersect( &F, &G ) << endl;
      cerr << "Tra_intersect( L, K ) = " << Tra_intersect( &L, &K ) << endl;
      cerr << "Tra_intersect( M, N ) = " << Tra_intersect( &M, &N ) << endl << endl;

      cerr << "Tra_intersect( A, D ) = " << Tra_intersect( &A, &D ) << endl;
      cerr << "Tra_intersect( H, K ) = " << Tra_intersect( &H, &K ) << endl;
      cerr << "Tra_intersect( M, I ) = " << Tra_intersect( &M, &I ) << endl << endl;

      cerr << "Tra_intersect( A, A ) = " << Tra_intersect( &A, &A ) << endl;
}


void TDecCDT::DetachEdge( PTEdge E )
{   
   #ifdef DEBUG4
      DEBUG4 << "TDecCDT::DetachEdge( " << *E << " )\n";
   #endif 
   
   if( E->Marked( CONSTRAINED ) )
   {
      nConstrInTRI--;
      check(nConstrInTRI < 0, "TDecCDT::DetachEdge < nConstrInTRI < 0 >");
      for (int v=0; v<2; v++ )
      {
         E->EV[v]->nIncConstr--;
         check(E->EV[v]->nIncConstr < 0, "TDecCDT::DetachEdge < nIncConstr < 0 >");
      }
		E->UnMark(CONSTRAINED);  // added again
		#ifdef array_Constraints
		   for( int c = 0; c < nConstrInFile; c++ )
		   	if( Constraints[c] == E )
		   	   { Constraints[c] = NULL; break; }
		#endif
   }
   TDecimDelaunay::DetachEdge(E);
}

//
//  Re-implemented to insert the new constraint in case VertexToRemove
//  is incident in two constraints before computing and re-triangulating 
//  the region of influence of VertexToRemove.
//
void TDecCDT::RemoveVertex()
{

   #ifdef DEBUG12
      DEBUG12 << "\nTDecCDT::RemoveVertex( " << *VertexToRemove << " )\n";
   #endif 

   nIncConstr = VertexToRemove->nIncConstr;
   // number of constraints incident in VertexToRemove

   #ifdef ROBUST
      check( nIncConstr < 0 || nIncConstr > 2, "TDecCDT::RemoveVertex(), called on a vertex with wrong constrDegree" );
   #endif

   if( nIncConstr == 0 )  { TDecimDelaunay::RemoveVertex(); }
   else
   {  
      PTEdge OldConstr1 = NULL,  // incident constraints in VertexToRemove
             OldConstr2 = NULL;

      TDoubleList<PTEdge> IncConstrL;
      // list of constraint edges incident in VertexToRemove

      if( nIncConstr == 2 ) 
      // first add new constraint, then re-triangulate region of influence
      // of VertexToRemove.
      {
         PTVertex NewConstrVtx[2];      // vertices of new constraint
         FindConstraintsInVE( VertexToRemove, IncConstrL );

		   #ifdef ROBUST
			   check( IncConstrL.Lenght() != 2, "TDecCDT::RemoveVertex(), <a> " );
		   #endif

         OldConstr1 = IncConstrL.RemoveHead();
         OldConstr2 = IncConstrL.RemoveLast();

         NewConstrVtx[0] = Other_v_in_e( OldConstr1, VertexToRemove );
         NewConstrVtx[1] = Other_v_in_e( OldConstr2, VertexToRemove );
         // Add_Constraint( NewConstrVtx[0], NewConstrVtx[1] );
         // Tra_Add_Constraint( NewConstrVtx[0], NewConstrVtx[1] );
         TTra_Add_Constraint( NewConstrVtx[0], NewConstrVtx[1] );
      }
      // TDecimDelaunay::RemoveVertex(); 
      // by now call original version that does not extend optimization
      // (TDecCDT::)RIRPrepare( SwapEdgeQueue );                     or  (TDecCDT::)ExtRIRPrepare( SwapEdgeQueue );
      // (TDecCDT::)RIRInitialTrg( InflRegnBorder, SwapEdgeQueue );  or  (TDecCDT::)ExtRIRInitialTrg( InflRegnBorder, SwapEdgeQueue );
      // (TDecCDT::)RIRDelOptTrg( SwapEdgeQueue );                   or  (TDecCDT::)ExtRIRDelOptTrg( SwapEdgeQueue );
      // (TDecCDT::) DeleteInfluenceRegion ( InflRegnBorder );       or  (TDecCDT::)ExtDeleteInfluenceRegion ( InflRegnBorder );

      if( ! ExtActive )
         TDecimDelaunay::RemoveVertex();
      else
      {
         if( nIncConstr == 1 )
         {
            FindConstraintsInVE( VertexToRemove, IncConstrL );

		      #ifdef ROBUST
			      check( IncConstrL.Lenght() != 1, "TDecCDT::RemoveVertex(), <b> " );
		      #endif

            OldConstr1 = IncConstrL.RemoveHead();
         }

         // Version with extension of optimization beyond the region
         // of influence of VertexToRemove
         //
         boolean closed = CalcInfluenceRegion();
         ExtRIRIPrepare( SwapEdgeQueue, OldConstr1, OldConstr2, closed, OrigEdgList );
         ExtRIRInitialTrg( InflRegnBorder, SwapEdgeQueue );
         ExtRIRDelOptTrg( SwapEdgeQueue, OrigEdgList );
         ExtDeleteInfluenceRegion ( InflRegnBorder );
      }
   }

   RemovedVertexIndex.AddTail( VertexToRemove->VID );
   // VID of removed vertices
   nRemovedVertex++;
}




//
// Re-implemented to read the constraints.
//
// Call ReadConstraints that marks CONSTRAINED the constraint edges 
// and fills array Constraints with pointers to them.

void TDecCDT::ReadData( const char *infname )
{

   #ifdef DEBUG3
    DEBUG3 << "\nTDecCDT::ReadData()" << endl;
   #endif
    
   ifstream inFile;
   
   inFile.open( infname );
   
   check( !inFile, "TDecCDT::ReadData(), cannot open input file");

	TDecimDelaunay::ReadVertices( inFile );
   
   TDecimDelaunay::ReadTriangles( inFile );

	ReadConstraints( inFile );

	inFile.close();
}

//
// Read constraints from inout file.
//
// For each constraint:
//   - mark it as CONSTRAINED
//   - increment by one variable nIncConstr of the two endpoint vertices
//     of the constraint
// Accept also files with degenerate constraints (with two equal endpoints)
// or constraints that do not correspond to any edge of the input
// triangulation.
// If macro array_Constraints is defined, then allocate and fills array
// Constraints with pointers to the constraint edges.
//

void TDecCDT::ReadConstraints( ifstream & inFile ) {

   ReadingConstr = TRUE;
	
	int i;

   #ifdef DEBUG3
       DEBUG3 << "\nTDecCDT::ReadConstraints()" << endl;
   #endif
   
   check( (inFile.eof()), "TDecCDT::ReadConstraints(), unexpected End Of File");
   
   //
   // read number of constraints
   //

   inFile >> nConstrInFile;

   check( (nConstrInFile < 0), "input with less than 0 constraints" );

   #ifdef array_Constraints   
      Constraints = new PTEdge[nConstrInFile];
   
      check( (Constraints == NULL), "TDecCDT::ReadConstraints(), insufficient memory" );
   #endif
   
   //
   // read constraints, and put pointers to them in array Constraints
   //

   PTEdge TmpE = new TEdge( NULL, NULL );  // current constraint

   check( ( TmpE == NULL), "TDecCDT::ReadConstraints(), insufficient memory" );

   PTEdge matchingE;  // pointer to corresponding edge 

   int vidx[2]; // indices (VID) of vertices of this constraint

   nConstrInTRI = nConstrInFile;

   int nNotExistingConstr = 0;  
   // number of input  constraints that do not correspond to
   // any edge of the triangulation

   int nDuplicatedConstr = 0;   // number of duplicated constraints
   int nDegenerateConstr = 0;   // number of degenerate constraints

   for( i=0; i<nConstrInFile; i++ )
   {
      
      check( (inFile.eof()), "TDecCDT::ReadConstraints(), unexpected End Of File");
  
      //
      // read edge (indices of its two vertices)
      //
       
      inFile >> vidx[0] >> vidx[1];

      // the two vertices may teh the same (degenerate constraint)

      if( vidx[0] == vidx[1] )
      {
         nConstrInTRI--;
         check(nConstrInTRI < 0, "TDecCDT::FindConstraintInVE: nConstrInTRI < 0  while reading input file!?!?!");
         nDegenerateConstr++;

         #ifdef array_Constraints
            Constraints[i] = NULL;
         #endif

         if ( nDegenerateConstr < 7 )
            cerr << "\nWARNING: constraint C" << i << "( "
                 << vidx[0] << ", " << vidx[1] << " ) is degenerate\n";
         else if ( nDegenerateConstr == 7 )
            cerr << "\nWARNING: other degenerate constraints read\n";
         continue;
      }

		TmpE->EV[0] = (PTVertex) (Points[vidx[0]]);
      TmpE->EV[1] = (PTVertex) (Points[vidx[1]]);

      matchingE = FindConstraintInVE( (PTVertex) (Points[vidx[0]]) , TmpE );

      #ifdef array_Constraints
         Constraints[i] = matchingE;
      #endif

      if( matchingE == NULL )
      {
         nConstrInTRI--;
         check(nConstrInTRI < 0, "TDecCDT::FindConstraintInVE: nConstrInTRI < 0  while reading input file!?!?!");
         nNotExistingConstr++;
      }
      else
	   {
         if( matchingE->Marked( CONSTRAINED ) )
         {  
            nDuplicatedConstr++;
            nConstrInTRI--;
            #ifdef array_Constraints
               Constraints[i] = NULL;
            #endif
            check(nConstrInTRI < 0, "TDecCDT::FindConstraintInVE: nConstrInTRI < 0  while reading input file!?!?!");
         }
         else
         {
   	    matchingE->Mark( CONSTRAINED );
            // increment number of constraint incident in the two endpoints
            // of the constraint
            matchingE->EV[0]->nIncConstr++;  
	    matchingE->EV[1]->nIncConstr++;  
         }
     	}    
		 

       #ifdef DEBUG4
        DEBUG4 << "TDecCDT::ReadConstraints: read constraint " << *Constraints[i] << endl;
       #endif
              
       if ( i % 1000 == 0 ) cerr << "\rread " << i << " constraints" << flush;
       
   }
   
    cerr << "\rread " << nConstrInFile << " constraints:\n      " <<
            nNotExistingConstr << " do not match with any edge\n      " <<
            nDuplicatedConstr  << " are duplicated\n      " <<
            nDegenerateConstr  << " are degenerate\n      " <<
            nConstrInTRI << " are present in the triangulation\n";

   ReadingConstr = FALSE;
}


//
// Search, among the edges incident in V, the ones that match with
// filter, put them in a list and return their number. 
// "Matching" is implemented by function filterFun passed as a 
// parameter.
//

int TDecCDT::FindEdgesInVE( PTVertex V, PTEdge filter, TDoubleList<PTEdge> & filtered, PF1b2PTEdge filterFun)
{

   #ifdef DEBUG3
     DEBUG3 << "\nTDecCDT::FindEdgesInVE( " << *V << " , " << *filter << " )\n";
   #endif
   
	#ifdef ROBUST
		check( ! filtered.IsEmpty() , "TDecCDT::FindEdgesInVE, called on non empty list!" );
	#endif

   int i, j, nMatches = 0;

   PTTriangle TFirst;
   
   PTEdge   EFirst = V->VE[0];
   PTVertex VFirst = ( EFirst->EV[0] != V ? EFirst->EV[0] : EFirst->EV[1] );
   

   if ( EFirst->OnConvexHull() )
   {
       
       //
       // Must decide if it is right to start from EFirst to traverse
       // moving in counterclockwise order, all the region of influence
       // of V. In case starting from EFirst in counterclockwise order
       // we get out of the convex hull, then we must start from the other
       // edge adjacent to V, V->VE[1]
       //  
       
       TFirst = ( EFirst->ET[0] != NULL ? EFirst->ET[0] : EFirst->ET[1] );
  
       #ifdef ROBUST
          check( (TFirst == NULL), 
	        "TDecCDT::FindEdgesInVE(), <0> inconsistency detected");
       #endif
            
       PTVertex v[3];
       TFirst->GetTV( v[0], v[1], v[2] );

       
       // search for vertex of TFirst different from
       
       for( i=0; i<3; i++ )
          if ( v[i] != V && v[i] != VFirst ) break;
	 
	 
       #ifdef ROBUST
          check( (i>=3), "TDecCDT::FindEdgesInVE(), <1> inconsistency detected" );
       #endif
	  
	  
       if ( Geom::Turnxy( V, VFirst, v[i] ) != TURN_LEFT )
          EFirst = V->VE[1];
      
   } // end ...if( EFirst->OnConvexHull() )
   

   
   // triangle after EFirst, in counterclockwise order
   
   TFirst = ( EFirst->EV[0] == V ? EFirst->ET[0] : EFirst->ET[1] );
   
   #ifdef ROBUST
      check( ( TFirst==NULL), "TDecCDT::FindEdgesInVE(), <2> inconsistency detected" );
   #endif
   
   
   PTEdge ENext = EFirst;
   PTTriangle TNext;

   do   
   {
      if( filterFun( ENext,filter ) )
      // matching found between edges incident in V
      {
	 filtered.AddTail( ENext );
         nMatches++;
      }
      
      // go to the next one
      
      TNext = ( ENext->EV[0] == V ? ENext->ET[0] : ENext->ET[1] );
      
      if ( TNext != NULL )
      {
         // search for index of E in relation TE of TNext
	 
	 for( j=0; j<3; j++ )
	   if ( TNext->TE[j] == ENext ) break;
	   
	 #ifdef ROBUST
	    check( (j>=3), "TDecCDT::FindConstraintInVE(), <3> inconsistency detected" );
	 #endif

         ENext = TNext->TE[(j+2)%3]; // edge preceding E in TE of TNext
	 
      }
           
   } while ( TNext != NULL && ENext != EFirst );

   #ifdef ROBUST
      check( nMatches != filtered.Lenght(), "TDecCDT::FindConstraintInVE(), BUG nella funzione");
   #endif

   return nMatches;
}


//
// Search, among the edges incident in V, the ones that match with
// filter, put them in a list and return their number. 
// "Matching" is implemented by function TEdge::Match. 
//

int TDecCDT::FindEdgesInVE( PTVertex V, PTEdge filter, TDoubleList<PTEdge> & filtered) {


   #ifdef DEBUG3
     DEBUG3 << "\nTDecCDT::FindEdgesInVE( " << *V << " , " << *filter << " )\n";
   #endif
   
	#ifdef ROBUST
		check( ! filtered.IsEmpty() , "TDecCDT::FindEdgesInVE, called on non empty list!" );
	#endif

   int i, j, nMatches = 0;

   PTTriangle TFirst;
   
   PTEdge   EFirst = V->VE[0];
   PTVertex VFirst = ( EFirst->EV[0] != V ? EFirst->EV[0] : EFirst->EV[1] );

   if ( EFirst->OnConvexHull() )
   {
       
       //
       // Must decide if it is right to start from EFirst to traverse
       // moving in counterclockwise order, all the region of influence
       // of V. In case starting from EFirst in counterclockwise order
       // we get out of the convex hull, then we must start from the other
       // edge adjacent to V, V->VE[1]
       //  
       
       TFirst = ( EFirst->ET[0] != NULL ? EFirst->ET[0] : EFirst->ET[1] );
  
       #ifdef ROBUST
          check( (TFirst == NULL), 
	        "TDecCDT::FindEdgesInVE(), <0> inconsistency detected");
       #endif
            
       PTVertex v[3];
       TFirst->GetTV( v[0], v[1], v[2] );

       
       // search for verticex of TFirst different from
       
       for( i=0; i<3; i++ )
          if ( v[i] != V && v[i] != VFirst ) break;
	 
	 
       #ifdef ROBUST
          check( (i>=3), "TDecCDT::FindEdgesInVE(), <1> inconsistency detected" );
       #endif
	  
	  
       if ( Geom::Turnxy( V, VFirst, v[i] ) != TURN_LEFT )
          EFirst = V->VE[1];
      
   } // end ...if( EFirst->OnConvexHull() )
   

   
   // triangle after EFirst, in counterclockwise order
   
   TFirst = ( EFirst->EV[0] == V ? EFirst->ET[0] : EFirst->ET[1] );
   
   #ifdef ROBUST
      check( ( TFirst==NULL), "TDecCDT::FindEdgesInVE(), <2> inconsistency detected" );
   #endif
   
   
   PTEdge ENext = EFirst;
   PTTriangle TNext;

   do   
   {
      if( ENext->Match( filter ) )
      // matching found among the edges incident in V
      {
	 filtered.AddTail( ENext );
         nMatches++;
      }
      
      // go to tye next one
      
      TNext = ( ENext->EV[0] == V ? ENext->ET[0] : ENext->ET[1] );
      
      if ( TNext != NULL )
      {
         // search for index of E in TE relation of TNext
	 
	 for( j=0; j<3; j++ )
	   if ( TNext->TE[j] == ENext ) break;
	   
	 #ifdef ROBUST
	    check( (j>=3), "TDecCDT::FindConstraintInVE(), <3> inconsistency detected" );
	 #endif

         ENext = TNext->TE[(j+2)%3]; // edge preceding E in TE of TNext
	 
      }
           
   } while ( TNext != NULL && ENext != EFirst );

   #ifdef ROBUST
      check( nMatches != filtered.Lenght(), "TDecCDT::FindConstraintInVE(), BUG nella funzione");
   #endif

   return nMatches;
}



//
// Search among the edges incident in V the ones that are constraints,
// return in IncConstr a list of such edges and return their number.
//

int TDecCDT::FindConstraintsInVE( PTVertex V, TDoubleList<PTEdge> & IncConstr) {

   #ifdef DEBUG3
     DEBUG3 << "\nTDecCDT::FindConstraintsInVE( V" << V->VID << " )"<< endl;
   #endif
   
	#ifdef ROBUST
		check( ! IncConstr.IsEmpty() , "TDecCDT::FindConstraintsInVE, <a> inconsistency detected!" );
	#endif

   PTEdge matchConstr = new TEdge();
   matchConstr->Mark( CONSTRAINED );
   return FindEdgesInVE( V, matchConstr, IncConstr);
   // return FindEdgesInVE( V, matchConstr, IncConstr, EqualsMarkConstrained);
}

//
// Search among the edges incident in V one edge formed by the same
// vertices as Edge2find; if found, return a pointer to such edge,
// otherwise return NULL.
//

PTEdge TDecCDT::FindEdgeInVE( PTVertex V, PTEdge Edge2find ) {

   #ifdef DEBUG3
     DEBUG3 << "\nTDecCDT::FindEdgeInVE( V" << V->VID << ", " << *Edge2find << " )" << endl;
   #endif

   TDoubleList<PTEdge> MatchList;
   int nMatches = FindEdgesInVE( V, Edge2find, MatchList );
//   int nMatches = FindEdgesInVE( V, Edge2find, MatchList, EqualsEV );
   if( nMatches > 1 )
      error("TDecCDT::FindEdgeInVE, inconsistency: more than one matching edge!\n");
   
   if( nMatches == 1 )
   {
      PTEdge MatchingEdge = MatchList.GetHead();
      MatchingEdge->Mark( CONSTRAINED );
      return MatchingEdge;
   }

   if( nMatches == 0 )
      return NULL;
   else
      error("TDecCDT::FindConstraintInVE, BUG: error in function\n");
      
   return NULL; // must return something
}


//
// Search among the edges incident in V one edge formed by the same
// vertices as constr; if found, return a pointer to such edge,
// otherwise return NULL.
// Similar to TDecErrDelaunay::RecalVertexErrorApprox()
//

PTEdge TDecCDT::FindConstraintInVE( PTVertex V, PTEdge constr )
{

   #ifdef DEBUG3
     DEBUG3 << "\nTDecCDT::FindConstraintInVE( V" << V->VID << ", " << *constr << " )" << endl;
   #endif

   TDoubleList<PTEdge> MatchList;
   int nMatches = FindEdgesInVE( V, constr, MatchList );
//   int nMatches = FindEdgesInVE( V, constr, MatchList, EqualsEV );
   if( nMatches > 1 )
      error("TDecCDT::FindConstraintInVE, inconsistency: more than one matching edge!\n");
   
   if( nMatches == 1 )
   {
      PTEdge MatchingEdge = MatchList.GetHead();
      return MatchingEdge;
   }

   // if we arrive here then we have not found "constr" among the edges
   // incident in V, thus such constraint is not present in the
   // triangulation
   if( nMatches == 0  && ReadingConstr )
   {
      if ( (nConstrInFile - nConstrInTRI) < 7 )
         cerr << "\nWARNING: constraint " << constr->EV[0]->VID << " "
              << constr->EV[1]->VID << " does not match with any edge\n";
      else if ( (nConstrInFile - nConstrInTRI) == 7 )
         cerr << "\nWARNING: other constraints do not match with any edge\n";
      return NULL;
   }
   else
      error("TDecCDT::FindConstraintInVE, inconsistency: constr not found\n");

   return NULL; /*PAOLA*/
}



boolean TDecCDT::EqualsEV( PTEdge e1, PTEdge e2 ) {
	return (  ( e1->EV[0] == e2->EV[0] && e1->EV[1] == e2->EV[1] ) || 
		       ( e1->EV[0] == e2->EV[1] && e1->EV[1] == e2->EV[0] )     );
}


// --------------------------------------------------------------------------------
//  
//  void TDecCDT::WriteConstraints( const char * )
//
//  Save to file the constraint edges of the final CDT.
//  
//  Called by TDecCDT::WriteData1 and, unlike TDecCDT::WriteData2,
//  it uses array Constraints.
//

#include <stdio.h>
#include <stdlib.h>  // for perror
#include <string.h>  // fpr strerror
#ifdef CC_GCC
   #include <errno.h>
#endif

void TDecCDT::WriteConstraints( const char *outfname  ) 
{

	#define ios_base ios  // for compatibility VC5, GCC (non standard C++)

	#ifdef DEBUG3 
		DEBUG3 << "\nTDecCDT::WriteConstraints(" << outfname << ")\n";
	#endif

	ofstream outFile( outfname, ios_base::app );

	if( ! outFile )
	{
		fprintf( stderr, "TDecCDT::WriteConstraints(), <0> cannot open output file %s\n", outfname );
//		fprintf( stderr, "strerror says open failed: %s\n", strerror( errno ) );
	}
	check( !outFile, "" );

	outFile << nConstrInTRI << endl;;

	int e;
	for( e = 0; e < nConstrInFile; e++ ) {
		if(Constraints[e] != NULL)
		  outFile << Constraints[e]->EV[0]->VID << '\t' << Constraints[e]->EV[1]->VID << endl;

		if( e % 1000 == 0) cerr << "output " << e << " constraints\r";
	}

    cerr << "output " << nConstrInTRI << " constraints of " << nConstrInFile << "\n";

};



// --------------------------------------------------------------------------------
//  
//  void TDecCDT::WriteData2( const char * )
//
//  Similar to TDecimDelaunay::WriteData but use array Points and manage
//  constraints. Save final CDT to file, including constraint edges.
//  
//  Version of TDecCDT::WriteData that does NOT use array Constraints but
//  uses CONSTRAINED marks on edges of the triangulation.
//

#ifndef array_Constraints

void TDecCDT::WriteData( const char *outfname )
{
    int iv, it, ie, i;

     #ifdef DEBUG
      DEBUG << "\nTTriangulation::WriteData(" << outfname << ")" << endl;
     #endif // DEBUG

    //
    // create array to contian output vertices, triangles and constraints
    //
    
    // PAOLA: number of really existing vertices
    int nVrt = 0;

    PTVertex *VtxArray;   
    PTTriangle *TrgArray;
    PTEdge * EdgArray;

cerr << "TDecCDT::WriteData: nPts=" << nPts <<
" nTrg=" << nTrg << " nConstrInTRI=" << nConstrInTRI << endl;

    VtxArray = new PTVertex[ nPts ];
    TrgArray = new PTTriangle[ nTrg ];
    EdgArray = new PTEdge[ nConstrInTRI ];

    check( ( VtxArray == NULL || TrgArray == NULL || EdgArray == NULL),
           "TTriangulation::WriteData(), insufficient memory" );
	  
    for( iv=0; iv<nPts; iv++ ) VtxArray[iv] = NULL;
    for( it=0; it<nTrg; it++ ) TrgArray[it] = NULL;
    for( ie=0; ie<nConstrInTRI; ie++ ) EdgArray[ie] = NULL;
    
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
    ie = 0;

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

       // ...edges of CurTrg
       int i;
       for( i=0; i<3; i++ )
       {
          if( ! CurTrg->TE[i]->Marked( VISITED ) && CurTrg->TE[i]->Marked( CONSTRAINED )) // ogni edge e' comune a due triangoli quindi
          { // it may have been already visited
            #ifdef ROBUST
               check( ie >= nConstrInTRI, "TDecCDT::WriteData2(), constraints inconsistency detected" );
            #endif

            EdgArray[ie] = CurTrg->TE[i];
            EdgArray[ie]->Mark( VISITED );
            ie++;
          }
       }
       
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
     
cerr << "True " << it << " triangles" << endl;
cerr << "True " << iv << " vertices" << endl;
cerr << "True " << ie << " constraints" << endl;

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

       CT->UnMark( VISITED ); // smarchiamo i triangoli
    }

    cerr << "output " << nTrg << " triangles" << endl;

    //
    // ...finally edges
    //
    check( ie != nConstrInTRI, "TDecCDT::WiteData2(), <3> inconsistency detected");

    outFile << nConstrInTRI << endl;

    for( ie=0; ie<nConstrInTRI; ie++ )
    {
       int vIndx[2];

       vIndx[0] = EdgArray[ie]->EV[0]->VID;
       vIndx[1] = EdgArray[ie]->EV[1]->VID;

       outFile << vIndx[0] << '\t' << vIndx[1] << endl;
       
       if ( ie % 1000 == 0 ) cerr << "output " << ie << " triangles\r";
       EdgArray[ie]->UnMark( VISITED ); // smarchiamo i lati
    }
    
    cerr << "output " << nConstrInTRI << " constraints of " << nConstrInFile << "\n";       
    
    outFile.close();
    
}
#endif


//
// TDecCDT::Tra_Add_Constraint(PTVertex v0, PTVertex v1)
//
// Insert a constraint in the triangulation, assuming that it is admissible.
// A constraint is admissible if it does not contain any vertex inside it,
// and it does not intersect any other constraints (except in a common
// endpoint).
// Modify triangulation in such a way that it contains an edge of endpoints
// v0 and v1 and mark such edge as CONSTRAINED.
//

void TDecCDT::Tra_Add_Constraint( PTVertex v0, PTVertex v1 )
{
   #ifdef DEBUG8
      DEBUG8 << "TDecCDT::Tra_Add_Constraint( V" << v0->VID << " , V" << v1->VID << " )\n";
   #endif 

   // new constraint to be inserted
   #ifdef _GC_ON
	  PTEdge newConstr = GC::NewEdge( v0, v1 );
	#else
	  PTEdge newConstr = new TEdge( v0, v1 );
	#endif // _GC_ON

	check( (newConstr == NULL ),"TDecCDT::Tra_Add_Constraint , insufficient memory" );


   TDoubleList<PTEdge> Right_Border, Left_Border; 
   // right and left border w.r.t. the edge directed from v0 to v1

   //
   // if edge already exists and is not marked as CONSTRAINED,
   // then mark it and exit
   // (function CalcEdgInflRegn marks it when it finds it)
   //

   PTEdge origEdg;

   boolean constructed_borders = CalcEdgInflRegn( newConstr, Right_Border, Left_Border, origEdg );

   if( ! constructed_borders )
   {
      newConstr->Mark(CONSTRAINED);
      v0->nIncConstr++;
      v1->nIncConstr++;
      nConstrInTRI++;
      #ifdef array_Constraints
         for( int c = 0; c < nConstrInFile; c++ )
		      if( Constraints[c] == NULL )
	         { Constraints[c] = newConstr; break; }
		#endif
		
      return;
   }

   //
   // else: edge is inside and the two lists Right_Border and Left_Border
   // contain the closed sequence of edges bounding the two regions of
   // influence of  newConstr to its right and left side, respectively.
   //

   #ifdef ROBUST
      CheckEdgeInflBorder( newConstr, Right_Border );
      CheckEdgeInflBorder( newConstr, Left_Border );
   #endif
   
   // RetriangulateInfluenceRegion( TRUE, Right_Border );

   // DO NOT understand why after RetriangulateInfluenceRegion we expect
   // that newConstr is marked as INFL_BORDER
   cerr << "newConstr->Marked( INFL_BORDER ) = " << newConstr->Marked( INFL_BORDER ) << endl;
   Pause( "It should be always true" );
   newConstr->Mark( INFL_BORDER );
   // it should be already markes, tested with Pause
   
   // RetriangulateInfluenceRegion( TRUE, Left_Border );

   newConstr->Mark(CONSTRAINED);
	v0->nIncConstr++;
   v1->nIncConstr++;
   nConstrInTRI++;  
   #ifdef array_Constraints
		for( int c = 0; c < nConstrInFile; c++ )
		   if( Constraints[c] == NULL )
		      { Constraints[c] = newConstr; break; }
	#endif

   #ifdef ROBUST
		PTEdge ConstrEdge = new TEdge( v0, v1 );
                // do not use new_e because maybe already deleted
      
      PTEdge tmpE1 = FindEdgeInVE( v0, ConstrEdge );  // crash if not found
      PTEdge tmpE2 = FindEdgeInVE( v1, ConstrEdge );  // crash if not found
      check( tmpE1 != tmpE2, "TDecCDT::Tra_Add_Constraint(), <tmpE1 != tmpE2>");
      check( ! tmpE1->Marked( CONSTRAINED ), "TDecCDT::Tra_Add_Constraint(), <! tmpE1->Marked( CONSTRAINED )>");
      TDoubleList<PTEdge> ConstrList;
      check( FindConstraintsInVE( v0, ConstrList ) != (v0->nIncConstr), "TDecCDT::Tra_Add_Constraint(), V1->nIncConstr inconsistent");
      ConstrList.ClearList();
      check( FindConstraintsInVE( v1, ConstrList ) != (v1->nIncConstr), "TDecCDT::Tra_Add_Constraint(), V2->nIncConstr inconsistent");
      ConstrList.ClearList();
   #endif
}


// -----------------------------------------------------------------------------
//  
//   int TDecCDT::CalcEdgInflRegn( PTEdge newConstr, TDoubleLisr<PTEdge> Right_Border,
//                                                   TDoubleLisr<PTEdge> Left_Border )
//
//   Compute region of influence of new constraint to be inserted into the
//   current triangulation, i.e., the lists of boundary edges of the two
//   regions to the right and to the left of newConstr. Such lists will
//   contain the edges in counterclockwise order.
//
//   INPUT: newConstr
//
//   OUTPUT: fill the doubly linked lists Right_Border and Left_Border:
//   such lists contain the right region and the left region of newConstr
//   (considering it as directed from newConstr->EV[0] to newConstr->EV[1]).
//   The edges of the two lists are marked as INFL_BORDER, and the edges
//   and triangles inside the two regions are marked as TO_DELETE.
//   The TTriangle pointer FirstTrgToDel is set to point one of the 
//   triangles of the region of influence Right_Border.
//
//   RETURN VALUE:
//         FALSE  if edge already exists (and is not marked as CONSTRAINED
//           - to be further investigated).
//           In this case the two border lists are not considered because
//           it is sufficient to mark the edge as CONSTRAINED.
//           It includes the case in which newConstr is on the convex hull
//           of the triangulation.
//         TRUE   if edge newConstr is internal and the two border lists  
//         are filled correctly.
//

boolean TDecCDT::CalcEdgInflRegn( PTEdge newConstr, TDoubleList<PTEdge> &Right_Border,
                                  TDoubleList<PTEdge> &Left_Border, PTEdge & origEdg   )
{
   #ifdef DEBUG7
      DEBUG7 << "\n\nTDecCDT::CalcEdgInflRegn( " << *newConstr << " )\n";
   #endif 

   #ifdef ROBUST
      check( ! Right_Border.IsEmpty() , "TDecCDT::CalcEdgInflRegn(..) : called on Right_Border not empty" );
      check( ! Left_Border.IsEmpty() , "TDecCDT::CalcEdgInflRegn(..) : called on Left_Border not empty" );
   #endif

   int intersection;
   // result of intersection of newConstr with the edges of examined
   // triangles

   PTVertex v0 = newConstr->EV[0],  // vertices of new constraint 
            v1 = newConstr->EV[1];

   //
   //  examine v0
   //

   // first and last triangle of VT( v0 )
   PTTriangle firstTrg = Next_CCW_t_of_e_by_v( FirstEdgeInVE( v0 ) , v0 );
   PTTriangle lastTrg = LastTrgInVE( v0 );

   PTTriangle currTrg = NULL;  // current triangle of VT( v0 )
   PTEdge oppositeEdg = NULL;  // edge of currTrg opposite to v0


   PTTriangle adj_v0_Trg = NULL, adj_v1_Trg = NULL;
   // triangles adjacent to v0/v1 and properly intersecting newConstr

   // scan the VT relation of v0 in counterclockwise order searching
   // for a triangle Trg such that the edge of Trg opposite to v0
   // intersects newConstr
   for( currTrg = firstTrg; /*currTrg != lastTrg*/; currTrg = Next_CCW_t_of_t_by_v( currTrg, v0 ) )
   {
      oppositeEdg = Opposite_e_of_t_by_v( currTrg, v0 );

      intersection = Tra_intersect( newConstr, oppositeEdg );

      if( intersection == PROPER_INTER )
      {
         #ifdef ROBUST
            check( oppositeEdg->Marked( CONSTRAINED ), "TDecCDT::CalcEdgInflRegn, <1> inconsistency detected" );
         #endif

         adj_v0_Trg = currTrg;

         #ifdef DEBUG7
            DEBUG7 << "   adj_v0_Trg = " << *adj_v0_Trg << endl;
         #endif

         break;
      }
      else
      if( intersection == ONLY_1CV )
      // edge newConstr that we are going to insert is the same as 
      // an edge already present in the triangulation
      {

         // currTrg may be firstTrg: in this case the edge candidate to be
         // the same as newConstr may be the one preeceding it in 
         // CCW w.r.t. oppositeEdg.
         PTEdge nextEdg = Next_CCW_e_of_t_by_e( currTrg, oppositeEdg ),
                prevEdg = Prev_CCW_e_of_t_by_e( currTrg, oppositeEdg );

         // set origEdg to point to the edge of the current triangulation
         // which is the same as newConstr, and return FALSE
         if( EqualsEV( newConstr, nextEdg ) ) // currTrg is not firstTrg
            origEdg = nextEdg;
         else                            // currTrg is firstTrg
         {
            origEdg = prevEdg;
            #ifdef ROBUST
               if( ! EqualsEV( newConstr, prevEdg ) )
               {
                  cerr << "TDecCDT::CalcEdgInflRegn( " << *newConstr << ", ... ) : the new edge "
                       << *newConstr << "\nintersects edge " << *oppositeEdg << " in only 1 common vertex;\nbut neither "
                       << *nextEdg << " nor " << *prevEdg << " matches  the new edge\n";
                  Pause("");
               }
            #endif
         }
         return FALSE;
      }

      check( currTrg == lastTrg, "TDecCDT::CalcEdgInflRegn, triangle intersecting V0-V1 not found in VT( V0 )" );
   }


   //
   //  examine v1 ( as v0 )
   //
   firstTrg = Next_CCW_t_of_e_by_v( FirstEdgeInVE( v1 ) , v1 );
   lastTrg = LastTrgInVE( v1 );;

   for( currTrg = firstTrg; /*currTrg != lastTrg*/; currTrg = Next_CCW_t_of_t_by_v( currTrg, v1 ) )
   {
      oppositeEdg = Opposite_e_of_t_by_v( currTrg, v1 );

      intersection = Tra_intersect( newConstr, oppositeEdg );

      if( intersection == PROPER_INTER )
      {
         #ifdef ROBUST
            check( oppositeEdg->Marked( CONSTRAINED ), "TDecCDT::CalcEdgInflRegn, <2> inconsistency detected" );
         #endif

         adj_v1_Trg = currTrg;

         #ifdef DEBUG7
            DEBUG7 << "   adj_v1_Trg = " << *adj_v1_Trg << endl;
         #endif

         break;
      }
      else
      {
         #ifdef ROBUST
            // if newE did not exixted already, we would have realized it
            // while examining v0
            check( intersection == ONLY_1CV, "TDecCDT::CalEdgInflRegn(..) : <2> inconsistency detected\n" );
            if( intersection != NO_INTER ) Pause( "TDecCDT::CalEdgInflRegn(..) : <3> inconsistency detected\n" );
         #else
            ;
         #endif
      }
      check( currTrg == lastTrg, "TDecCDT::CalcEdgInflRegn, triangle intersecting V0-V1 not found in VT( V1 )" );
   }

   //
   // Examine triangles intersected by newCostr between v0 and v1:
   // start from adj_v0_Trg and navigate through adjacency links until
   // we arrive at adj_v1_Trg; 
   // - encountered triangles must be left marked as TO_DELETE,
   // - edges intersecting newConstr must be marked as TO_DELETE,
   // - edges NOT intersecting newConstr must be marked as INFL_BORDER 
   //   and inserted in the corresponding boundary list.
   //
   // Insert boundary edges on the right of newConstr (directed from v0 to
   // v1) to the end of Right_Border, and insert biundary edges on the left
   // to the beginning of Left_Border:
   // this guarantees that the resulting lists are sorted counterclockwise.
   //

      // Examine triangles adj_v0_Trg and adj_v1_Trg:
      // init lists Right_Border and Left_Border with newConstr and
      // with boundary edges incident in adj_v0_Trg and adj_v1_Trg.

   Right_Border.AddTail( Left_e_of_t_by_v( adj_v1_Trg, v1 ) );
   Right_Border.AddTail( newConstr );
   Right_Border.AddTail( Right_e_of_t_by_v( adj_v0_Trg, v0 ) );

   Left_Border.AddHead( Right_e_of_t_by_v( adj_v1_Trg, v1 ) );
   Left_Border.AddHead( newConstr );
   Left_Border.AddHead( Left_e_of_t_by_v( adj_v0_Trg, v0 ) );

      // mark initial boundary edges and initial internal triangles
   newConstr->Mark( INFL_BORDER );
   Right_Border.GetHead()->Mark( INFL_BORDER );
   Right_Border.GetLast()->Mark( INFL_BORDER );
   Left_Border.GetHead()->Mark( INFL_BORDER );
   Left_Border.GetLast()->Mark( INFL_BORDER );
   adj_v0_Trg->Mark( TO_DELETE );
   adj_v1_Trg->Mark( TO_DELETE );

      // now examine triangles between adj_v0_Trg and adj_v1_Trg

   PTTriangle prevTrg;
   // triangle adjacent to current triangle currTrg along edge 
   // intersectingEdg properly intersecting newConstr

   prevTrg = adj_v0_Trg;
   
   oppositeEdg = Opposite_e_of_t_by_v( adj_v0_Trg , v0 );

   PTEdge intersectingEdg = oppositeEdg;
   // current edge intersecting newConstr

   currTrg = Opposite_t_of_t_by_e( adj_v0_Trg , intersectingEdg ); 
   // currently examined triangle

   // intersectingEdg and currTrg are always internal to the region:
   // they must be marked as TO_DELETE
   intersectingEdg->Mark( TO_DELETE );
   currTrg->Mark( TO_DELETE );

   // positions of the three edges in currTrg->TE
   int intersectingEdg_pos, nextintersectingEdg_pos, notintersectingEdg_pos;  

 
   while( currTrg != adj_v1_Trg )
   {
      // search among the edges of currTrg the positions of the three edges
      //  - intersectingEdg  (certainly internal, it must be marked TO_DELETE),
      //  - the other edge, different from intersectingEdg, properly
      //    intersecting newConstr
      //  - the edge not intersecting newConstr (it is on the boundary
      //    of the region, it must be marked as INFL_BORDER).
      for( int i = 0; i < 3; i++ )
      {
         if( currTrg->TE[i] == intersectingEdg )
         {
            intersectingEdg_pos = i;
         }

         // currTrg->TE[i] is the other edge of currTrg, different from
         // intersectingEdg, that intersects newConstr, and thus it must 
         // not be inserted in one of the two lists
         else if( (intersection = Tra_intersect( newConstr , currTrg->TE[i])) )
         {
            nextintersectingEdg_pos = i;

            #ifdef ROBUST
               if( intersection != PROPER_INTER  ||  currTrg->TE[i]->Marked( CONSTRAINED ) )
               {   cerr << "ERROR: TDecCDT::CalcEdgInflRegn( " << *newConstr << ", ...), <c> inconsistency\n"; exit(-1); }
            #endif
         }
         else  
         // currTrg->TE[i] is the edge of currTrg that does NOT intersect
         // newConstr, and thus it MUST be inserted in one of the two lists
         {
            notintersectingEdg_pos = i;
         }
      } // ...for

      #ifdef ROBUST
         if( intersectingEdg_pos == notintersectingEdg_pos  ||  intersectingEdg_pos == nextintersectingEdg_pos )
         {   cerr << "ERROR: TDecCDT::CalcEdgInflRegn( " << *newConstr << ", ...), <c> inconsistency\n"; exit(-1); }
      #endif


      // Insert the edge that does not intersect newConstr in one of the 
      // two lists:
      // - if it follows (in CCW order) intersectingEdg in the TE relation
      //   of currTrg, then it is on the right of  newConstr,
      // - otherwise it is on the left of  newConstr.
      if( notintersectingEdg_pos == (intersectingEdg_pos+1)%3 )
         // edge right of newConstr
         Right_Border.AddTail( currTrg->TE[notintersectingEdg_pos ] );
      else                                                      
         // edge left of newConstr
         Left_Border.AddHead( currTrg->TE[notintersectingEdg_pos ] );

      // mark edge that does not intersect newConstr
      currTrg->TE[notintersectingEdg_pos]->Mark( INFL_BORDER );

      // update intersectingEdg and currTrg, and mark them

      intersectingEdg = currTrg->TE[nextintersectingEdg_pos];
      currTrg = Opposite_t_of_t_by_e( currTrg, intersectingEdg );

      intersectingEdg->Mark( TO_DELETE );
      currTrg->Mark( TO_DELETE );

   } // ...while


   //
   // set FirstTrgToDel to point one of the old triangles of the region of
   // influence of newConstr: the one in VT(v0) intersecting newConstr:
   // adj_v0_Trg  
   //
   FirstTrgToDel = adj_v0_Trg;


   // adjust the two listst in such a way that their first element is
   // newConstr
   //
   Right_Border.AddTail( Right_Border.RemoveHead() );
   Left_Border.AddHead( Left_Border.RemoveLast() );
   Left_Border.AddHead( Left_Border.RemoveLast() );

   return TRUE;
}


//
// Classify the removability of a vertex V based on the number nIncConstr 
// of constraints incident in it, and based on the value of variables
// AllowFeaturesDel and AllowChainBrk. 
//
// INPUT: a vertex V
//
// OUTPUT: TRUE if  nIncConstr == 0
//           OR ( nIncConstr == 1 AND AllowFeaturesDel == true )
//           OR ( nIncConstr == 2  AND  AllowChainBrk==true   AND 
//              the new constraint c that should replace the two constraints
//              is admissible )
//           OR ( nIncConstr == 2  AND  AllowChainBrk==false  AND 
//              c is admissible and if not already a constraint in the
//              triangulation )
//         FALSE  otherwise. 
//
// ATTENTION: if edges represent iso-lines, they are formed by closed 
// chains, otherwise they are features and thus it makes sense to 
// simplify them but not to eliminate them completely.
// If V has two incidemnt constraints C1 and C2 incidenti and newConstr
// is the same as an edge E of the triangulation, if E is not a constraint
// then return true.
// This choice is motivated by the fact that in case newConstr is the same
// as an existing constraint C, then the sequence of constraints
// < C C1 C2 > is a polygon and we do not want it to disappear.
//
// ALGORITHM:
//    if( nIncConstr  == 0  ||  ( nIncConstr  == 1 && AllowFeaturesDel == true ) )
//         return TRUE
//    else if( V->nIncConstr == 2 )
//      {  
//         V1 = Other_v_in_e( IncConstr.GetHead(), V );  V2 = Other_v_in_e( IncConstr.GetLast(), V );
//         if(    ( AllowChainBrk==true  &&  ( NEW_EDGE || EXIST || EXIST_CONSTRAINED ) == IsEdgeAdmissible( V1, V2 ) )
//             || ( AllowChainBrk==false  && ( NEW_EDGE || EXIST ) == IsEdgeAdmissible( V1, V2 ) )    )
//             return TRUE;
//      }
//    else return FALSE;
//

boolean TDecCDT::OkConstrDegree( PTVertex V )
{
   #ifdef DEBUG10
      DEBUG10 << "\n\nTDecCDT::OkConstrDegree( " << *V << " )\n";
   #endif 

   #ifdef ROBUST
      TDoubleList<PTEdge> IncConstr;
      check( V->nIncConstr != FindConstraintsInVE( V, IncConstr ) ,
             "TDecCDT::OkConstrDegree(V) : inconsistency: value V->nIncConstr does not match the CDT" );
      check( V->nIncConstr < 0 , "TDecCDT::OkConstrDegree(V) : constrDegree(V) < 0" );
   #endif

//   if( TDecimDelaunay::InitialPhase )  // in initial phase RemoveVertex is never called
//      nIncConstr = V->nIncConstr;

   if( V->nIncConstr == 0   ||   ( V->nIncConstr == 1  &&  AllowFeaturesDel == TRUE )   )
      return TRUE;
   else
      if( V->nIncConstr == 2 )
      {
         PTVertex V1 = Other_v_in_e( IncConstr.GetHead(), V );
         PTVertex V2 = Other_v_in_e( IncConstr.GetLast(), V );
         int isadmissible = IsEdgeAdmissible( V1, V2 );

         if(  ( AllowChainBrk == TRUE  &&  isadmissible >= EXIST_CONSTRAINED ) ||
              ( AllowChainBrk == FALSE &&  isadmissible >= EXIST             )     )
            return TRUE;
         else return FALSE;
      }
      else return FALSE;
}


// -----------------------------------------------------------------------------
//  
//   boolean TDecCDT::ReCheckVertex( PTVertex V )
//
//   Check if vertex V is removable after the re-triangulation of the
//   region of influence it belongs to the boumndary of.
//   Re-implemented to check also ConstraintDegee(V).
//

boolean TDecCDT::ReCheckVertex( PTVertex V )
{
   return( IsVtxElim(V) && OkDegree( CalcDegree(V) ) && OkConstrDegree( V ) );
}


int TDecCDT::IsEdgeAdmissible( PTVertex v0, PTVertex v1 )
{
   #ifdef DEBUG11
      DEBUG11 << "\n\nTDecCDT::IsEdgeAdmissible( " << *v0 << ", " << *v1 <<  " )\n";
   #endif

	PTEdge newE = new TEdge( v0, v1 );

   
   int intersection;
   // result of intersection of newE with edges of examined triangles

   //
   //  examine v0
   //

   // first and last triangle of VT( v0 )
   PTTriangle firstTrg = Next_CCW_t_of_e_by_v( FirstEdgeInVE( v0 ) , v0 );
   PTTriangle lastTrg = LastTrgInVE( v0 );

   PTTriangle currTrg;  // current triangle of VT( v0 )
   PTEdge oppositeEdg;  // edge of currTrg opposite to v0


   PTTriangle adj_v0_Trg, adj_v1_Trg;
   // adjacent triangles to v0/v1 and properly intersecting newE

   // scandisco the VT relation of v0 in counterclockwise order searching
   // for a triangle Trg such that the edge of Trg opposite to v0
   // intersects E
   for( currTrg = firstTrg; /*currTrg != lastTrg*/; currTrg = Next_CCW_t_of_t_by_v( currTrg, v0 ) )
   {
      oppositeEdg = Opposite_e_of_t_by_v( currTrg, v0 );

      intersection = Tra_intersect( newE, oppositeEdg );

      if( intersection == PROPER_INTER )
      {
         if( oppositeEdg->Marked( CONSTRAINED ) )
            return NOT_ADMISSIBLE;

         adj_v0_Trg = currTrg;
         break;
      }
      else
      if( intersection == ONLY_1CV )
      // edge newE( V0, V1 ) that we want to insert is the same as edge
      // E( Va, Vb ) already present in the current triangulation
      { 

         PTEdge matchingEdg;  // edge that is the same as newE

         // currTrg may be firstTrg: in this case the edge candidate
         // to be the same as newE can be the one preceding it in
         // CCW order w.r.t. oppositeEdg
         PTEdge nextEdg = Next_CCW_e_of_t_by_e( currTrg, oppositeEdg ),
                prevEdg = Prev_CCW_e_of_t_by_e( currTrg, oppositeEdg );

         if( EqualsEV( newE, nextEdg ) ) // currTrg is not firstTrg
            matchingEdg = nextEdg;
         else                            // currTrg is firstTrg
         {
            matchingEdg = prevEdg;
            #ifdef ROBUST
               if( ! EqualsEV( newE, prevEdg ) )
               {
                  cerr << "TDecCDT::IsEdgeAdmissible( " << *v0 << ", " << *v1 << " ) : the new edge "
                       << *newE << "\nintersects edge " << *oppositeEdg << " in only 1 common vertex;\nbut neither "
                       << *nextEdg << " nor " << *prevEdg << " matches  the new edge\n";
                  Pause("");
               }
            #endif
         }

         // if the edge of the current triangulation which is the same as
         // newE then return EXIST_CONSTRAINED, else return EXIST
         if( matchingEdg->Marked( CONSTRAINED ) )
            return EXIST_CONSTRAINED;
         else
            return EXIST;
      }
      else
      if( intersection == NOT_PROPER )
      // l'edge newE( V0, V1 ) that we want to insert contains a vertex
      // of the current triangulation inside it
      {                                 
         return NOT_ADMISSIBLE;          
      }
      #ifdef ROBUST
         if( intersection != NO_INTER ) 
            Pause( "TDecCDT::IsEdgeAdmissible(..) : <2> inconsistency detected\n" );
 
      #endif
      
      if( currTrg == lastTrg )
      { 
         cerr << "ERRORE: TDecCDT::IsEdgeAdmissible: triangle intersecting " << *newE
              << " not found in VT( " << *v0 << " )\n";
         Pause("");
      }
   }

   //
   //  examine v1
   //
   
   firstTrg = Next_CCW_t_of_e_by_v( FirstEdgeInVE( v1 ) , v1 );
   lastTrg = LastTrgInVE( v1 );;

   for( currTrg = firstTrg; /*currTrg != lastTrg*/; currTrg = Next_CCW_t_of_t_by_v( currTrg, v1 ) )
   {
      oppositeEdg = Opposite_e_of_t_by_v( currTrg, v1 );

      intersection = Tra_intersect( newE, oppositeEdg );

      if( intersection == PROPER_INTER )
      {
         if( oppositeEdg->Marked( CONSTRAINED ) )
            return NOT_ADMISSIBLE;

         adj_v1_Trg = currTrg;
         break;
      }
      else
      if( intersection == NOT_PROPER )
      // edge newE( V0, V1 ) that we want to insert contains a vertex of
      // the current triangulation inside it
      {
         return NOT_ADMISSIBLE;          
      }
      #ifdef ROBUST
         if( intersection != NO_INTER ) 
            Pause( "TDecCDT::IsEdgeAdmissible(..) : <3> inconsistency detected\n" );
 
      #endif
      
      if( currTrg == lastTrg )
      { 
         cerr << "ERRORE: TDecCDT::IsEdgeAdmissible: triangle intersecting " << *newE
              << " not found in VT( " << *v1 << " )\n";
         Pause("");
      }
   }

   //
   // examine triangles intersected by newCostr  between v0 and v1:
   // start from adj_v1_Trg and navigate through adjacency links until
   // we arrive at adj_v0_Trg
   //

   PTTriangle prevTrg;
   // triangle adjacent to the current triangle currTrg along edge
   // intersectingEdg properly intersecting newE

   prevTrg = adj_v1_Trg;

   PTEdge intersectingEdg = oppositeEdg;

   currTrg = Opposite_t_of_t_by_e( currTrg, intersectingEdg );

 
   while( currTrg != adj_v0_Trg )
   {
      // search among the edges of currTrg the one properly intersecting
      // newE
      for( int i = 0; i < 3; i++ )
      {
         if( currTrg->TE[i] != intersectingEdg )
         {
            if( (intersection = Tra_intersect( newE , currTrg->TE[i] )) )
            {
               if( intersection == PROPER_INTER )
               {
                  if( ! currTrg->TE[i]->Marked( CONSTRAINED ) )
                  {
                     intersectingEdg = currTrg->TE[i];
                     currTrg = Opposite_t_of_t_by_e( currTrg, intersectingEdg );
                     break;
                  }
                  else return NOT_ADMISSIBLE;
               }
               else return NOT_ADMISSIBLE;
            }
         }
      }
   }
   return NOT_EXISTENT;
}


// -----------------------------------------------------------------------------
//  
//   void TDecCDT::UpdateStep()
//
//   Perform one update step of the current triangulation by eliminating 
//   one vertex.
//   The steps are: find the next vertex to be removed (the minimum
//   verted in tree ElimVtxTree) by calling NextVertex();
//   deletion of such vertex by calling RemoveVertex().
//
//   Re-defined temporarily to avoid using (and writing) function
//   RechekConstrAdjVtxs.
//

void TDecCDT::UpdateStep()
{
   #ifdef DEBUG
      DEBUG << "\nTDecCDT::UpdateStep()" << endl;
     DEBUG << "ntrg = " << nTrg << endl;
   #endif
   
   //
   // TDestroyDelaunay::NextVertex(), put in VertexToRemove the pointer
   // to the next vertex to be removed
   // 
   
   NextVertex();
   
   //
   // now vertex VertexToRemove is removed from the triangulation
   //
   if( OkConstrDegree( VertexToRemove ) )
      RemoveVertex();
   else
   {  cerr << "TDecCDT::UpdateStep() : NextVertex() has selected " << *VertexToRemove
           << " which is not OkConstrDegree\n"; Pause(""); }
}



void TDecCDT::Edge2Constraint( PTEdge E )
{
   nConstrInTRI++;
   E->EV[0]->nIncConstr++;
   E->EV[1]->nIncConstr++;
   E->Mark(CONSTRAINED);

   #ifdef array_Constraints
      for( int c = 0; c < nConstrInFile; c++ )
         if( Constraints[c] == NULL )
         {
            Constraints[c] = E;
            break;
         }
   #endif
}

void TDecCDT::AdjustBordersVertices_VE( TDoubleList<PTEdge> Border )
{
   //
   // Adjiust VE relations for all endpoint vertices of edges in Border.
   // Set such relations, which may point to internal edges of the region
   // to be removed, in such a way that every vertex always points
   // to the two edges in Border adjacent to the vertex.
   //
   
   TDoubleListIterator<PTEdge> IRegnIter( &Border );
   IRegnIter.Restart();
   
   PTEdge ECurr;
   PTEdge EPrev = Border.GetLast();
   
   PTVertex VCurr = NULL;
   
   while ( ! IRegnIter.EndOfList() )
   {
      ECurr = IRegnIter.Current()->object;
      
      // VCurr is the common vertex of EPrev and ECurr, i.e., the
      // second vertex of EPrev in CCW order w.r.t. VertexToRemove
      
      if ( EPrev->EV[0] == ECurr->EV[0] || EPrev->EV[0] == ECurr->EV[1] )
        VCurr = EPrev->EV[0];
      else
        VCurr = EPrev->EV[1];
      
      #ifdef ROBUST
         check( (ECurr->EV[0]!=VCurr && ECurr->EV[1]!=VCurr),
	      "TDecCDT::AdjustBordersVertices_VE(..), <1> inconsistency detected" ); 
      #endif
      
      
      // set VE relations of VCurr to point to ECurr/EPrev, if they
      // point to internal edges of the region to be deleted

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
      
      // else ... !E0Del && !E1Del... leave VE unchanged

     if (CheckIsolatedPoint(VCurr)) 
       cout << "AdjustBordersVertices_VE" << endl;

      // go to next edge
      
      EPrev = ECurr;
      IRegnIter.GoNext();
      
      
   } // end ... while ( ! EndOfList() )
  
}


void TDecCDT::Delete_Constr_InfluenceRegion( PTEdge newConstr )
{

   #ifdef DEBUG
     DEBUG << "TDecCDT::Delete_Constr_InfluenceRegion()" << endl;
   #endif
    
   //
   // delete old triangles
   //
   
     if ( ! TDecimDelaunay::InitialPhase ) TDecimDelaunay::MT_KillInterference();
    
   //
   // Delete newConstr from Left_Border and Right_Border, and append the
   // two lists in InflRegnBorder, to form a unique chain of edges
   // bounding the region of influence of newConstr
   //

   check( newConstr != Right_Border.RemoveHead(), "TDecCDT::Delete_Constr_InfluenceRegion(..) : <1> inconsistency detected\n" );
   check( newConstr != Left_Border.RemoveHead(), "TDecCDT::Delete_Constr_InfluenceRegion(..) : <2> inconsistency detected\n" );
   #ifdef ROBUST
      check( ! InflRegnBorder.IsEmpty(), "TDecCDT::Delete_Constr_InfluenceRegion(..) : <3> inconsistency detected\n" );
      int nEdges = Right_Border.Lenght() + Left_Border.Lenght();
   #endif
   InflRegnBorder.Append( Right_Border );
   InflRegnBorder.Append( Left_Border );
   #ifdef ROBUST
      check( nEdges != InflRegnBorder.Lenght(), "TDecCDT::Delete_Constr_InfluenceRegion(..) : <4> inconsistency detected\n" );
      check( 0 != Right_Border.Lenght(), "TDecCDT::Delete_Constr_InfluenceRegion(..) : <5> inconsistency detected\n" );
      check( 0 != Left_Border.Lenght(), "TDecCDT::Delete_Constr_InfluenceRegion(..) : <6> inconsistency detected\n" );
   #endif


   //
   // call old version of DeleteInfluenceRegion() to delete the triangles
   //
   
   TDestroyDelaunay::DeleteInfluenceRegion();


   // 
   // TDestroyDelaunay::DeleteInfluenceRegion() scans the list
   // InflRegnBorder and reconsiders the removability of each vertex 
   // V encountered (VertexToRemove may be among them):
   //   - delete it from ElimVtxTree,
   //   - do RecheckVertex to test removability
   //   - if vertex is removable then insert it in ElimVtxTree
   //  Obviously we do not want that VertexToRemove is inserted in
   //  ElimVtxTree: then, if it has been inserted, we remove it.
   //
   if( ElimVtxTree.IsIn( VertexToRemove ) )
      ElimVtxTree.Remove( VertexToRemove );

   //
   // Reposition "Detached" points from old triangles into new triangles,
   // in order to compute their error correctly in MT_AddComponent().
   //
//   TDecimDelaunay::RepositionPoint( VertexToRemove );
     // here we have not re-triangulated the region of influence of
     // VertexToRemove: it is not to be repositioned
                                                       
   while( ! TDecimDelaunay::DetachedPoints.IsEmpty() )
   {
       PTPoint p = TDecimDelaunay::DetachedPoints.RemoveHead();
       TDecimDelaunay::RepositionPoint( p );       
   } 

   //
   // now add new triangles (the ones marked as NEW_TRIANGLE) and unmark
   // them
   //
   
   if ( ! TDecimDelaunay::InitialPhase ) MT_AddComponent();


   //
   // remove mark INFL_BORDER from edges in InflRegnBorder
   //
      
   TDoubleListIterator<PTEdge> IRegnIter( &InflRegnBorder );
   IRegnIter.Restart();
   while ( ! IRegnIter.EndOfList() )
   {
       IRegnIter.Current()->object->UnMark( INFL_BORDER );
       IRegnIter.GoNext();
   }
      
   // empty InflRegnBorder, is it no longer needed
   
   InflRegnBorder.ClearList();


   //
   // remove mark INFL_BORDER from newConstr
   //
   newConstr->UnMark( INFL_BORDER );

}


void TDecCDT::TTra_Add_Constraint( PTVertex v0, PTVertex v1 )
{
   #ifdef DEBUG13
      cerr << "TDecCDT::TTra_Add_Constraint( V" << v0->VID << " , V" << v1->VID << " )\n";
      cerr << "Now constr num = " << nConstrInTRI << endl;
   #endif 

   // new constraint to be inserted
   #ifdef _GC_ON
	  PTEdge newConstr = GC::NewEdge( v0, v1 );
	#else
	  PTEdge newConstr = new TEdge( v0, v1 );
	#endif // _GC_ON

	check( (newConstr == NULL ),"TDecCDT::TTra_Add_Constraint , insufficient memory" );


   PTEdge origEdg;

   if( ! CalcEdgInflRegn( newConstr, Right_Border, Left_Border, origEdg ) )
   {  
      delete newConstr;  
      if( ! origEdg->Marked( CONSTRAINED ) )
         Edge2Constraint( origEdg );
   }
   else
   {
      Edge2Constraint( newConstr );
      AdjustBordersVertices_VE( Right_Border );
      AdjustBordersVertices_VE( Left_Border );
      RIRInitialTrg( Right_Border, SwapEdgeQueue );
      RIRInitialTrg( Left_Border, SwapEdgeQueue );
      RIRDelOptTrg();
      Delete_Constr_InfluenceRegion( newConstr );
   }
}


void TDecCDT::Constraint2Edge( PTEdge E )
{
   nConstrInTRI--;
   E->EV[0]->nIncConstr--;
   E->EV[1]->nIncConstr--;
   E->UnMark(CONSTRAINED);

   #ifdef array_Constraints
      for( int c = 0; c < nConstrInFile; c++ )
         if( Constraints[c] == E )
         {
            Constraints[c] = NULL;
            break;
         }
   #endif
}


//
// Let c be a constraint incident in VertexToRemove.
// Insert in SwapEdgeQueue and mark as SWAP the non-constraint edges
// which do not belong to the convex hull, belong to InflRegnBorder 
// and are incident in c.
//
void TDecCDT::ExtRIRIPrepare(TDoubleList<PTEdge> & SwapEdgeQueue, PTEdge c,
                             TDoubleList<TEdge> &OrigEdgList)
{
   PTEdge e;
   for( int i = 0; i <= 1; i++ )
   {
      if( c->ET[i] != NULL )
      {
         e = Opposite_e_of_t_by_v( c->ET[i], VertexToRemove);
         if( ! ( e->Marked(CONSTRAINED) || e->OnConvexHull() )  )
         {
            SwapEdgeQueue.AddHead(e);
            e->Mark(SWAP_EDGE_QUEUE);
//          OrigEdgList.AddHead(*e);
         }
      }
   }
}


//
// Insert in SwapEdgeQueue and mark as SWAP the non-constraint edges
// which do not belong to the convex hull, belong to InflRegnBorder 
// and are incident in an original constraint incident in VertexToRemove.
// Let c1 and c2 be such constraint(s).
// 
// Convert c1 and c2 to normal edges.
//
// Put into OrigEdgList a copy of the edges present in InflRegnBorder
// and mark them as COPIED.
//
// In case the boundary InflRegnBorder of the influence region is not
// a closed sequence of edges:
// create a new edge connecting the first and the last boundary vertex
// and insert it in InflRegnBorder; in the new triangulation,
// such edge will replace the two edges of the convex hull incident 
// in VertexToRemove (not present in InflRegnBorder).
//
// For each edge E belonging to InflRegnBorder, for each vertex V of E:
// adjust VE relations in such a way that, if they pointed to edges 
// internal to the region of influence, now they point to edges of its
// boundary.
//
void TDecCDT::ExtRIRIPrepare(TDoubleList<PTEdge> & SwapEdgeQueue, PTEdge c1, PTEdge c2,
                             boolean & closed, TDoubleList<TEdge> &OrigEdgList)
{
   if( c1 != NULL )
   {
       ExtRIRIPrepare (SwapEdgeQueue, c1, OrigEdgList);
       Constraint2Edge(c1);
   }
   if( c2 != NULL )
   {
       ExtRIRIPrepare (SwapEdgeQueue, c2, OrigEdgList);
       Constraint2Edge(c2);
   }

   // Put a copy of the edges of InflRegnBorder in OrigEdgList, and mark
   // them as COPIED.
   // Do not worry about parameter "closed": an edge of the convex hull
   // is never considered for the swap during the optimization.
   TDoubleListIterator<PTEdge> InflRegnIter( &InflRegnBorder );
   InflRegnIter.Restart();
   while ( ! InflRegnIter.EndOfList() )
   {
       OrigEdgList.AddHead( *(InflRegnIter.Current()->object) );
       InflRegnIter.Current()->object->Mark(COPIED);
       InflRegnIter.GoNext();
   }


   // Create a new edge joining the first and the last boundary vertex of
   // the region of influence and insert it in InflRegnBorder.
   // In the new triangulation such edge will replace  the two edges of 
   // the convex hull incident in VertexToRemove (not present in
   // InflRegnBorder).
   //
   // Step 1: If InflRegnBorder is open create a new edge joining its first
   // and the last boundary vertex; the Et relations of such edge are 
   // still NULL. Make sure not to use them before they are set.
   //
   
   PTEdge NewEdge; 
   // edge added to the triangulation and to InflRegnBorder
      
   //
   // find VFirst and VLast, first and the last vertex of
   // tyhe sequence of edges in InflRegnBorder
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
               "TDecCDT::ExtRIRPrepare(), <1> inconsistency detected" );
   #endif // ROBUST   
   
   if ( !closed )
   {
       // create new edge and put it into the list
       
       NewEdge = new TEdge( VFirst, VLast );
       check( (NewEdge == NULL), "TDestroyDelaunay::RetriangulateInfluenceRegion(), insufficient memory" );
       
       
       // temporarily adjast VE relations of VFirst and VLast (now they
       // point to edges of the hull to be deleted)
       
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
       
if (CheckIsolatedPoint(VFirst) || CheckIsolatedPoint(VLast))
   cout << "ExtRIRIPrepare (1)" << endl;
       
       #ifdef DEBUG
          DEBUG << "Create edge E" << NewEdge->EID 
	        << "( V" << VFirst->VID << ", V" << VLast->VID << " )" << endl;
       #endif
       
       InflRegnBorder.AddTail( NewEdge );
       NewEdge->Mark( INFL_BORDER );
       
       
       #ifdef ROBUST
         // now set closed = TRUE and VLast = VFirst (probably we do not
         // use such variables any more, but do this for consistency
	 closed = TRUE;
	 VLast  = VFirst;
       #endif // ROBUST       
       
   }

   //
   // Step 2: Adjiust VE relations for all endpoint vertices of edges in
   // InflRegnBorder. Set such relations in such a way that each vertex
   // points to the two edges in InflRegnBorder adjacent to it.
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
      // vertex of EPrev in CCW order w.r.t. VertexToRemove 
      
      if ( EPrev->EV[0] == ECurr->EV[0] || EPrev->EV[0] == ECurr->EV[1] )
        VCurr = EPrev->EV[0];
      else
        VCurr = EPrev->EV[1];
      
      #ifdef ROBUST
         check( (ECurr->EV[0]!=VCurr && ECurr->EV[1]!=VCurr),
	      "TDestroyDelaunay::RetriangulateInfluenceRegion(), <2> inconsistency detected" ); 
      #endif
      
      
      // set VE relations of VCurr to point to ECurr/EPrev, if they point
      // to internal edges of the influence region to be deleted

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
      
      // else ... !E0Del && !E1Del... leave VE relations unchanged

      if (CheckIsolatedPoint(VCurr))
         cout << "ExtRIRIPrepare (2)" << endl;
      
      // go to next edge
      
      EPrev = ECurr;
      IRegnIter.GoNext();
      
      
   } // end ... while ( ! EndOfList() )   
}


//
// Equal to TDestroyDelaunay::ExtRIRInitialTrg with the difference that
// it marks every newly created edge as NEW_EDGE.
//

void TDecCDT::ExtRIRInitialTrg( TDoubleList<PTEdge> & InflRegnBorder,
                                      TDoubleList<PTEdge> & SwapEdgeQueue   )
{

    int e;

    #ifdef DEBUG
       DEBUG << "TDestroyDelaunay::RIRInitialTrg()" << endl;
    #endif
    
    PTTriangle NewTriangle;
    
    //
    // Make a copy of list InflRegnBorder. This second list is needed 
    // because we make modifications to it, while we need to have the
    // content of InflRegnBorder intact when we call function
    // RIRDelOptTrg().
    // Use also a different mark, INFL_BORDER_AUX, to mark edges in 
    // this second list.
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
    // Main loop to re-triangulate the star-shaped region whode
    // boundary is formed by the edges in InflRegnBorder/InflRegnAux
    //
    
    TDoubleListIterator<PTEdge> IAuxIter( &InflRegnAux );
    
    IAuxIter.Restart();
    
    // current and next edge on InflRegnAux
    PTEdge ECurr = IAuxIter.Current()->object;
    PTEdge ENext;

    // next three vertices (in CCW order) on InflRegnAux
    PTVertex v0, v1, v2;

    while( InflRegnAux.Lenght() > 3 )
    {
 
       IAuxIter.GoNext();
       if ( IAuxIter.EndOfList() ) IAuxIter.Restart(); 
       // manage list as if it were circular
 
       // next edge (after ECurr)
       ENext = IAuxIter.Current()->object;
       
       
       // v0, v1, v2 are the endpoints of ECurr and ENext (v1 is common)
       // sorted counterclockwise
              
       v0 = ECurr->EV[0];
       v1 = ECurr->EV[1];

       if ( ENext->EV[0] == v0 || ENext->EV[1] == v0 )
       {
          // v0 is the common vertex, swap v0 and v1
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
           
   	   // create triangle v0, v1, v2, the boundary edges are
	   // ECurr, ENext, and a new edge (v2, v0)

	   PTEdge NewEdge = new TEdge( v2, v0 );
	   check( (NewEdge == NULL), "TDestroyDelaunay::RIRInitialTrg(), insufficient memory" );
	   
	   NewTriangle = new TTriangle( ECurr, ENext, NewEdge );
	   	   
	   // since NewEdge = (v2, v0), NewTriangle is to the left 
           // of NewEdge
	   
	   NewEdge->ET[0] = NewTriangle;
	   NewEdge->ET[1] = NULL;


	   // enqueue new edge for RIRDelOptTrg()
	   
	   SwapEdgeQueue.AddHead( NewEdge );
	   NewEdge->Mark( SWAP_EDGE_QUEUE );

      // mark new edge as NEW_EDGE to distinguish it, during optimization,
      // from original edges external to the region of influence which
      // are reached by propagation of optimization

      NewEdge->Mark( NEW_EDGE );
	   
	   // ET relations of  ECurr and ENext from the "internal" side
	   // to the remaining part (not yet triangulated) of the region 
	   // of influence are set to NewTriangle. Recall that 
           // ECurr = (v0, v1) and ENext = (v1, v2).

	   if ( ECurr->EV[0] == v0 && ECurr->EV[1] == v1 )
              // NewTriangle to the left of ECurr
	      ECurr->ET[0] = NewTriangle;
	   else 
	   if ( ECurr->EV[0] == v1 && ECurr->EV[1] == v0 )
              // NewTriangle to the right of ECurr
	      ECurr->ET[1] = NewTriangle;
	   else
	      error( "TDestroyDelaunay::RIRInitialTrg(), <> inconsistency detected" );
	      
	      
	   // same thing for ENext
	   
	   if ( ENext->EV[0] == v1 && ENext->EV[1] == v2 ) 
              // NewTriangle to the left of ECurr
	      ENext->ET[0] = NewTriangle;
	   else 
	   if ( ENext->EV[0] == v2 && ENext->EV[1] == v1 )
              // NewTriangle to the right of ECurr
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
	      
		    
	   
	   // Remove ECurr, ENext from list InflRegnAux (and remove mark
	   // INFL_BORDER_AUX), and replace them with NewEdge (and mark it
           // as INFL_BORDER_AUX). Recall that IAuxIter points to the
	   // node containing ENext. 
	   // Since we have managed the list as a circular list, we have
           // to pay attention when  ECurr is the last element and
	   // ENext is the first one.
	   
	   ECurr->UnMark( INFL_BORDER_AUX );
	   ENext->UnMark( INFL_BORDER_AUX );
	   
	   if ( IAuxIter.StartOfList() ) // IAuxIter points to first of list
	   {
	      // ECurr is the last element of the list,  ENext the first one
	      
	      InflRegnAux.RemoveHead(); // remove ENext
	      InflRegnAux.RemoveLast(); // remove ECurr
	      
	      NewEdge->Mark( INFL_BORDER_AUX ); // add NewEdge
	      InflRegnAux.AddHead( NewEdge );
	      
	      IAuxIter.Restart();
	         
	   }
	   else // ! IAuxIter.BeginOfList()
	   {
	      // remove ECurr
	      InflRegnAux.RemoveBefore( IAuxIter.Current() ); 

              // put NewEdge after ENext, and mark it
	      NewEdge->Mark( INFL_BORDER_AUX );
	      InflRegnAux.AddAfter( IAuxIter.Current(), NewEdge );
	      
	      // go to next edge (NewEdge)	      
	      IAuxIter.GoNext();
	      
	      // delete ENext (now the one before NewCurrent())
	      InflRegnAux.RemoveBefore( IAuxIter.Current() );	      
	   }
	   	   
       }
       
              
       ECurr = IAuxIter.Current()->object;
	          
    } // end ...while( InflRegnAux.Lenght() > 3 )
    

    
    //
    // now there are three edges left in lust InflRegnAux, create last
    // triangle adjacent to such edges
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
    
    
    // adjust ET relations of the three edges

    PTVertex v[3];
        
    NewTriangle->GetTV( v[0], v[1], v[2] );
    
    
    // based on the implementation of GetTV, V[i] and V[i+1 % 3] are the
    // vertices of E[i] in counterclockwise order
    
    for( e=0; e<3; e++ )
    {
       if ( E[e]->EV[0] == v[e] && E[e]->EV[1] == v[(e+1)%3] )
          E[e]->ET[0] = NewTriangle; // left of E[e]
       else
       if ( E[e]->EV[1] == v[e] && E[e]->EV[0] == v[(e+1)%3] )
          E[e]->ET[1] = NewTriangle; // right of E[e]
       else
          error( "TDestroyDelaunay::RIRInitialTrg(), <> inconsistency detected" );
    }
    
    // add last triangle
            
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


//
// Called also on edges outside the influence region:
//  - if OldT0/1 are original, do not call Detach but simply mark them as
//    TO_DELETE;
//  - intermediate edges are marked as NEW_EDGE 
//  - adjust VE relations of the two endpoints of edge E removed from
//    the triangulation:
//    E may be outside the region of influence and VE* of one of its
//    endpoints may point to E; thus we set such pointer to one edge
//    of the quadrilateral.
//
PTEdge TDecCDT::ExtSwapEdge( PTEdge &E, TDoubleList<TEdge> &OrigEdgList )
{

    #ifdef DEBUG
       DEBUG << "TDecCDT::ExtSwapEdge(E" << E->EID << ")" << endl;
    #endif
    
    PTTriangle OldT0 = E->ET[0];
    PTTriangle OldT1 = E->ET[1];
   
    #ifdef ROBUST
       check( (OldT0 == NULL || OldT1 == NULL), "TDestroyDelaunay::SwapEdge() <1> inconsistency detected" );
    #endif

    PTEdge E0a, E0b, E1a, E1b;
    GetQuadBorder( E, E0a, E0b, E1a, E1b );    
    
    // vo0 and vo1 the vertices of  OldT0 and OldT1, respectively, that are
    // not endpoints of E (vo0 is common to E0a and E0b, while
    // vo1 is common to E1a and  E1b)
        
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

    PTEdge NewE = new TEdge( vo0, vo1 ); // ...NB: diretto da vo0 a vo1
    check( (NewE == NULL), "TDecCDT::ExtSwapEdge(), insufficient memory" );


    // mark the new edge as NEW_EDGE in order to optimize Ext
    //
    NewE->Mark( NEW_EDGE );

    PTTriangle NewT0 = new TTriangle( NewE, E1b, E0a );
    PTTriangle NewT1 = new TTriangle( NewE, E0b, E1a );
    check( (NewT0 == NULL || NewT1 == NULL), "TDestroyDelaunay::SwapEdge(), insufficient memory" );
    
    // adjust Et relations of  NewE, E0a, E0b, E1a, E1b
    //
    // Some of E0a, ... E1b may be an original edge; such edge may then 
    // be involved in the propagation of optimization and be marked as
    // TO_DELETE.
    // The problem is that here we set its ET to a new, non-original
    // and intermediate triangle, that could later be eliminated,
    // while usually an edge is marked as TO_DELETE by function CalcInflRegn
    // (that marks only edges inside the influence region and adjacent to
    // original triangles marked as TO_DELETE).
    // This creates problems in TDelaunayBase::DeleteInfluenceRegion that, 
    // when considers an edge E marked as TO_DELETE and decides whether
    // to include the adjacent triangles of E among the triangles
    // on which to call DetachTriangle, does not expect that E may have
    // one of the two triangles in its ET relation set to NULL.
    // This problem is solved by function TDecCDT::ExtDelBaseDelInflRegn(),
    // equal to TDelaunayBase::DeleteInfluenceRegion, but solving this
    // problem; it is called by ExtDeleteInfluenceRegion instead of
    // TDelaunayBase::DeleteInfluenceRegion.
    //
    
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

   // We must djust the VE relations of the 2 vertices of edge E that is
   // removed from the triangulation:
   // E may be outside the region of influence and thus the VE relation 
   // of one of its vertices may point to E.
   // Thus we set such relation to point to one edge of the quadrilateral.
   //

   // first vertex of E is incident in edges E0b and E1a
   if( E->EV[0]->VE[0] == E ) 
      // if VE[0] of first vertex of E points to E, set VE* of such
      // vertex to E0b or E1a 
   { 
      E->EV[0]->VE[1] == E0b ? E->EV[0]->VE[0] = E1a : E->EV[0]->VE[0] = E0b;
   }
   if( E->EV[0]->VE[1] == E )
      // if VE[1] of first vertex of E points to E, set VE* of such
      // vertex to E0b or E1a 
   {     
      E->EV[0]->VE[0] == E0b ? E->EV[0]->VE[1] = E1a : E->EV[0]->VE[1] = E0b;
   }

if (CheckIsolatedPoint(E->EV[0]))
   cout << "ExtSwapEdge (1)" << endl;

   //
   // second vertex of E is incident in edges E1b and E0a
   if( E->EV[1]->VE[0] == E )
      // if VE[0] of second vertex of E points to E, set VE* of such
      // vertex to E1b or E0a
   {  
      E->EV[1]->VE[1] == E1b ? E->EV[1]->VE[0] = E0a : E->EV[1]->VE[0] = E1b;
   }
   if( E->EV[1]->VE[1] == E )
      // if VE[1] of second vertex of E points to E, set VE* of such
      // vertex to E1b or E0a
   { 
      E->EV[1]->VE[0] == E1b ? E->EV[1]->VE[1] = E0a : E->EV[1]->VE[1] = E1b;
   }

if (CheckIsolatedPoint(E->EV[1]))
   cout << "ExtSwapEdge (2)" << endl;

   // edge E may be original, in this case we must not delete it directly
   //
   if( E->Marked( NEW_EDGE ) ) // E is not an original edge
   {
      DetachEdge( E );
      delete(E);
      E = NULL;
   }
   else   // E is an original edge
   {
      E->Mark( TO_DELETE );
      ReconstructET( E, OrigEdgList );
      // restore original E->ET relation and remove from OrigEdgList
      // the copy of E
      E->UnMark( COPIED );
   }
    
   // triangles OldT0/OldT1 may be original, in this case we must not 
   // delete them directly
   //
   if( OldT0->Marked( NEW_TRIANGLE ) )
   {
      DetachTriangle( OldT0 );
      delete( OldT0 );
   }
   else
      OldT0->Mark( TO_DELETE );

   if( OldT1->Marked( NEW_TRIANGLE ) )
   {
      DetachTriangle( OldT1 );
      delete( OldT1 );
   }
   else
      OldT1->Mark( TO_DELETE );

   
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


//
// Equal to TDestroyDelaunay::RIRDelOptTrg with the following differences:
//    - to swap a diagonal call ExtSwapEdge instead of SwapEdge;
//    - insert in SwapEdgeQueue also edges outside the region of influence
//      provided that they are not constraints and the are not on the 
//      convex hull
//
void TDecCDT::ExtRIRDelOptTrg( TDoubleList<PTEdge> & SwapEdgeQueue,TDoubleList<TEdge> &OrigEdgList )
{
   int e;
   // PTEdge pe;

   #ifdef DEBUG
      DEBUG << "TDecCDT::ExtRIRDelOptTrg()" << endl;
   #endif

   PTEdge ECurr, NewEdge, E[4];

   while( !SwapEdgeQueue.IsEmpty() )
   {
      ECurr = SwapEdgeQueue.RemoveHead();
      ECurr->UnMark( SWAP_EDGE_QUEUE );
   
      if ( EdgeToSwap( ECurr ) )
      {
         // copy in OrigEdgList the boundary edges of the quadrilateral
         // formed by the two triangles adjacent to ECurr that will be put
         // in SwapEdgeQueue after the diagonal swap and that are original
         // and not yet present in OrigEdgList (marked as COPIED)
          
         GetQuadBorder( ECurr, E[0], E[1], E[2], E[3] );
       
         for( e=0; e<4; e++ )
         {
            if (    !E[e]->Marked( SWAP_EDGE_QUEUE ) && !E[e]->Marked( CONSTRAINED ) && !E[e]->OnConvexHull()
                 && !E[e]->Marked(NEW_EDGE) && !E[e]->Marked(COPIED) )
            {
               // pe = new TEdge(*E[e]);
               OrigEdgList.AddHead( *E[e] );
               E[e]->Mark(COPIED);
            }
         }

         // swap edge
       
         NewEdge = ExtSwapEdge( ECurr, OrigEdgList );
       
         // enqueue the boundary edges of the quadrilateral formed
         // by the two triangles adjacent to NewEdge
          
         GetQuadBorder( NewEdge, E[0], E[1], E[2], E[3] ); // 2DO: it should not be necessary
       
         for( e=0; e<4; e++ )
         {
            if ( !E[e]->Marked( SWAP_EDGE_QUEUE ) && !E[e]->Marked( CONSTRAINED )
                                            && !E[e]->OnConvexHull() )
            {
               E[e]->Mark( SWAP_EDGE_QUEUE );
               SwapEdgeQueue.AddTail( E[e] );
            }
         }
      } // end ...if( EdgeToSwap )
   } // end ...while( !SwapEdgeQueue.IsEmpty() )

   // End of optimization. Empty list OrigEdgList that contains a copy
   // of the original edges considered but not swapped.
   // ATTENTION: the COPIED mark must be removed from such edges.
   //            This is done by function ExtDelBaseDelInflRegn called by
   //            ExtDeleteInfluenceRegion instead of
   //            TDelaunayBase::DeleteInfluenceRegion() .
   //
   #ifdef ROBUST
      TDoubleListIterator<TEdge> OrigEdgIter( &OrigEdgList );
      OrigEdgIter.Restart();
      PTEdge edg;
      while ( ! OrigEdgIter.EndOfList() )
      {
         edg = &( OrigEdgIter.Current()->object );
         if( edg->Marked( TO_DELETE ) || edg->Marked( NEW_EDGE ) )
            Pause("BUG: list OrigEdgList contains a non-original edge or an adge not marked as TO_DELETE");

         OrigEdgIter.GoNext();
      }
   #endif
   OrigEdgList.ClearList();
}


// 
// Similar to TDecimDelaunay::DeleteInfluenceRegion with the following 
// differences:
//    - to recheck vertices involved in the re-triangulation do not use
//      list InflRegnBorder (that contains only those of the influence 
//      region) but use function RecheckOptimizedRegion.
//    - call ExtMT_AddComponent that adds MT_AddComponent the following 
//      operations: delete mark RECHECKED from vertices and triangles
//      marked by function RecheckOptimizedRegion; delete mark NEW_EDGE
//      from final new edges built during re-triangulation.
//
void TDecCDT::ExtDeleteInfluenceRegion( TDoubleList<PTEdge> & InflRegnBorder )
{

   #ifdef DEBUG
     DEBUG << "TDecimDelaunay::DeleteInfluenceRegion()" << endl;
   #endif
    
   //
   // delete old triangles
   //
   
   if ( ! InitialPhase ) MT_KillInterference();
    

   // recheck vertices involved in the re-triangulation
   //
   RecheckOptimizedRegion();

   // Traverse triangles marked as TO_DELETE through adjacency navigation
   // starting from FirstTrgToDel.
   // Remove mark TO_DELETE from triangles and mark edges as COPIED.
   // Insert vertices present in the PointList's of triangles and edges
   // in DetachedPoints and eliminate them from the point lists by using
   // DetachTriangle/Edge.
   //
   //  TDelaunayBase::DeleteInfluenceRegion();
   ExtDelBaseDelInflRegn();

   //
   // Reposition points "Detached" from deleted triangles into new
   // triangles, in such a way that we can compute errors correctly
   // in MT_AddComponent().
   //
   
   RepositionPoint( VertexToRemove );
   
   while( ! DetachedPoints.IsEmpty() )
   {
       PTPoint p = DetachedPoints.RemoveHead();
       RepositionPoint( p );       
   } 

   //
   // Now add new triangles (those marked as NEW_TRIANGLE) and unmark them.
   // In addition to MT_AddComponent we do:
   // - delete mark NEW_EDGE from marked edges belonging to TE(t)
   // - delete mark RECHEKED from marked vertices belonging to TV(t).
   // - delete mark RECHEKED from visited triangles
   //   
   if ( ! InitialPhase ) ExtMT_AddComponent();

   // delete INFL_BORDER mark from edges in InflRegnBorder and free 
   // list InflRegnBorder
   //
   Del_InflRegnBorder( InflRegnBorder );

}


void TDecCDT::RecheckOptimizedRegion()
{

   #ifdef DEBUG
     DEBUG << endl << "MT: delete V" << VertexToRemove->VID 
          << ": " << VertexToRemove->x << ", " << VertexToRemove->y 
          << " => ERR : " << VertexToRemove->Error << endl;
   #endif


   if ( FirstTrgToDel != NULL )
   {
   
      #ifdef ROBUST
         check( (! FirstTrgToDel->Marked( TO_DELETE ) ),
            "TDecCDT::RecheckOptimizedRegion(), inconsistency detected" );
      #endif

      TDoubleList<PTTriangle> Triangles;
      
      FirstTrgToDel->Mark( RECHECKED );
      Triangles.AddHead( FirstTrgToDel );
      
      while( !Triangles.IsEmpty() )
      {
          PTTriangle T = Triangles.RemoveHead();
          PTTriangle TT[3];
          PTVertex TV[3];
      
          // cannot remove mark TO_DELETE, it is needed for
          // TBuildDelaunay::DeleteInfluenceRegion()
      
          T->GetTT( TT[0], TT[1], TT[2] );
          T->GetTV( TV[0], TV[1], TV[2] );
      
          //
          // operations needed for the extension of optimization
          //
          int i;
          for( i = 0; i <= 2; i++ )
          {
             if( ! TV[i]->Marked( RECHECKED ) )
             {
                if( ElimVtxTree.IsIn( TV[i] ) )
                   ElimVtxTree.Remove( TV[i] );

                if( ReCheckVertex( TV[i] ) )
                   ElimVtxTree.Insert( TV[i] );

                TV[i]->Mark( RECHECKED );
             }

          }

          // enqueue triangles adjacent to T and not yet visited
          //
          for( int t=0; t<3; t++ )
             if ( TT[t] != NULL && TT[t]->Marked( TO_DELETE ) 
                  && ! TT[t]->Marked( RECHECKED ) )
                  {  
                     TT[t]->Mark( RECHECKED );
                     Triangles.AddTail( TT[t] );
                  }		     
      }
      
   }
   // VertexToRemove may have been re-inserted in ElimVtxTree: remove it
   //
   if( ElimVtxTree.IsIn( VertexToRemove ) )
      ElimVtxTree.Remove( VertexToRemove );
}


// Equal to MT_AddComponent but add what follows,
// for each visited triangle t:
// - delete mark NEW_EDGE from marked edges belonging to TE(t)
// - delete mark RECHEKED from marked vertices belonging to TV(t)
// - delete mark RECHEKED from any visited triangle t
//   
void TDecCDT::ExtMT_AddComponent()
{

   #ifdef MT_DEBUG
      MT_DEBUG << "TDecCDT::ExtMT_AddComponent() record retriangulation after deleting V" << VertexToRemove->VID << endl;
   #endif

   #ifdef ROBUST
      check( (FirstTriangle == NULL || !FirstTriangle->Marked( NEW_TRIANGLE )),
             "TDecCDT::ExtMT_AddComponent(), inconsistency detected" );
   #endif
     
   TDoubleList<PTTriangle> Triangles;
     
   FirstTriangle->UnMark( NEW_TRIANGLE );
   Triangles.AddHead( FirstTriangle );
   
   while ( ! Triangles.IsEmpty() )
   {
      PTTriangle T = Triangles.RemoveHead();
      PTTriangle TT[3];
      PTVertex   TV[3];

      T->GetTV( TV[0], TV[1], TV[2] );
      T->GetTT( TT[0], TT[1], TT[2] );

      // operations of version Ext
      //

      T->UnMark( RECHECKED );

      int i;
      for( i = 0; i <= 2; i++ )
      {
         TV[i]->UnMark( RECHECKED );
         T->TE[i]->UnMark( NEW_EDGE );
      }

      double Error = 0.0;  
      Error = Error;  //  needed to avoid warning
      if ( ! T->PointList.IsEmpty() ) Error = T->PointList.GetHead()->Error;
      
      MT->MakeTriangle( T );

      // enqueue adjacent triangles that are not yet in queue
      //
      for( int t=0; t<3; t++ )
         if ( TT[t] != NULL && TT[t]->Marked( NEW_TRIANGLE ) )
         {
            TT[t]->UnMark( NEW_TRIANGLE );
            Triangles.AddTail( TT[t] );
         }
      
   }
      
   MT->MeshOk();
  
}



// Remove mark INFL_BORDER from edges in InflRegnBorder and free such list
//
void TDecCDT::Del_InflRegnBorder( TDoubleList<PTEdge> & InflRegnBorder )
{
   TDoubleListIterator<PTEdge> IRegnIter( &InflRegnBorder );
   IRegnIter.Restart();
   while ( ! IRegnIter.EndOfList() )
   {
       IRegnIter.Current()->object->UnMark( INFL_BORDER );
       IRegnIter.GoNext();
   }
      
   // empty InflRegnBorder, it is no longer needed
   
   InflRegnBorder.ClearList();
}


// Traverse original influence region, insert its points in 
// DetachedPoints, and delete its edges and triangles.
//
// Traverse triangles marked as TO_DELETE through adjacency navigation 
// starting from FirstTrgToDel.
// Remove mark TO_DELETE from triangles and mark edges as COPIED.
// Insert vertices present in the PointList's of triangles and edges
// into DetachedPoints and remove such lists by calling DetachTriangle/Edge
//
    // Some of E0a, ... E1b may be an original edge; such edge may then 
    // be involved in the propagation of optimization and be marked as
    // TO_DELETE.
    // The problem is that here we set its ET to a new, non-original
    // and intermediate triangle, that could later be eliminated,
    // while usually an edge is marked as TO_DELETE by function CalcInflRegn
    // (that marks only edges inside the influence region and adjacent to
    // original triangles marked as TO_DELETE).
    // This creates problems in TDelaunayBase::DeleteInfluenceRegion that, 
    // when considers an edge E marked as TO_DELETE and decides whether
    // to include the adjacent triangles of E among the triangles
    // on which to call DetachTriangle, does not expect that E may have
    // one of the two triangles in its ET relation set to NULL.
    // This problem is solved by function TDecCDT::ExtDelBaseDelInflRegn(),
    // equal to TDelaunayBase::DeleteInfluenceRegion, but solving this
    // problem; it is called by ExtDeleteInfluenceRegion instead of
    // TDelaunayBase::DeleteInfluenceRegion.
//
// NOTE: indeed, the solution consists of reconstructing the region 
// marked as TO_DELETE in such a way that all (internal) edges of such
// region have their Et relation correctly set, i.e., containing the 
// two triangles of the region (original edges marked TO_DELETE)
// adjacent to it. Such reconstruction must be done at the end of
// optimization and before visiting the region.
// Reconstruction can be done in two ways:
//  - record in a list a copy of all original edges inserted in
//    SwapEdgeQueue; when we mark an original edge E as TO_DELETE
//    inside function ExtSwapEdge, we extract the copy of E from the 
//    list and restore original adjacency relations by copying then
//    from the copy.
//  - after the end of optimization, visit the region marked as TO_DELETE 
//    and reconstruct ET relations for the modified edges.
//
void TDecCDT::ExtDelBaseDelInflRegn()
{
   #ifdef DEBUG
      DEBUG << "\nTDecCDT::ExtDelBaseDelInflRegn( "
            << "FirstTrgToDel = ";
      if ( FirstTrgToDel == NULL ) DEBUG << "NULL"; else 
      DEBUG << "T" << FirstTrgToDel->TID;
      DEBUG << " )" << endl;
   #endif // DEBUG

   #ifdef DEBUG15
      DEBUG15 << "\nTDecCDT::ExtDelBaseDelInflRegn : visit the following entities:\n";
   #endif

   // if no triangles to be deleted
   if ( FirstTrgToDel == NULL ) return; 
     
   TDoubleList<PTTriangle> TrgToDel;
   PTTriangle CurrTrg, NextTrg;
   PTEdge NextEdg;

   TrgToDel.AddHead( FirstTrgToDel );
     
   while ( ! TrgToDel.IsEmpty() )
   {
      CurrTrg = TrgToDel.RemoveHead();

      #ifdef DEBUG15
         DEBUG15 << "T" << CurrTrg->TID << endl;
      #endif
      
      #ifdef ROBUST
         check( (CurrTrg == NULL), "TDelaunayBase::DeleteInfluenceRegion(), <1> this should not happen");
      #endif
        
      for( int e=0; e<3; e++ )
      {      
         NextEdg = CurrTrg->TE[e];

         //
         // edge may have been deleted in previous iterations of the
         // while loop, we need to check
         //
         
         if ( NextEdg != NULL ) 
         {
            #ifdef DEBUG15
               DEBUG15 << *NextEdg << endl;
            #endif
      
            // smarco la marctura COPIED se presente
            if( NextEdg->Marked(COPIED) )
            {
               #ifdef DEBUG15
                  DEBUG15 << "Unmark copied " << *NextEdg << endl;
               #endif
               NextEdg->UnMark(COPIED);
            }

            if ( NextEdg->Marked(TO_DELETE) )
            {
               // NextTrg = ( (NextEdg->ET[0] == CurrTrg || NextEdg->ET[0] == NULL) ? NextEdg->ET[1] : NextEdg->ET[0] );
               if( NextEdg->OnConvexHull() )
                  NextTrg = NULL;
               else
               {
                  check( NextEdg->ET[0] == NULL || NextEdg->ET[1] == NULL, "ExtDelBaseDelInflRegn: <a> inconsistency detected\n" );
                  NextTrg = ( NextEdg->ET[0] == CurrTrg ? NextEdg->ET[1] : NextEdg->ET[0] ); 
               }    
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
      // remove CurrTrg from the triangulation, destructor of TTriangle
      // sets to NULL the pointers to CurrTrg contained in edges adjacent
      // to it
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
      DEBUG << "exit TDecCDT::ExtDelBaseDelInflRegn" << endl;
   #endif 
}


//
// Restore original E->ET relation and remove the copy of E from
// OrigEdgList.
// Given an original edge E that must be removed from the triangulation
// (marked as TO_DELETE) and the list OrigEdgList containing the copies
// of original edges involved in the optimization,
// extract from the list the copy of E and reconstruct from it the
// ET relation of E (by copying from the copy).
//
void TDecCDT::ReconstructET( PTEdge E, TDoubleList<TEdge> &OrigEdgList )
{
   TDoubleListIterator<TEdge> OrigEdgIter( &OrigEdgList );
   OrigEdgIter.Restart();
   PTEdge edg;
   while ( ! OrigEdgIter.EndOfList() )
   {
      edg = &( OrigEdgIter.Current()->object );
      if( edg->EID == E->EID )
      {
         E->ET[0] = edg->ET[0];
         E->ET[1] = edg->ET[1];
         OrigEdgList.RemoveAt( OrigEdgIter.Current() );
         return;
      }
      OrigEdgIter.GoNext();
   }
   // if we arrive here, we have not found the edge
   cerr << "TDecCDT::ReconstructET(" << *E << " , OrigEdgList): edge not found in list\n";
   Pause("");
}

//----------------------------------------------------------------------------
//
//  void TDecCDT ::ConvertData(int *vNum, int *tNum, int *eNum,
//                             float **vData, int **tData, int **eData)
//

void TDecCDT::ConvertData(int *vNum, int *tNum, int *eNum,
                          float **vData, int **tData, int **eData)
{
    int iv, it, ie, i;

     #ifdef DEBUG
      DEBUG << "\nTDecCDT::ConvertData(" << outfname << ")" << endl;
     #endif // DEBUG

    //
    // create arrays to contain output vertices, triangles and constraints
    //
    
    // PAOLA: number or really existing vertices
    int nVrt = 0;

    PTVertex *VtxArray;   
    PTTriangle *TrgArray;
    PTEdge * EdgArray;

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
    EdgArray = new PTEdge[ nConstrInTRI ];
    for( ie=0; ie<nConstrInTRI; ie++ ) EdgArray[ie] = NULL;

    check( ( VtxArray == NULL || TrgArray == NULL ),
           "TDecCDT::ConvertData(), insufficient memory" );
	  
    for( iv=0; iv<nPts; iv++ ) VtxArray[iv] = NULL;
    for( it=0; it<nTrg; it++ ) TrgArray[it] = NULL;
   
    TDoubleList<PTTriangle> Triangles;
   
    check( (FirstTriangle==NULL), "TDecCDT::ConvertData(), No triangles?");
   
    Triangles.AddHead( FirstTriangle );
    FirstTriangle->Mark( VISITED );
   
    PTVertex v[3];
   
    PTTriangle NextTrg[3];
    PTTriangle CurTrg;
   
    it = 0;
    iv = 0;
    ie = 0;

    while( !Triangles.IsEmpty() )
    {
       
       CurTrg = Triangles.RemoveHead();
       
       TrgArray[it++] = CurTrg;

       // ...vertices of CurTrg
      
       CurTrg->GetTV( v[0], v[1], v[2] );
      
       #ifdef ROBUST
          check( (v[0] == NULL || v[1] == NULL || v[2] == NULL),
              "TDecCDT::ConvertData(), inconsistency detected");
       #endif
      
       // ...adjacent triangles
       
       CurTrg->GetTT( NextTrg[0], NextTrg[1], NextTrg[2] );

       // ...edges of CurTrg
       int i;
       for( i=0; i<3; i++ )
       {
          if( ! CurTrg->TE[i]->Marked( VISITED ) && CurTrg->TE[i]->Marked( CONSTRAINED )) // ogni edge e' comune a due triangoli quindi
          {        // it may have been already visited
            #ifdef ROBUST
               check( ie >= nConstrInTRI, "TDecCDT::ConvertData(), constraints inconsistency detected" );
            #endif

            EdgArray[ie] = CurTrg->TE[i];
            EdgArray[ie]->Mark( VISITED );
            ie++;
          }
       }
       
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
    
if (it != nTrg) printf("it = %d diverso da nTrg = %d\n",it, nTrg); 
    check( it != nTrg, "TDecCDT::ConvertData(), <1> inconsistency detected");
    
    
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
              "TDecCDT::ConvertData(), error in vertex number");
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
       
       // PAOLA: eliminata la parte che serviva a stampare PS
       // e stampa solo i vertici non nulli.
      
       if ( CV != NULL )
       {  (*vData)[3*i] = CV->x;
          (*vData)[3*i+1] = CV->y;
          (*vData)[3*i+2] = CV->z;
          i++;
       }      
       if ( iv % 1000 == 0 ) cerr << "output " << iv << " vertices\r";

    }
    
    //
    // ...then triangles
    //
    
    for( it=0; it<nTrg; it++ )
    {
       PTTriangle CT = TrgArray[it];

       check( (CT == NULL), "TDecCDT::ConvertData(), <2> inconsistency detected");

       CT->GetTV( v[0], v[1], v[2] );
       
       (*tData)[3*it] = v[0]->VID;
       (*tData)[3*it+1] = v[1]->VID;
       (*tData)[3*it+2] = v[2]->VID;
       
       if ( it % 1000 == 0 ) cerr << "output " << it << " triangles\r";

       CT->UnMark( VISITED ); // smarchiamo i triangoli
    }

    cerr << "output " << nTrg << " triangles" << endl;

    //
    // ...finally edges
    //

    (*eNum) = nConstrInTRI;
    *eData = (int*)malloc(2*nConstrInTRI*sizeof(int));

    check( ie != nConstrInTRI, "TDecCDT::ConvertData(), <3> inconsistency detected");

    for( ie=0; ie<nConstrInTRI; ie++ )
    {
       int vIndx[2];

       vIndx[0] = EdgArray[ie]->EV[0]->VID;
       vIndx[1] = EdgArray[ie]->EV[1]->VID;

       (*eData)[2*ie] = vIndx[0];
       (*eData)[2*ie+1] = vIndx[1];
       
       if ( ie % 1000 == 0 ) cerr << "output " << ie << " triangles\r";
       EdgArray[ie]->UnMark( VISITED ); // smarchiamo i lati
    }
    
    cerr << "output " << nConstrInTRI << " constraints of " << nConstrInFile << "\n";       
    
}
