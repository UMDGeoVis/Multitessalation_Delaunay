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

// ---------------------------------------------------------------------------
//
//   file   : utils.cpp
//   author : Alessio Calcagno and Christian Melchiorre
// 
//   Some routins for utility purposes (for instance, swap two variables,
//   etc...), used in the Delaunay triangulation programs.
//

#include "defs.h"
#include "ttriang.h"
#include "utils.h"
#include "geom.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

///////////////////////////////////////////////////////////////////////////
//
//  Methods to navigate the triangulation and check consistency of
//  topological relations.
//
///////////////////////////////////////////////////////////////////////////

//
// Test if topological relations are defined.
// 

boolean Is_EV_Defined ( PTEdge e ) {

   if( e->EV[0] != NULL && e->EV[1] != NULL )  return TRUE;
   else return FALSE;
}

boolean Is_ET_Defined ( PTEdge e ) {

   if( e->ET[0] != NULL && e->ET[1] != NULL )  return TRUE;
   else return FALSE;
}

boolean Is_v_in_e( PTVertex V, PTEdge E)
{
   if( E->EV[0] == V || E->EV[1] == V ) return TRUE;
   else                                 return FALSE;
}

boolean Is_v_in_t( PTVertex V, PTTriangle T)
{
   if( Is_v_in_e( V, T->TE[0] ) || Is_v_in_e( V, T->TE[1] ) )  return TRUE;
   else                                                        return FALSE;
}

//
//  NAVIGATION around a vertex or vertex-based navigation
//

PTVertex Other_v_in_e ( PTEdge e, PTVertex v ) {

  return  e->EV[0] == v ?  e->EV[1] :  e->EV[0] ;  
}


PTEdge Other_e_in_VE ( PTEdge E, PTVertex V )
{
  return  V->VE[0] == E ?  V->VE[1] :  V->VE[0] ;  
}


PTEdge Left_e_of_t_by_v( PTTriangle T, PTVertex V)
{
   for( int i = 0; i < 3; i++ )
      if( Is_v_in_e( V, T->TE[i] ) && Is_v_in_e( V, T->TE[(i+1)%3] ) ) 
         // in VC6 a wrong address is passed instead of
         // T->TE[(i+1)%3] when i = 2 (don't know why)
         return T->TE[i];                                             
   return 0; // otherwise in VC5: not all control paths return a value
}

PTEdge Right_e_of_t_by_v( PTTriangle T, PTVertex V)
{
   for( int i = 0; i < 3; i++ )
      if( Is_v_in_e( V, T->TE[i] ) && Is_v_in_e( V, T->TE[(i+1)%3] ) )
         return T->TE[(i+1)%3];
   return 0;
}

PTEdge  Prev_CCW_e_of_e_by_v( PTEdge E, PTVertex V )
{
   return  Next_CCW_e_of_t_by_e( Prev_CCW_t_of_e_by_v( E, V ), E );
}

PTEdge  Opposite_e_of_t_by_v( PTTriangle T, PTVertex V )
{
   for( int i = 0; i < 3 ; i++ )
      if( ! Is_v_in_e( V, T->TE[i] )  ) return T->TE[i];
   return 0; // otherwise in VC5: not all control paths return a value
}

PTTriangle  Opposite_t_of_t_by_e( PTTriangle T, PTEdge E )
{
   #ifdef ROBUST
      if( T->TE[0] != E  &&  T->TE[1] != E  &&  T->TE[2] != E )
      {  cerr << "ERROR: Opposite_t_of_t_by_e( " << *T << " , " << *E << " ), edge not in TE of triangle\n";
         exit(-1); }
   #endif

   if( E->ET[0] == T )  return E->ET[1];
   else                 return E->ET[0]; 
}

PTEdge  Next_CCW_e_of_e_by_v( PTEdge E, PTVertex V )
{
   PTTriangle nextTrg = Next_CCW_t_of_e_by_v( E, V );

   if( nextTrg == NULL )  return NULL;
   else                   return Prev_CCW_e_of_t_by_e( nextTrg, E );
}

PTTriangle  Next_CCW_t_of_e_by_v( PTEdge E, PTVertex V )
{
   if( E->EV[0] == V )  return E->ET[0];
   else                 return E->ET[1];
}

PTTriangle  Prev_CCW_t_of_e_by_v( PTEdge E, PTVertex V )
{
   if( E->EV[0] == V )  return E->ET[1];
   else                 return E->ET[0];
}

PTEdge  Prev_CCW_e_of_t_by_e( PTTriangle T, PTEdge E )
{
   for( int i = 0; i < 3; i++ )
   {
      if( T->TE[i] == E )  return T->TE[(i+2)%3]; // in VC5 non so perche' ma ritorna un indirizzo sbagliato
   }
   return 0; // otherwise in VC5: not all control paths return a value
}

PTEdge  Next_CCW_e_of_t_by_e( PTTriangle T, PTEdge E )
{
   for( int i = 0; i < 3; i++ )
   {
      if( T->TE[i] == E )  return T->TE[(i+1)%3];
   }
   return 0; // otherwise in VC5: not all control paths return a value
}

PTTriangle Next_CCW_t_of_t_by_v( PTTriangle T, PTVertex V )
{
   return Opposite_t_of_t_by_e( T, Left_e_of_t_by_v( T, V ) );
}

PTTriangle  LastTrgInVE( PTVertex V )  // deve eseere cambiata in LastTrgInVT
{
     if( V->VE[0]->OnConvexHull() )  return Prev_CCW_t_of_e_by_v( LastEdgeInVE( V ) , V );
     else                            return Next_CCW_t_of_e_by_v( LastEdgeInVE( V ) , V );
}

//
// If V is on the convex hull, return the last edge incident in V scanning
// them in cloackwise order starting from the infinite face,
// otherwise return the first edge in cloackwise order starting from
// VE[0], i.e., Prev_CCW_e_of_e_by_v( VE[0] ) .
//
PTEdge  LastEdgeInVE( PTVertex V )
{
     if( V->VE[0]->OnConvexHull() )  return Other_e_in_VE( FirstEdgeInVE( V ), V );
     else                            return Prev_CCW_e_of_e_by_v( V->VE[0], V );
}

//
// If V is on the convex hull, return the first edge incident in V scanning
// them in counterclockwise order starting from the infinite face,
// otherwise return VE[0].
//
PTEdge FirstEdgeInVE( PTVertex V )
{

   #ifdef DEBUG3
     DEBUG3 << "\nTDecCDT::FirstEdgeInVE( " << *V << endl;
   #endif

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

       // search for the vertex of TFirst different from 
       
       int i;
       for( i=0; i<3; i++ )
          if ( v[i] != V && v[i] != VFirst ) break;
	 
	 
       #ifdef ROBUST
          check( (i>=3), "TDecCDT::FindEdgesInVE(), <1> inconsistency detected" );
       #endif
	  
       if ( Geom::Turnxy( V, VFirst, v[i] ) != TURN_LEFT )
          EFirst = V->VE[1];

       return EFirst;
      
   } // end ...if( EFirst->OnConvexHull() )
   else
   return V->VE[0];
}

/* PAOLA JAN. 2000 */ int CheckIsolatedPoint( PTVertex V ) 
{
  if ( (V->VE[0]==NULL) && (V->VE[1]==NULL) )
  {
     cout << "Isolated point " << V->VID << endl;
     return 1;
  }
  return 0;
}

// Check if character c appears in the string letters
//
int IsInString(char c, char * letters)
{
  int i = 0;
  while ( letters[i]!='\0' )
  {  if (c==letters[i]) return 1; else i++; }
  return 0;
}

// Show the prompt and ask the user for one character that must belong
// to the given string allowedletters
//
char AskLetter(char * prompt, char * allowedletters)
{
  char c;
  int i = 0;
  do
  {
      cerr << prompt << " (";
      while ( allowedletters[i]!='\0' )
      {  cerr << allowedletters[i];
         if ( allowedletters[++i]!='\0' ) cerr << "/";
      }
      cerr << ") :" << endl;
      cin  >> c;
  }
  while ( !IsInString(c,allowedletters) );
  return c;
}

// Search the i-th argument string for a character that must belong   
// to the given string of allowedletters
//
char ParseLetter(int argc, char ** argv, int i, char * allowedletters)
{
  char c;
  if (i>=argc)
  {  cerr << "Bad format string: missing mandatory argument" << endl;
     exit(-2);
  }
  if (argv[i][1]!='\0') /* argument string is not a single char */
  {  cerr << "Bad format string: invalid argument " << argv[i] << endl;
     exit(-2);
  }
  c = argv[i][0]; /* first (and only)  character of argument string */  
  if ( !IsInString(c,allowedletters) )
  {  cerr << "Bad format string: expected one of " << allowedletters << endl;
     exit(-2);
  }
  return c;
}

// Show the prompt and ask the user for an integer which must be
// at least as large as lowerbound
//
int AskIntegerAtLeast(char * prompt, int lowerbound)
{
  int v;
  do
  {
      cerr << prompt << "(integer>=" << lowerbound << ")";
      cin  >> v;
  }
  while (v<lowerbound);
  return v;
}

// Search the i-th argument string for an integer which must be
// at least as large as lowerbound
//
int ParseIntegerAtLeast(int argc, char ** argv, int i, int lowerbound)
{
  int v;
  if (i>=argc)
  {  cerr << "Bad format string: missing mandatory argument" << endl;
     exit(-2);
  }
  if ( sscanf(argv[i], "%d", &v) != 1 )
  {  cerr << "Bad format string: invalid argument " << argv[i] << endl;
     exit(-2);
  }
  if (v<lowerbound)
  {  cerr << "Bad format string: invalid argument " << argv[i] << endl;
     exit(-2);
  }
  return v;
}

// Show the prompt and ask the user for a positive float
//
float AskPositiveFloat(char * prompt)
{
  float v;
  int i = 0;
  do
  {
      cerr << prompt << "(positive float)";
      cin  >> v;
  }
  while ( !(v>0.0) );
  return v;
}

// Search the i-th argument string for a positive float
//
float ParsePositiveFloat(int argc, char ** argv, int i)
{
  float v;
  if (i>=argc)
  {  cerr << "Bad format string: missing mandatory argument" << endl;
     exit(-2);
  }
  if ( sscanf(argv[i], "%f", &v) != 1 )
  {  cerr << "Bad format string: invalid argument " << argv[i] << endl;
     exit(-2);
  }
  if ( v<0.0 )
  {  cerr << "Bad format string: invalid argument " << argv[i] << endl;
     exit(-2);
  }
  return v;
}
