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

// ----------------------------------------------------------------------
//
//   file   : utils.h
//   author : Christian Melchiorre
//
//   Some routins for utility purposes (for instance, swap two variables,
//   etc...), used in the Delaunay triangulation programs.
//


#ifndef _UTILS_H
#define _UTILS_H

#include "defs.h"
#include "ttriang.h" // For PTPoint

// ---------------------------------------------------------------------------
//
// Return minimum / maximum between two values, versions for int/double
//

inline int Min( int x, int y )
{ return( (x<y) ? x : y ); };
inline double Min( double x, double y ) { return( (x<y) ? x : y ); };

inline int Max( int x, int y ) { return( (x>y) ? x : y ); };
inline double Max( double x, double y ) { return( (x>y) ? x : y ); };

// ---------------------------------------------------------------------------
//
// Swap the values of two variables. Only for PTPoint (the only used one)
//

inline void Swap( PTPoint& v0, PTPoint& v1 )
{
   PTPoint vTmp = v0;
   v0 = v1;
   v1 = vTmp;
}

extern boolean Is_EV_Defined( PTEdge e );
extern boolean Is_ET_Defined( PTEdge e );
extern boolean Is_v_in_e( PTVertex V, PTEdge E);
extern boolean Is_v_in_t( PTVertex V, PTTriangle T);

extern PTTriangle  Opposite_t_of_t_by_e( PTTriangle T, PTEdge E );
extern PTEdge      Opposite_e_of_t_by_v( PTTriangle T, PTVertex V );

extern PTVertex    Other_v_in_e ( PTEdge e, PTVertex v );
extern PTEdge      Other_e_in_VE ( PTEdge E, PTVertex v );
extern PTEdge      Left_e_of_t_by_v( PTTriangle T, PTVertex V);
extern PTEdge Right_e_of_t_by_v( PTTriangle T, PTVertex V);
extern PTEdge      Prev_CCW_e_of_e_by_v( PTEdge E, PTVertex V );
extern PTEdge      Next_CCW_e_of_e_by_v( PTEdge E, PTVertex V );
extern PTTriangle  Next_CCW_t_of_e_by_v( PTEdge E, PTVertex V );
extern PTTriangle  Prev_CCW_t_of_e_by_v( PTEdge E, PTVertex V );
extern PTEdge      Prev_CCW_e_of_t_by_e( PTTriangle T, PTEdge E );
extern PTEdge      Next_CCW_e_of_t_by_e( PTTriangle T, PTEdge E );
extern PTTriangle  Next_CCW_t_of_t_by_v( PTTriangle T, PTVertex V );
extern PTTriangle  LastTrgInVE( PTVertex V );
extern PTEdge      LastEdgeInVE( PTVertex V );
extern PTEdge      FirstEdgeInVE( PTVertex V );

/* PAOLA GENNAIO 2000 */ extern int CheckIsolatedPoint( PTVertex V );

//
// Functions for supporting interactive mode of the programs
//

// Show the prompt and ask the user for one character that must belong
// to the given string of allowedletters
extern char AskLetter(char * prompt, char * allowedletters);

// Search the i-th argument string for a character that must belong   
// to the given string of allowedletters
extern char ParseLetter(int argc, char ** argv, int i, char * allowedletters);

// Show the prompt and ask the user for an integer which must be
// at least as large as lowerbound
extern int AskIntegerAtLeast(char * prompt, int lowerbound);

// Search the i-th argument string for an integer which must be
// at least as large as lowerbound
extern int ParseIntegerAtLeast(int argc, char ** argv, int i, int lowerbound);

// Show the prompt and ask the user for a positive float
extern float AskPositiveFloat(char * prompt);

// Search the i-th argument string for a positive float
extern float ParsePositiveFloat(int argc, char ** argv, int i);

#endif // _UTILS_H
