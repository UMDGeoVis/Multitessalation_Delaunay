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
//   file   : defs.h
//   author : Christian Melchiorre
//
//   General definitions
//

#ifndef _DEFS_H
#define _DEFS_H


//
// Macros for conditioned compilation
//
// OUTPUT : show for instance the number of entities processed so far
//
// array_Constraints, mark_CONSTRAINED : decimate CDT by using the array
//             "Constraints" or by marking edges as CONSTRAINED.
//
// CC_GCC, CC_SILICON, CC_VISUAL5 : one of them must be defined,
//             according to the compiler.
//
// DEBUG     : show several messages on cout during computation
//             (script show redirects such output on file "output")
// MT_TRACER : enable generation of an MT
// MT_DEBUG  : as DEBUG, but just show which functions of the MT building
//             interface are called.
// ROBUST    : perform several consistency checks for problems that should
//             not arise, but ... one can never know
// PS_OUTPUT : output triangulation contains just x and y, no z,
//            of its vertices


// #define array_Constraints  1
#define mark_CONSTRAINED  2

// #define PAOLO 1

//#define CC_GCC  1
// #define CC_SILICON  2
 #define CC_VISUAL5  3



#define ROBUST  2


//typedef int bool;
#define boolean int

#define CDT_SPIRITO 0

const boolean TRUE = 1;
const boolean FALSE = 0;


//
// These constants are used in the point location process
//

const int PL_UNDEFINED = -1;    // before doing point location
const int PL_EXTERNAL  = 0;     // point outside convex hull
const int PL_TRIANGLE  = 1;     // point inside a triangle
const int PL_EDGE      = 2;     // point inside an edge
const int PL_VERTEX    = 3;     // point coincident with a vertex

//
// These constants define the method for re-calculating the error
// associated with a vertex (exact or approximated).
// They are used in class TDecErrDelaunay
//


const int RECALC_APPROX = 1;
const int RECALC_EXACT  = 2;

#endif // _DEFS_H
