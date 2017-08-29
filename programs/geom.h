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
//   file   : geom.h
//   author : Christian Melchiorre
//
//   Definition of a class collecting all utility static functions 
//   for geometric calculations (for instance, test if three points 
//   define a left turn...).
//
//   Collecting all such functions in a class makes the code clearer.
//    For instance, when we read:
//
//       Geom::CalcDet( ... );
//
//   we can easily find out the module in which the function has been 
//   defined.
//


#ifndef _GEOM_H
#define _GEOM_H

#include "defs.h"
#include "ttriang.h"

//
// Frequeltly used macros
//

#define MAX(x,y)  ((x) < (y) ? (y) : (x))
#define MIN(x,y)  ((x) < (y) ? (x) : (y))

//
// Constants for Geom::Turnxy(...)
//

const int TURN_LEFT  = +1;
const int TURN_RIGHT = -1;
const int ALIGNED    = 0;


// ----------------------------------------------------------------------
//
//  class Geom
//


class Geom
{
   private:
     
      static double EQ_TOLL; // CHANGE IT!!!

      // The default constructor is private. This makes it impossible
      // to create instances of this class.
      Geom() { error("private constructor Geom::Geom() called" ); };
      
   public:

      // Functions to change / check the value of constant EQ_TOLL.
      static void SetTolleranceValue( double newToll );	
      static double GetTolleranceValue();

      // Compute sign of a double and a float number.
      static int Sign( double val );	
      static int Sign( float val );
      
      // Compute absolute value of a double and a float number.
      static double Abs( double val ) ;  
      static float Abs( float val );
	
      // Test if two double numbers are equal/larger/smaller 
      // up to the tolerance contained in constant EQ_TOLL:
      static boolean EqDouble( double val0, double val1 );
      static boolean GtDouble( double val0, double val1 );
      static boolean LtDouble( double val0, double val1 );
 	

      // Compute determinant of a 3x3 matrix.    
      static double Det3x3( double a00, double a01, double a02,
                            double a10, double a11, double a12,
                            double a20, double a21, double a22 );

      // Optimized version for the frequent case in which the third
      // column contains all 1's.
      static double Det3x3( double a00, double a01,
                            double a10, double a11,
			    double a20, double a21 );

      // Version with vertices as parameters.
      static double Det3x3( PTVertex v1,  PTVertex v2,  PTVertex v3 );

      // Test if three points, in the given order, define a 
      // right turn (return +1), a left turn (return -1), or are
      // aligned (return 0).
      // Consider points in the plane (z coordinate is ignored).
      static int Turnxy( PTPoint P0, PTPoint P1, PTPoint P2 );
      
      // Return true iff the three points are aligned, i.e., the
      // determinant is 0 up to the tolerance contained in EQ_TOLL.
      // Consider points in the plane (z coordinate is ignored).
      static boolean Alignedxy( PTPoint P0, PTPoint P1, PTPoint P2 );

      // Return distance between two points.
      // Consider points in the plane (z coordinate is ignored).
      static double Distancexy( PTPoint p0, PTPoint p1 );
      static double Distancexy(double x0, double y0, double x1, double y1 );
      
      
      // Test, by using turns, if point p falls inside or on the boundary
      // of triangle with vertices p0, p1, p2. Vertices p0, p1, p2 must
      // be given in counterclockwise order.

      static boolean InTrianglexy( PTPoint p0, PTPoint p1, PTPoint p2, PTPoint p )
      {
         #ifdef ROBUST 
	     check( (Geom::Turnxy( p0, p1, p2 ) != TURN_LEFT),
	       "Geom::InTrianglexy(), points are not in ccw order" );
	 #endif
	 
         return(
	         Geom::Turnxy( p0, p1, p ) != TURN_RIGHT
		 &&
		 Geom::Turnxy( p1, p2, p ) != TURN_RIGHT
		 &&
		 Geom::Turnxy( p2, p0, p ) != TURN_RIGHT
	       );
      }
       

      // Return, in the last three parameters (passed as references),
      // the coordinates of the center and the radius of the circle
      // psssing through the points given in the first three parameters.
      // Consider only xy coordinates, z coordinate is ignored.
      static void CalcCirclexy( PTPoint p0, PTPoint p1, PTPoint p2,
                                double &CX, double &CY, double &CR );
    

      // Return TRUE if the three distict vertices are aligned.
      static int Collinear( PTVertex v0, PTVertex v1, PTVertex v2);

      // Given a point and an edge/triangle, the following functions
      // compute the z coordinate of the vertical projection ot the point
      // on the given edge/triangle.      
      static double Edgez( PTEdge E, PTPoint P );
      static double Trianglez( PTTriangle T, PTPoint P );
      static double Trianglez( PTPoint, PTPoint, PTPoint, PTPoint );
      
      // Compute the error of a point w.r.t. the edge/triangle in which 
      // the vertical projection ot the point falls, i.e., the 
      // difference between the z coordinate of the point and the 
      // z coordinate of its projection.
      static double CalcError( PTEdge E, PTPoint P );
      static double CalcError( PTTriangle T, PTPoint P );
      static double Dist(PTTriangle TCurr, PTPoint P);

};

#endif // _GEOM_H
