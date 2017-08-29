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
//   file   : geom.cpp
//   author : Christian Melchiorre
//
//   Definition of a class collecting all utility static functions 
//   for geometric calculations (for instance, test if three points 
//   define a left turn...).
//
//   Collecting all such functions in a class makes the code clearer.
//   For instance, when we read:
//
//       Geom::CalcDet( ... );
//
//   we can easily find out the module in which the function has been 
//   defined.
//

#include "defs.h"
#include "geom.h"
#include "ttriang.h" 

#include <math.h> // ...for sqrt()

double Geom::EQ_TOLL = 1E-12;

// -----------------------------------------------------------------
//  
// void Geom::SetTolleranceValue( double newToll )
// double Geom::GetTolleranceValue()
//
// Functions to change/check the value of constant EQ_TOLL.
//

void Geom::SetTolleranceValue( double newToll )
        { EQ_TOLL = newToll; }
	
double Geom::GetTolleranceValue()
        { return(EQ_TOLL); }
	

// -----------------------------------------------------------------
//
// int Sign( double/float val )
//
// Compute the sign of a double and a float number.
//

int Geom::Sign( double val )
        { return( (val<0 ? -1 : (val>0 ? 1 : 0)) ); }
	
int Geom::Sign( float val )
        { return( (val<0 ? -1 : (val>0 ? 1 : 0)) ); }
      
      
// -----------------------------------------------------------------   
//
// double/float Abs ( double/float val )
//
// Compute the absolute valuen of a double and a float number.
//

double Geom::Abs( double val ) 
        { return( (val<0 ? -val : (val == -0.0 ? 0.0 : val)) ); }
  
float Geom::Abs( float val ) 
        { return( (val<0 ? -val : (val == -0.0 ? 0.0 : val)) ); }
	
	
// -----------------------------------------------------------------
//
// boolean Eq/Gt/LtDouble( double, double )
//
// Test if two double numbers are equal/larger/smaller 
// up to the tolerance contained in constant EQ_TOLL.
//

boolean Geom::EqDouble( double val0, double val1 )
        { return( Abs(val0-val1) <= EQ_TOLL ); }

boolean Geom::GtDouble( double val0, double val1 )
        { return( (!EqDouble(val0, val1)) && val0 > val1 ); };

boolean Geom::LtDouble( double val0, double val1 )
        { return( (!EqDouble(val0, val1)) && val0 < val1 ); };


// -----------------------------------------------------------------    
//
// double Geom3x3(...)
//
// Compute determinant of a 3x3 matrix.    
//

double Geom::Det3x3( double a00, double a01, double a02,
                     double a10, double a11, double a12,
                     double a20, double a21, double a22 )
{
   
   double det = (
                    (a00 * a11 * a22) +
	            (a01 * a12 * a20) +
	            (a02 * a10 * a21) -
	            (a20 * a11 * a02) -
	            (a21 * a12 * a00) -
	            (a22 * a10 * a01)
                );	

   if ( Abs(det) > EQ_TOLL ) return(det); else return(0.0);
}


//
// Optimized version for the frequent case in which the third
// column contains all 1's.
//

double Geom::Det3x3( double a00, double a01,
                     double a10, double a11,
                     double a20, double a21 )
{
   double det = (
                 (a00 * a11) +
		 (a01 * a20) +
		 (a10 * a21) -
		 (a20 * a11) -
		 (a21 * a00) -
		 (a10 * a01)
	        );
	      
   if ( Abs(det) > EQ_TOLL ) return(det); else return(0.0);
 		    
}
	
  
//
// Version with vertices as parameters.
//
    
double Geom::Det3x3( PTVertex v1,  PTVertex v2,  PTVertex v3 )
{

   double a00 = v1->x;
   double a01 = v1->y;
   double a02 = v1->z;
   double a10 = v2->x;
   double a11 = v2->y;
   double a12 = v2->z;
   double a20 = v3->x;
   double a21 = v3->y;
   double a22 = v3->z;
   
   double det = (
                    (a00 * a11 * a22) +
	            (a01 * a12 * a20) +
	            (a02 * a10 * a21) -
	            (a20 * a11 * a02) -
	            (a21 * a12 * a00) -
	            (a22 * a10 * a01)
                );	

   if ( Abs(det) > EQ_TOLL ) return(det); else return(0.0);
}


// -----------------------------------------------------------------
//
// int Geom::Turnxy( PTPoint, PTPoint, PTPoint )
//
// Test if three points, in the given order, define a right turn
// (return +1), a left turn (return -1), or are aligned (return 0).
// Consider points in the plane (z coordinate is ignored).
//

int Geom::Turnxy( PTPoint P0, PTPoint P1, PTPoint P2 )
{ 
    return( 
           Sign( 
	        Det3x3( P0->x, P0->y,
	                P1->x, P1->y,
		        P2->x, P2->y )
	       )
	  );
}
  
      
// -----------------------------------------------------------------
//
// boolean Geom::Alignedxy( PTPoint, PTPoint, PTPoint )
//
// Return true iff the three points are aligned, i.e., the
// determinant is 0 up to the tolerance contained in EQ_TOLL.
// Consider points in the plane (z coordinate is ignored).
//

boolean Geom::Alignedxy( PTPoint P0, PTPoint P1, PTPoint P2 )
{ 
   return( Turnxy( P0, P1, P2 ) == ALIGNED );  // ALIGNED == 0
}


// -----------------------------------------------------------------
//
// double Geom::Distancexy( PTPoint, PTPoint )
//
// Return distance between two points.
// Consider points in the plane (z coordinate is ignored).
//
   
double Geom::Distancexy(double x0, double y0, double x1, double y1 )
{
   double dx = Geom::Abs( x0 - x1 );
   double dy = Geom::Abs( y0 - y1 );

   return( sqrt( (dx*dx)+(dy*dy) ) );
}

      
double Geom::Distancexy( PTPoint p0, PTPoint p1 )
{
   return( Geom::Distancexy( p0->x, p0->y, p1->x, p1->y ) );	 
}
     
      
// -----------------------------------------------------------------
//
// void Geom::CalcCirclexy( PTPoint, PTPoint, PTPoint, 
//                          double&, double&, double& )
//
// Return, in the last three parameters (passed as references),
// the coordinates of the center and the radius of the circle passing
// through the points given in the first three parameters.
// Consider only xy coordinates, z coordinate is ignored.
//
   
void Geom::CalcCirclexy( PTPoint p0, PTPoint p1, PTPoint p2,
                         double &CX, double &CY, double &CR )
{

   double A = p1->x - p0->x;
   double B = p1->y - p0->y;
   double C = p2->x - p0->x;
   double D = p2->y - p0->y;
   
   double E = ( A * ( p0->x + p1->x ) ) + ( B * ( p0->y + p1->y ) );
   double F = ( C * ( p0->x + p2->x ) ) + ( D * ( p0->y + p2->y ) );
   
   double G = 2.0 * ( ( A * ( p2->y - p1->y ) ) - ( B * ( p2->x - p1->x ) ) );

   if (EqDouble(G, 0.0)) cerr << "\n Points " << *p0 << "  " << *p1 << "  " << *p2 << endl;
   check( (EqDouble(G, 0.0)), "Geom::CalcCirclexy(), points are collinear");
   
   CX = ( (D*E) - (B*F) ) / G;
   CY = ( (A*F) - (C*E) ) / G;
   
   double dx = ( p0->x - CX );
   double dy = ( p0->y - CY );
   
   CR = sqrt( (dx * dx) + (dy * dy) );
   
/* 
   ----------------------------------------------------------------------
   
    Found on internet...

    Subject 1.04: How do I generate a circle through three points?

    Let the three given points be a, b, c.  Use _0 and _1 to represent
    x and y coordinates. The coordinates of the center p=(p_0,p_1) of
    the circle determined by a, b, and c are:

        A = b_0 - a_0;
        B = b_1 - a_1;
        C = c_0 - a_0;
        D = c_1 - a_1;
    
        E = A*(a_0 + b_0) + B*(a_1 + b_1);
        F = C*(a_0 + c_0) + D*(a_1 + c_1);
    
        G = 2.0*(A*(c_1 - b_1)-B*(c_0 - b_0));
    
        p_0 = (D*E - B*F) / G;
        p_1 = (A*F - C*E) / G;
 
    If G is zero then the three points are collinear and no finite-radius
    circle through them exists.  Otherwise, the radius of the circle is:

            r^2 = (a_0 - p_0)^2 + (a_1 - p_1)^2

   ----------------------------------------------------------------------
*/
   
}


//
// Return TRUE if the three distict vertices are aligned.
//

int Geom::Collinear( PTVertex v0, PTVertex v1, PTVertex v2)
{   
   double g =  (v1->x - v0->x) * ( v2->y - v1->y );
   g  -=  (v1->y - v0->y) * ( v2->x - v1->x ) ;
   return EqDouble(g, 0.0);
}


// -----------------------------------------------------------------
//
//  double Geom::Edgez( PTEdge E, PTPoint P )
//
//  Compute the heigth of the vertical projection of P over edge E.
//  Computation is based on the parametric formula:
//  P(t) = P0 + (P0-P1)t
//

double Geom::Edgez( PTEdge E, PTPoint P )
{

   double x0 = E->EV[0]->x;
   double y0 = E->EV[0]->y;
   double z0 = E->EV[0]->z;
   
   double x1 = E->EV[1]->x;
   double y1 = E->EV[1]->y;
   double z1 = E->EV[1]->z;

   boolean EqX = Geom::EqDouble( x0, x1 );
   boolean EqY = Geom::EqDouble( y0, y1 );

   
   #ifdef ROBUST
       
       //
       // chack that point really falls inside edge
       //
       
       double minx = ( x0 < x1 ? x0 : x1 ); 
       double maxx = ( x0 >= x1 ? x0 : x1 );
       double miny = ( y0 < y1 ? y0 : y1 );
       double maxy = ( y0 >= y1 ? y0 : y1 );
       check( ( P->x < minx || P->y < miny || P->x > maxx || P->y > maxy ||
       		!Geom::Alignedxy( P, E->EV[0], E->EV[1] ) ),
		"Geom::Edgez(), the point does not lie inside the edge" );
		
   #endif // ROBUST
   
   
   if ( EqX && EqY )
   {
       error( "Geom::Edgez(), edge's extreme points coincident" );   
       return(0.0); // ...to avoid warning from compiler
   }   
   else if ( !EqX )
   {
       return( z0 + ( ((z1 - z0) * (P->x - x0)) / (x1 - x0) ) );
   }
   else
   {
       return( z0 + ( ((z1 - z0) * (P->y - y0)) / (y1 - y0) ) );
   }
      
}




// -----------------------------------------------------------------
//
//  double Geom::Trianglez( PTTriangle T, PTPoint P )
//
//  Compute the heigth of the vertical projection of P over triangle T.
//  Computation is based on setting the determinant of a 4x4 matrix to zero.
//  No check is made to make sure that P really falls inside T, as for
//  Edgez().
//

double Geom::Trianglez( PTTriangle T, PTPoint P )
{
   PTVertex V0, V1, V2;
   T->GetTV( V0, V1, V2 );
   return( Geom::Trianglez( V0, V1, V2, P ) );
}

// -----------------------------------------------------------------
//
//  double Geom::Trianglez( PTPoint v0, PTPoint v1, PTPoint v2, PTPoint P )
//
//  Compute the heigth of the vertical projection of P over the triangle 
//  formed by vertices v0, v1, v2.
//  Computation is based on setting the determinant of a 4x4 matrix to zero.
//  No check is made to make sure that P really falls inside T, as for
//  Edgez().
//

double Geom::Trianglez( PTPoint V0, PTPoint V1, PTPoint V2, PTPoint P )
{

   double x0 = V0->x;
   double y0 = V0->y;
   double z0 = V0->z;
   
   double x1 = V1->x;
   double y1 = V1->y;
   double z1 = V1->z;
   
   double x2 = V2->x;
   double y2 = V2->y;
   double z2 = V2->z;
   
   #ifdef ROBUST
   
       //
       // check that point really falls inside triangle
       //
       
       check( ( Geom::Turnxy( V0, V1, P ) == TURN_RIGHT ||
                Geom::Turnxy( V1, V2, P ) == TURN_RIGHT ||
                Geom::Turnxy( V2, V0, P ) == TURN_RIGHT ),
            "Geom::Trianglez(), the point does not lie inside the triangle" );    
   
   #endif // ROBUST

   double a1 = y0 * (z1 - z2);
          a1 -= z0 * (y1 - y2);
          a1 += (y1 * z2 - y2 * z1);
   double a2 = x0 * (z1 - z2);
          a2 -= z0 * (x1 - x2);
          a2 += (x1 * z2 - x2 * z1);
   double a3 = x0 * (y1 * z2 - y2 * z1);
          a3 -= y0 * (x1 * z2 - x2 * z1);
          a3 += z0 * (x1 * y2 - x2 * y1);
   double a = P->x * a1 - P->y * a2 - a3;

// PAOLA: it must be done step by step, as above. If we do it all together
// sometimes it returns NaN, and I don't understand why!
//   double a = P->x * (y0 * (z1 - z2) - z0 * (y1 - y2) + (y1 * z2 - y2 * z1)) - 
//              P->y * (x0 * (z1 - z2) - z0 * (x1 - x2) + (x1 * z2 - x2 * z1)) - 
//              (x0 * (y1 * z2 - y2 * z1) - y0 * (x1 * z2 - x2 * z1) + z0 * (x1 * y2 - x2 * y1));

   double b = x0 * (y1 - y2);
          b -= y0 * (x1 - x2);
          b += (x1 * y2 - x2 * y1);

// PAOLA
//cerr << "Geom::Trianglez a1=" << a1 << " a2=" << a2 << " a3=" << a3
//     << " a=" << a << " b=" << b << " ris=" << ( -(a / b) ) << endl;
// PAOLA

//   #ifdef ROBUST
      // ...I don't know what it means b=0...
      check( (b == 0), "Geom::Trianglez(), this should not happen ( b == 0 )" );
//   #endif      
      	   
   return( -( a / b ) );
   
}


// -----------------------------------------------------------------
//
//  double Geom::CalcError( PTTriangle/PTEdge, PTPoint )
//
//  Compute the error of a point w.r.t. the triangle/edge that 
//  contains the vertical projection of the point inside it.
//  Such error is the difference between the z coordinate of the 
//  point and the z coordinate of its projection on the triangle/edge.
//

double Geom::CalcError( PTEdge E, PTPoint P )
{
   #ifdef DEBUG
        DEBUG << endl << endl << "Errore del punto ( " << *P << " ) nello spigolo E"
              <<  E->EID << endl
              <<  (*(E->EV[0])) << " / " << (*(E->EV[1]))
              << endl << " ERR: " << Geom::Abs( P->z - Geom::Edgez( E, P ) )
              << endl;
   #endif

   return( Geom::Abs( P->z - Geom::Edgez( E, P ) ) );
}

double Geom::CalcError( PTTriangle T, PTPoint P )
{
   #ifdef DEBUG
   
        PTVertex v0, v1, v2;
	T->GetTV( v0, v1, v2 );
        DEBUG << endl << endl <<  "Error of point ( " << *P << " ) in triangle T"
              <<  T->TID << endl
	      << (*(PTPoint)v0) << " / " << (*(PTPoint)v1) << " / " << (*(PTPoint)v2)
              << endl << " ERR: " << Geom::Abs( P->z - Geom::Trianglez( T, P ) )
              << endl << endl;
   #endif

   return( Geom::Abs( P->z - Geom::Trianglez( T, P ) ) );
}

double Geom::Dist(PTTriangle TCurr, PTPoint P)
{
    int i;
	double V1[3],V2[3],W[3];
	PTVertex V[3];

	TCurr->GetTV( V[0], V[1], V[2] );

	V1[0]=V[1]->x-V[0]->x;
	V1[1]=V[1]->y-V[0]->y;
	V1[2]=V[1]->z-V[0]->z;

	V2[0]=V[2]->x-V[0]->x;
	V2[1]=V[2]->y-V[0]->y;
	V2[2]=V[2]->z-V[0]->z;

	W[0]=P->x-V[0]->x;
	W[1]=P->y-V[0]->y;
	W[2]=P->z-V[0]->z;

	double A=V1[1]*V2[2]-V1[2]*V2[1];
	double B=V1[2]*V2[0]-V1[0]*V2[2];
	double C=V1[0]*V2[1]-V1[1]*V2[0];
	double norm=sqrt(A*A+B*B+C*C);		

	return (Geom::Abs(A*W[0]+B*W[1]+C*W[2])/norm);
}
