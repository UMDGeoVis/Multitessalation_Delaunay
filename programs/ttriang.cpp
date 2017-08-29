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
//   file   : ttriang.cpp
//   author : Christian Melchiorre
//
//   Definitions of basic objects used in a triangulation.
//

#define TREDIM

#include "defs.h"
#include "error.h"
#include "geom.h"
#include "ttriang.h"

int TVertex::NextVID = 0;
int TEdge::NextEID = 0;
int TTriangle::NextTID = 0;

// -----------------------------------------------------------------------------
//
//  methods of class TPoint
//



// -----------------------------------------------------------------------------
//
//  boolean Equals/Equalsxy( PTPoint )
//
//  Compare the coordinates of two points and check if they are equal.
//  Also version tha considers only x and y coordinates.
//

boolean TPoint::Equals( PTPoint p ) 
{ 
   return( Geom::EqDouble( x, p->x ) &&
           Geom::EqDouble( y, p->y ) &&
           Geom::EqDouble( z, p->z )
         );
};  

boolean TPoint::Equalsxy( PTPoint p ) 
{ 
   return( Geom::EqDouble( x, p->x ) &&
           Geom::EqDouble( y, p->y )
         );
};  
  
// -----------------------------------------------------------------------------
//
//  Input/output operators for class TPoint.
//

ostream& operator<< ( ostream& os, RTPoint p )
{
  #ifdef TREDIM
  os << p.x << " " << p.y << " " << p.z;
  #else
  os << p.x << " " << p.y;
  #endif
  return(os);
}

istream& operator>> ( istream& is, RTPoint p )
{
  #ifdef TREDIM
  is >> p.x >> p.y >> p.z;
  #else
  is >> p.x >> p.y; p.z=0.0;
  #endif
  return(is);
}


// -----------------------------------------------------------------------------
//
//  Methods of class TVertex.
//

// -----------------------------------------------------------------------------
//
//  Build the list EdgeList of edges adjacent to this vertex.
//

void TVertex::GetVE()
{	PTEdge FirstEdge = this->VE[0];
	PTEdge CurrentEdge = this->VE[0];	
	PTTriangle CurrentTriangle = NULL;
	
        #ifdef ROBUST
          	check( (FirstEdge == NULL), "TVertex::GetVE: inconsistency!" );
        #endif

//	#ifdef DEBUG
//	SKIPPO
//	DEBUG << "Edges del vertice "<< this->VID<<": "<<*FirstEdge;
//	#endif

	EdgeList.ClearList();
	if (FirstEdge->OnConvexHull()) EdgeList.AddHead(FirstEdge); 
    
	do
	{	CurrentTriangle = (CurrentTriangle==CurrentEdge->ET[0]) ? CurrentEdge->ET[1] : CurrentEdge->ET[0];
		//cerco il prox edge e lo aggiungo in lista
		for (int i=0; i<3; i++)
			if (CurrentTriangle->TE[i] != CurrentEdge)
				if ((CurrentTriangle->TE[i]->EV[0] == this) ||
				    (CurrentTriangle->TE[i]->EV[1] == this))
					{	CurrentEdge = CurrentTriangle->TE[i];
						break;
					}
//        	#ifdef DEBUG
//                DEBUG << ", "<<*CurrentEdge;
//        	#endif
				
		EdgeList.AddHead(CurrentEdge);
		

	}
	while ((CurrentEdge != FirstEdge) && (!(CurrentEdge->OnConvexHull())) );
//	#ifdef DEBUG
//	DEBUG <<endl;
//	#endif

}


// -----------------------------------------------------------------------------
//
//  Constructors
//


TVertex::TVertex( double xi=0, double yi=0, double zi=0 )
   : TPoint( xi, yi, zi ), EdgeList ()
{
	nIncConstr = 0;
   VID = NextVID++;
   VE[0] = VE[1] = NULL;
#ifdef MT_TRACER
   MTIdx = 0;
#endif
}

TVertex::TVertex( PTPoint p )
   : TPoint( p->x, p->y, p->z ), EdgeList ()
{
	nIncConstr = 0;
   VID = NextVID++;
   VE[0] = VE[1] = NULL;
#ifdef MT_TRACER
   MTIdx = 0;
#endif
   Error = p->Error;

   #ifdef DEBUG
     DEBUG << "constructed TVertex( V" << VID << ": " 
		<< x << ", " << y << ", " << z << " )" << endl;
   #endif 

}
   

TVertex::TVertex( RTPoint p )
   : TPoint( p ), EdgeList ()
{
	nIncConstr = 0;
   VID = NextVID++;
   VE[0] = VE[1] = NULL;
#ifdef MT_TRACER
   MTIdx = 0;
#endif
}



// -----------------------------------------------------------------------------
//
//   Destructor of class TVertex
//
//   Make sure that possible pointers to this in edges adjacent to it
//   are set to NULL. Empty the list of adjacent edges.
//


TVertex::~TVertex()
{
   for ( int i=0; i<2; i++ )
   {
      if ( VE[i] != NULL )
      {
         if ( VE[i]->EV[0] == this ) VE[i]->EV[0] = NULL;
	 else
	 if ( VE[i]->EV[1] == this ) VE[i]->EV[1] = NULL;
      }
   }
   EdgeList.ClearList();
}
 

// -----------------------------------------------------------------------------
//
//  Output operator for class TVertex
//

ostream& operator<< ( ostream& os, RTVertex V )
{
  os << "V" << V.VID << "[ " << (*((PTPoint)&V)) << " ]";
  return( os );
}


// -----------------------------------------------------------------------------
//
//   Methods class TEdge.
//


// -----------------------------------------------------------------------------
//
//   Constructors.
//


TEdge::TEdge( PTVertex v0, PTVertex v1 ) : PointList()
{

   EID = NextEID++;

   EV[0] = v0;
   EV[1] = v1;
   
   ET[0] = NULL;
   ET[1] = NULL;
   ToVertex=NULL;
   Error=0;

};


TEdge::TEdge() : PointList()
{

   EID = NextEID++;

   EV[0] = NULL;
   EV[1] = NULL;
   
   ET[0] = NULL;
   ET[1] = NULL;

   ToVertex=NULL;
   Error=0;
   this->Mark( NULLMARK );

};


// copy constructor //
TEdge::TEdge(TEdge &e)
{
   *this = e;
   PointList = TList<PTPoint>(); // 2DO: must be replaced with copy of the list
   // EID = e.EID;
   // ET[0] = e.ET[0];
   // ET[1] = e.ET[1];
   // EV[0] = e.EV[0];
   // EV[1] = e.EV[1];
   // NextEID = e.NextEID;
   // PointList = e.PointList;
   // mark = e.MarkValue();
}

// -----------------------------------------------------------------------------
//
//  void TEdge::AddPoint( PTPoint )
//
//  Insert the given point into the PointList of this edge.
//  Compute the real error of the point in the current triangulation:
//  if such error is larger than the error of the point at the head
//  of the PointList, then insert this point at the beginning;
//  else insert it in second position.
//

void TEdge::AddPoint( PTPoint PointToAdd )
{
    double NewError = Geom::CalcError( this, PointToAdd );
    PointToAdd->Error = NewError;

    double OldPtErr = ( PointList.IsEmpty() ? 0.0 : PointList.GetHead()->Error );

    //
    // If the point being inserted is the one that causes the max error,
    // then insert it at the head of PointList; otherwise, insert it
    // in second position.
    //
    
    if ( NewError >= OldPtErr )
    {    
       PointList.AddHead( PointToAdd );
    }
    else
    {
       TListIterator<PTPoint> PLIter( &PointList );
       PLIter.Restart();
       PointList.AddAfter( PLIter.Current(), PointToAdd );
    }
}



//
// Return TRUE if this edge is matching with filter or
// if they point to the same edge object/edge.
// To "match with" means that all comparable fields of filter that are
// non-zero must be equal, i.e.:
// -  if in filter field EV[?] and ET[?] are not NULL, then they must
//    be equal up the commutativity of EV and ET.
// -  field mark, instead, must be compared only with active mark
//    values in filter.
//    (e.g., if filter is marked only CONSTRAINED and TO_DELETE, then
//    all edges with both such mark values "pass" the filter).
//
boolean TEdge::Match( PTEdge filter ) {
      

   if( this == filter ) return TRUE;  // they are the same object

   //
   // filter fields EV
   //

   if( filter->EV[0] != NULL || filter->EV[1] != NULL )
   {  
      // filter.EV : either they are both NULL or none of them is 
      // null and they are different
      check( filter->EV[0] == NULL || filter->EV[1] == NULL || filter->EV[0] == filter->EV[1],
             "TEdge::Match,<1>  filter->EV  error");

      // degenerate edges are not accepted (by now, but should reason about it)
      check( this->EV[0] ==  this->EV[1], "TEdge::Match, calling edge has degenerate EV");
      
      // if we arrive here, filter->EV are different and also this->EV,
      // thus it is enough to check that filter->EV are equal to at 
      // least one of the elements of this->EV fields
      if( filter->EV[0] != this->EV[0] && filter->EV[1] != this->EV[0] ) return FALSE;
      // if( filter->EV[1] != this->EV[0] && filter->EV[1] != this->EV[1] ) return FALSE;
      if( filter->EV[0] != this->EV[1] && filter->EV[1] != this->EV[1] ) return FALSE;
   }

   //
   // filter fields ET: differently from EV, we accept that one field is
   // NULL and the other one is non-NULL, because it may be that the edge
   // is on the convex hull or is on the boundary of a region of influence
   // and its ET has not yet been completely set.
   // Therefore, accept more edges w.r.t. the filter of EV since we are
   // not able to discriminate between
   // - "comparing only one of the two ET fields" and
   // - "comparing both fields, searching for a boundary edge"
   // But, in summary, it is not interesting to "compare only one of the
   // two ET fields", therefore the rationale for ET changes:
   // if filter.ET are both NULL, it means that we don't whant to 
   // compare such fields; otherwise we compare both.
   // In this way we are able to search for boundary edges (which is a
   // useless thing if we search for those on the convex hull since they
   // should already be found with VE*, but one never knows).
   //

   if( filter->ET[0] == NULL && filter->ET[1] == NULL ) ;
   else // filter.ET :  at least one is not NULL
   {  
      // filter.ET : either they are both NULL or they are different
      check( filter->ET[0] != NULL && filter->ET[1] != NULL && filter->ET[0] == filter->ET[1],
             "TEdge::Match,<2>  filter->EV  error");

      // filter.ET : equal values are not accepted (by now, but should
      // reason about it), but accept that ET has not yet been set
      check( this->ET[0] != NULL  &&  this->ET[1] != NULL  &&  this->ET[0] == this->ET[1],
             "TEdge::Match, calling edhe has degenerate ET");

/* Version that compares both ET fields */

      // if we arrive here, filter->ET are different and also this->ET, thus
      // it is enough to check that they are equal to at least one of
      // the elements of this->EV
      if( filter->ET[0] != this->ET[0] && filter->ET[1] != this->ET[0] ) return FALSE;
      if( filter->ET[1] != this->ET[0] && filter->ET[1] != this->ET[1] ) return FALSE;


/* Version that allows to "compare only one ot the ET fields" 

      // if we arrive here, filter->ET are different and also this->ET, thus
      // if filter->EV are not null
      // it is enough to check that they are equal to at least one of
      // the elements of this->EV
      if( filter->ET[0] != NULL && filter->ET[0] != this->ET[0] && filter->ET[1] != this->ET[0] )
         return FALSE;
      if( filter->ET[1] != NULL && filter->ET[1] != this->ET[0] && filter->ET[1] != this->ET[1] )
         return FALSE;
*/
   }

   //
   // filter mark field
   //

   for( MARKTYPE m = FIRST_MARK;  m <= LAST_MARK;  m *= 2 )
      if( filter->Marked( m )  && ! this->Marked( m ) ) return FALSE;

   // if we arrive here, this has passed all tests of the filter

   return TRUE;
}
 

boolean TEdge::OldMatch( PTEdge filter ) {

   if( this == filter ) return TRUE;  // sono lo stesso oggetto
   if( filter->EV[0] != NULL   &&   filter->EV[0] != this->EV[0] ) return FALSE;
   if( filter->EV[1] != NULL   &&   filter->EV[1] != this->EV[1] ) return FALSE;
   if( filter->ET[0] != NULL   &&   filter->ET[0] != this->ET[0] ) return FALSE;
   if( filter->ET[1] != NULL   &&   filter->ET[1] != this->ET[1] ) return FALSE;
   for( MARKTYPE m = FIRST_MARK;  m <= LAST_MARK;  m *= 2 )
      if( filter->Marked( m )  && ! this->Marked( m ) ) return FALSE;
   return TRUE;
}




// -----------------------------------------------------------------------------
//  
//   Destructor of class TEdge
//
//   Make sure that possible pointers to this in vertices or triangles 
//   adjacent to it are set to NULL. Empty the PointList.
//

TEdge::~TEdge()
{
   for( int i=0; i<2; i++ )
   {
      if ( this->EV[i] != NULL )
      {
         if ( this->EV[i]->VE[0] == this ) this->EV[i]->VE[0] = NULL;
	 else
	 if ( this->EV[i]->VE[1] == this ) this->EV[i]->VE[1] = NULL;
      }
      
      if ( this->ET[i] != NULL )
      {
         for ( int e=0; e<3; e++ )
	   if ( this->ET[i]->TE[e] == this )
	   {
	       this->ET[i]->TE[e] = NULL;
	       break;
	   }
      }
      
   }
   
   PointList.ClearList();
}



// -----------------------------------------------------------------------------
//
//  Ouput operator for class TEdge
//

ostream& operator<< ( ostream& os, RTEdge E )
{
   os << "E" << E.EID << "[ V";
   if( E.EV[0] != NULL ) os << E.EV[0]->VID; else os << " NULL";
   os << ", V";
   if( E.EV[1] != NULL ) os << E.EV[1]->VID; else os << " NULL";
   os << ", T";
   if( E.ET[0] != NULL ) os << E.ET[0]->TID; else os << " NULL";
   os << ", T";
   if( E.ET[1] != NULL ) os << E.ET[1]->TID; else os << " NULL";
   os << ", mark " << E.MarkValue() << " ]";
   return( os );
}




// -----------------------------------------------------------------------------
//
//   Methods class TTriangle.
//


// -----------------------------------------------------------------------------
//
//   Constructor.
//


TTriangle::TTriangle( PTEdge e0, PTEdge e1, PTEdge e2 ) : PointList()
{

   TID = NextTID++;

   TE[0] = e0;
   TE[1] = e1;
   TE[2] = e2;

   if (e0 && e1 && e2) CalcCircle();//PAOLA 22 FEB. 2001
//PAOLA 22 FEB. 2001   CalcCircle();

}


// -----------------------------------------------------------------------------
//
//   void TTriangle::CalcCircle()
//
// Compute the coordinates and the radius of the circum-circle.
//

void TTriangle::CalcCircle()
{
   PTVertex v0, v1, v2;
   this->GetTV( v0, v1, v2 );
   
   // ...piccolo controllo
   check ( (v0==NULL || v1==NULL || v2==NULL),
   	"TTriangle::CalcCircle(), inconsistency detected" );
   
   Geom::CalcCirclexy( v0, v1, v2, InCircleX, InCircleY, InCircleRad );

}



// -----------------------------------------------------------------------------
//
//   boolean TTriangle::InCircle( PTPoint p )
//
//   Return TRUE iff point p is inside circum-circle of this triangle.
//

boolean TTriangle::InCircle( PTPoint p )
{
   double Dist = Geom::Distancexy( p->x, p->y, InCircleX, InCircleY );
   
   return( Geom::LtDouble( Dist, InCircleRad ) );
}



// -----------------------------------------------------------------------------
//
//   void TTriangle::GetTV( PTVertex &v0, PTVertex &v1, PTVertex &v2 )
//
//   Compute Triangle-Vertex relation by using Triangle-Edge and
//   Edge-Vertex relations.
//   The three vertex pointers are returned, in counterclockwise order, 
//   in the three parameters passed by reference.
//
//   Note: If we call:
//
//      T.GetTV( v0, v1, v2 );
//
//    for all i=0..2, T.TE[i] is the edge of endpoints  v(i), v((i+1)%2)
//

void TTriangle::GetTV( PTVertex &v0, PTVertex &v1, PTVertex &v2 )
{

   check( (TE[0] == NULL || TE[1] == NULL || TE[2] == NULL),
          "TTriangle::GetTV(), inconsistent data structure detected" );
  
   //
   // v0 and v1 are the two endpoints of this->TE[0]
   //
   
   v0 = TE[0]->EV[0];
   v1 = TE[0]->EV[1];
   
   //
   // v2 is the vertex of this->TE[1] different from v0 and v1
   //
   
   v2 = ( TE[1]->EV[0] != v0 && TE[1]->EV[0] != v1 ? TE[1]->EV[0] : TE[1]->EV[1] );

   //
   // not sure that triplet v0, v1, v2 is in counterclockwise order,
   // therefore check and possibly swap two vertices to change the order
   //  
   
   if ( Geom::Turnxy( v0, v1, v2 ) != TURN_LEFT ) 
   {
     PTVertex vTmp = v0;
     v0 = v1;
     v1 = vTmp;
   }

}



// -----------------------------------------------------------------------------
//  
//  void TTriangle::GetTT( PTTriangle &t0, PTTriangle &t1, PTTriangle &t2 )
//
//  Compute Triangle-Triangle relation starting from Triangle-Edge and
//  Edge-Triangle relations.
//
//  In the first version, the three triangle pointers are returned,
//  in counterclockwise order, in the three parameters passed by 
//  reference.
//
//  Note: If we call:
//    
//     T.GetTT( t0, t1, t2 )
//
//  then, for all i=0..2, T.TE[i] is the common edge of this and ti
//

 
void TTriangle::GetTT( PTTriangle &t0, PTTriangle &t1, PTTriangle &t2 )
{
   check( (TE[0] == NULL || TE[1] == NULL || TE[2] == NULL),
          "TTriangle::GetTT(), inconsistent data structure detected" );
	  
   t0 = ( TE[0]->ET[0] != this ? TE[0]->ET[0] : TE[0]->ET[1] );
   t1 = ( TE[1]->ET[0] != this ? TE[1]->ET[0] : TE[1]->ET[1] );
   t2 = ( TE[2]->ET[0] != this ? TE[2]->ET[0] : TE[2]->ET[1] );
   
   //
   // since TE[0], TE[1], TE[2] are in counterclockwise order, t0,t1,t2 
   // (computed in this way) are in counterclockwise order, too.
   // Therefore, we can exit the function without further checking.
   //
   
}


//
//  Versione taking a PTEdge as parameter: return only the
//  triangle adjacent to this along such edge.
//

PTTriangle TTriangle::GetTT( PTEdge e )
{
   for( int i=0; i<3; i++ )
   {
     if ( TE[i] == e )
         return(  TE[i]->ET[0] != this ? TE[i]->ET[0] : TE[i]->ET[1] );
   }
     
   error( "TTriangle::GetTT( PTEdge ), edge is not adjacent to this" );

   return(NULL); // ...to avoid warning from the compiler
}



//
//  Return only the adjacent triangle along i-th edge (i = 0..2)
//

PTTriangle TTriangle::GetTT( int i )
{
   check( ( i<0 || i>=3 ), "TTriangle::GetTT( int ), wrong index number" );
   
   return( TE[i]->ET[0] != this ? TE[i]->ET[0] : TE[i]->ET[1] );
}


 
// -----------------------------------------------------------------------------
//
//  void TTriangle::AddPoint( PTPoint )
//
//  Insert the given point in the PointList of (PTTriangle)this.
//  Compute the real error ot the point in the current triangulation:
//  if such error is larger than the error of the point at the head
//  of the PointList, then insert this point at the beginning;
//  else insert it in second position.
//

void TTriangle::AddPoint( PTPoint PointToAdd )
{
    double NewError = Geom::CalcError( this, PointToAdd );

    PointToAdd->Error = NewError;

    double OldPtErr = ( PointList.IsEmpty() ? 0.0 : PointList.GetHead()->Error );

    //
    // If the point being inserted is the one that causes the max error,
    // then insert it at the head of PointList; otherwise, insert it
    // in second position.
    //
    
    if ( NewError >= OldPtErr )
    {    
       PointList.AddHead( PointToAdd );
    }
    else
    {
       TListIterator<PTPoint> PLIter( &PointList );
       PLIter.Restart();
       PointList.AddAfter( PLIter.Current(), PointToAdd );
    }
}



// -----------------------------------------------------------------------------
//  
//   double TTriangle::GetError()
//  
//   Return the maximum error among the points in the PointList of this
//   triangle and of the edges adjacent to it.
//

double TTriangle::GetError()
{
   double errors[4];
   int i; /* PAOLA */

   for( i=0; i<3; i++ )
     if ( TE[i]->PointList.IsEmpty() ) 
       errors[i] = 0.0;
     else
       errors[i] = TE[i]->PointList.GetHead()->Error;

   if ( PointList.IsEmpty() ) errors[3] = 0.0 ;
    else errors[3] = PointList.GetHead()->Error;
    
   double error = errors[3];
   for( i=0; i<3; i++ )
      if ( errors[i] > error ) error = errors[i];
     
   return( error );

}


// -----------------------------------------------------------------------------
//  
//   void TTriangle::~TTriangle()
//
//   Destructor of TTriangle.
//   Make sure that possible entities pointing to this do not point 
//   to it any more. Empty the PointList.
//

TTriangle::~TTriangle()
{
   
   for( int e=0; e<3; e++ )
   {
      if ( this->TE[e] != NULL )
      {
         if ( this->TE[e]->ET[0] == this ) this->TE[e]->ET[0] = NULL;
	 else
	 if ( this->TE[e]->ET[1] == this ) this->TE[e]->ET[1] = NULL;
      }
   } 
   
   PointList.ClearList();
   
}



// -----------------------------------------------------------------------------
//
//  Output operator for class TTriangle
//

ostream& operator<< ( ostream& os, RTTriangle T )
{
   PTVertex v0, v1, v2;

   T.GetTV( v0, v1, v2 );

   os << "T" << T.TID << "[ V" << v0->VID << ", V" << v1->VID << ", V" << v2->VID;
   os << " ][ E" << T.TE[0]->EID << ", E" << T.TE[1]->EID << ", E" << T.TE[2]->EID;
   os << " ]";

   return os;

}


// -------------------------------------------------------------------------
//
//  int compare( PTPoint, PTPoint )
//
//  Compare two points based on the maximum error. Return -1, 0, +1 
//  according to the result of the comparison.
//  This function is required for using type TBTree<PTPoint>.
//  

int compare( PTPoint P0, PTPoint P1 )
{
   
   if ( P0->Error < P1->Error )           // T0 < T1
      return -1;
   else if ( P0->Error > P1->Error ) 	  // T0 > T1
      return  1;
   else
   {
   
       //
       // NB: it is necessary to sort well two points. Value 0 must
       // be returned only if the two points are perfectly equal.
       // In case error are equal, therefore, sort the two points based
       // on their coordinates.
       //

       int sx = Geom::Sign( P0->x - P1->x );
       if ( sx != 0 )  return sx;

       sx = Geom::Sign( P0->y - P1->y );
       if ( sx != 0 )  return sx;

       return Geom::Sign( P0->z - P1->z );       
   }
 
}   
      
