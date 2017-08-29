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
//   file   : decrnddel.h
//   author : Christian Melchiorre
//
//   Definition of class TDecRndDelaunay, sub-class of TDecimDelaunay,
//   for decimation of a Delaunay triangulation with a random choice of
//   the vertex to be removed.
//


#ifndef _DECRNDDEL_H
#define _DECRNDDEL_H

#include "defs.h"
#include "ttriang.h"
#include "decdel.h"
#include "mttracer.h"

class TDecRndDelaunay;

typedef class TDecRndDelaunay *PTDecRndDelaunay;
typedef class TDecRndDelaunay &RTDecRndDelaunay;

class TDecRndDelaunay : virtual public TDecimDelaunay
{
   protected:
   /*
   #ifndef CC_SILICON
       virtual int srand48(int);
       
       virtual int drand48();
   #endif 
   */
       virtual void InitialTriangulation();

   public:
  
       TDecRndDelaunay( int iK, PMTTracer iMT ) : TDecimDelaunay( iK, iMT ) {};
};


#endif // _DECRNDDEL_H
