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

// -------------------------------------------------------------------------
//
//   file   : refrnddel.h
//   author : Christian Melchiorre
//
//   Definition of class TRefRndDelaunay, sub-class of TRefineDelaunay,
//   for building a Delaunay triangulation through refinement, with a
//   random choice of the point to be inserted.
//   The only difference With respect to the base class, is that points
//   not belonging to the convex hull are sorted based on a random key.
//   Thus, method InitialTriangulation() ius redefined: it performs
//   such random sorting after the operation that it performed in the 
//   base class. 
//


#ifndef _REFRNDDEL_H
#define _REFRNDDEL_H

#include "refdel.h"
#include "mttracer.h"


class TRefRndDelaunay : public TRefineDelaunay
{
   private:
   
      virtual void InitialTriangulation();
      
      void RandomSortPoints();

   public:

      TRefRndDelaunay ( PMTTracer iMT ) : TRefineDelaunay( iMT ) {};

};

#endif // _REFRNDDEL_H
