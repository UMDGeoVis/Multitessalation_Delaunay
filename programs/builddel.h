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
//   file   : builddel.h
//   author : Christian Melchiorre
//
//   Definition of class TBuildDelaunay, for the construcyion of a 
//   Delaunay triangulation starting from a list of input points
//   (construction through refinement).
//


#ifndef _BUILDDEL_H
#define _BUILDDEL_H



//
// include files
//

#include "defs.h"
#include "error.h"
#include "tdoublelist.h"

#include "ttriang.h"
#include "ttriangulation.h"
#include "basedel.h"



// ----------------------------------------------------------------------
//
//   class TBuildDelaunay
//


class TBuildDelaunay;

typedef class TBuildDelaunay *PTBuildDelaunay;
typedef class TBuildDelaunay &RTBuildDelaunay;


// Construction of a Delaunay triangulation given the input points.
// Initially, the convex hull of the point set is computed and 
// triangulated. Then a process of refinement will insert all other
// input points.

class TBuildDelaunay : virtual public TDelaunayBase
{

   protected:
 
      
      //
      // Methods
      //
      
      virtual void ReadData( const char * );

      virtual void InitialTriangulation();
       
      // There are no more updates to be done if the list of
      // points to be inserted is empty, i.e., if iNextPoint >= nPts.
      virtual boolean NoMoreUpdates()
        { return( iNextPoint >= nPts ); };

      
      // One update step.
      virtual void UpdateStep();
        virtual void NextPoint();
        virtual void InsertVertex();

      
      // Compute the influence region of a point to be inserted.
      virtual void CalcInfluenceRegion();
         virtual void InitInflRegnInternal();
         virtual void InitInflRegnExternal();
         virtual void CalcInflRegnMain();
      
      // Retriangulate the influence region.
      virtual void RetriangulateInfluenceRegion();
            

   public:

      TBuildDelaunay();
      
};


#endif // _BUILDDEL_H

