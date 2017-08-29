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
//   file   : basedel.h
//
//   Definition of class TDelaunayBase. Used as an abstract base class,
//   it defines common attributes and methods for all classes which
//   implement Delaunay triangulations.
//


#ifndef _BASEDEL_H
#define  _BASEDEL_H

#include "defs.h"
#include "ttriang.h"
#include "tdoublelist.h"
#include "ttriangulation.h"


class TDelaunayBase;

typedef class TDelaunayBase *PTDelaunayBase;
typedef class TDelaunayBase &RTDelaunayBase;

// Specific abstract base class for Delaunay triangulations.
// Abstract base class used to define common attributes and methods of
// classes which build Delaunay triangulations (such attributes and methods
// cannot be defined in TTriangulation because they are specific of
// the Delaunay criterion).

class TDelaunayBase : public TTriangulation
{
    protected:
      
      //
      // Status variables
      //
      
      //
      // Variables used in the process of determining the influence 
      // region and in its retriangulation.
      //
      
      // List containing the boundary edges of the influence region
      // of a vertex, sorted in counterclockwise order.
      TDoubleList<PTEdge> InflRegnBorder;
      
      // Pointer to an arbitrary triangle of a connected set of triangles
      // marked as TO_DELETE.
      // Used as an entry point for such triangle set, in order to
      // visit the set.
      PTTriangle FirstTrgToDel;
 
      // virtual void DeleteInfluenceRegion();
      // this function has been temporarily made public in order
      // to make it visible to TDecCDT::ExtDeleteInfluenceRegion
      // under the REDH52 compiler.

    public:
    
       virtual void DeleteInfluenceRegion();
    
       TDelaunayBase();
};

#endif  // _BASEDEL_H
