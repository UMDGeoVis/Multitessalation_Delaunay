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

// ----------------------------------------------------------------------------------------
//
//  file   : deberg.h
//  author : Christian Melchiorre
//
//  Definition of class DeBerg, which exports function SelectVertices.
//  Such function takes a set (structured as a tree) of removable vertices
//  in a triangulation, and returns a list of a subset of independent 
//  vertices in it.
//  Such function will be used by classes TDecRndDeBerg (sub-class of
//  TDecRndDelaunay) and TDecErrDeBerg (sub-class of TDecErrDelaunay).
//  rispettivamente, da TDecRndDelaunay e TDecErrDelaunay.
//  The solution of a static class avoids using multiple inheritance
//  with all problems connected to virtual base classes etc...
//


#ifndef _DEBERG_H
#define _DEBERG_H

#include "defs.h"
#include "ttriang.h"
#include "tbtree.h"
#include "tdoublelist.h"

class DeBerg
{
   private:
   
   public:
    
      static void MarkAllNeighbours( PTVertex, boolean *, int );
      static void SelectVertices( TBTree<PTVertex>&, TDoubleList<PTVertex> & );

};

#endif // _DEBERG_H
