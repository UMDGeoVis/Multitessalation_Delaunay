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
//  file   : markable.cpp
//   author : Christian Melchiorre
//
//  Implementation of class Markable used to define the behaviour of
//  objects which can be marked.
//

#include "defs.h"
#include "markable.h"

#include <iostream>

// -----------------------------------------------------------------------------
//
//  void Mark( MARKTYPE )
//
//  Mark this with value MarkValue
//

void Markable::Mark( MARKTYPE MarkValue )
{
    mark = mark | MarkValue;
}


// -----------------------------------------------------------------------------
//
//  void UnMark( MARKTYPE )
//
//  remove mark MarkValue from this
//

void Markable::UnMark( MARKTYPE MarkValue )
{
   mark &= ~MarkValue;
}


// -----------------------------------------------------------------------------
//
//  boolean Marked( MARKTYPE )
//
//  Check if this is marked with value MarkValue.
//


boolean Markable::Marked( MARKTYPE MarkValue )
{
   return( mark & MarkValue );
}



// -----------------------------------------------------------------------------
//
//  void MarkReset()
//
//  Remove all marks from this.
//

void Markable::MarkReset()
{
    mark = NULLMARK;
}



// -----------------------------------------------------------------------------
//
//  MARKTYPE MarkValue()
//
//  Return the content of field mark. Used only for debug purpose.
//

MARKTYPE Markable::MarkValue()
{
  return(mark);
}

