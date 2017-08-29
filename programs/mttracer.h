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

// -----------------------------------------------------------------------------------------
//
//  file   : mttracer.h
//  author : Christian Melchiorre
//
//  Definition of static class MTTracer. It acts as an interface between 
//  the Delaunay triangulation programs and the MT History Tracer for
//  building a multiresolution model (MT).
//  It implements the termination criteria for the process of triangulation
//  update, calls the functions of the MT library, and takes track of
//  the error of the current triangulation and of the mumber of updates 
//  performed so far.
//

#ifndef MTHIST_H
#define MTHIST_H

#include "defs.h"
#include "ttriang.h"
#include "tbtree.h"

#ifdef MT_TRACER
#include "fieldmt.h"
#else
/* it does not matter if we don't build the MT, they must be defined */
#define MT_REFINING 1
#define MT_COARSENING 2
#endif

#define ADDTRG      0
#define DELTRG      1

#define NORM_MAX    0   // maximum error among triangles
#define NORM_MED    1   // mean error among triangles
#define NORM_SQM    2   // squared mean error among triangles
 
#define NO_TERM     0   // no termination condition... (as long as updates arrive, axecute them)
#define TERM_NUPD   1   // termination based on number of updates
#define TERM_ERR    2   // termination based on error


class MTTracer;

typedef MTTracer *PMTTracer;
typedef MTTracer &RMTTracer;


#ifdef MT_TRACER
class MTTracer : public FieldWithErrorBuildingInterfaceClass
#else
class MTTracer
#endif
{
   private:
 
       // History type
       int Type; // MT_REFINING or MT_COARSENING
       
       //
       // Variables for maintaining the current triangulation error
       //
       
       double TotError;       // current triangulation error
       TBTree<double> Errors; // pnly used if norm = NORM_MAX       
       int nTrgs;             // number of triangles in the current triangulation

       //
       // Variables for termination
       //
              
       int Term;              // TERM_NUPD or TERM_ERR
       int Norm;              // type of norm
       int UpdLev;            // max number of updates
       double ErrLev;         // max/min error level       
       int nUpd;              // number of updates performed so far
      
   public:
   
       MTTracer(void);
       
       // Set termination contition.
       void SetTerminateCondition( int, int ); // terminazione per numero update
       void SetTerminateCondition( int, int, double ); // terminazione per livello di errore
       
       void StartHistory( int );
      
       // Record in MT History Tracer the insertion of a triangle
       // with its error.
       void MakeTriangle( PTTriangle );

       // Record in MT History Tracer the deletion of a triangle.
       void KillTriangle( PTTriangle );
      
       // Communicate to MT History Tracer that the update step is done
       // and the triangulation is in consistent state.
       void MeshOk();
      
       void EndHistory();
      
       // Update the triangulation error after the deletion or insertion 
       // of a triangle.
       void UpdateError( int, double );
       
       // Check if the termination condition is satisfied.
       boolean TerminateCondition();
       
};

#endif // MTHIST_H
