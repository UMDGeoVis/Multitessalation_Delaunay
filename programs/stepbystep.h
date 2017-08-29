
#ifndef _STEP_BY_STEP_H
#define _STEP_BY_STEP_H


#include "defs.h" // necessary for type 'boolean'
#include "ttriangulation.h"

class TTriangulation;


//
// Wrapper class to make accessible some private methods of class
// TTriangulation.
// These are the methods which allow a step by step execution from
// an application outside this library.
//


class StepByStep
{
   public:

      StepByStep(TTriangulation* triangulation);
      ~StepByStep();

      void ReadData( const char * );
      void WriteData( const char * );
      void ConvertData(unsigned int *vNum, unsigned int *tNum, unsigned int *eNum,
                       float **vData, int **tData, int **eData);
               
      //
      //! Preparation before starting to update the initial triangulation.
      /*! In refinement, compute an initial triangulation (made by the vertices
          of the convex hull), set up array Points and, in TRefErrDelaunay, 
          put into PtsErrTree the points  which are still to be inserted.
          In decimation, insert into ElimVtxTree the removable vertices 
          whose degree fullfills the requirements.
      */
      void InitialTriangulation();

      //! Boolean function telling whether further updates are necessary.
      /*! Boolean function telling whether further updates are necessary
          in order to obtain the desired precision or the desired number 
          of updates.
      */
      boolean NoMoreUpdates();

      //
      //! This function executes a signle update step.
      /*! In refinement, select a vertex to be inserted, and insert it
          into the triangulation; in decimation, select a vertex to be
          removed and remove it from the triangulation.
      */
      void UpdateStep();

      //
      // This function performs work needed after the end of the
      // construction of the triangulation. By default, nothing.
      //
      void EndTriangulation();

      
    protected:

      TTriangulation* triangulation;
};


#endif // _STEP_BY_STEP_H
