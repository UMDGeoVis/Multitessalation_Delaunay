
#include "stepbystep.h"
#include "ttriangulation.h"

StepByStep::StepByStep(TTriangulation* triangulation)
{
   this->triangulation = triangulation;
}


StepByStep::~StepByStep()
{
   //triangulation = 0;
}


void StepByStep::ReadData( const char *ch )
{
   triangulation->ReadData(ch);
}


void StepByStep::WriteData( const char *ch )
{
    triangulation->WriteData(ch);
}

void StepByStep::ConvertData(unsigned int *vNum, unsigned int *tNum,
                             unsigned int *eNum,
                             float **vData, int **tData, int **eData)
{
    int vv, tt, ee;
    triangulation->ConvertData(&vv,&tt,&ee,vData,tData,eData);
    *vNum = vv; *tNum = tt; *eNum = ee;
}

void StepByStep::InitialTriangulation()
{
    triangulation->InitialTriangulation();
}


boolean StepByStep::NoMoreUpdates()
{
   return  triangulation->NoMoreUpdates();
}


void StepByStep::UpdateStep()
{
    triangulation->UpdateStep();
    if (triangulation->NoMoreUpdates()) triangulation->PrepareToEnd();
}



void StepByStep::EndTriangulation()
{
    triangulation->EndTriangulation();
}
