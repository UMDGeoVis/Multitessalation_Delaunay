/* ----------------------------------------------------------------------- */
/*                              i_common.c                                 */
/* ----------------------------------------------------------------------- */

/*
User interface for terrain refinement. Common functions for both cases.
*/

/* ----------------------------------------------------------------------- */

#include "i_global.h"

/* ----------------------------------------------------------------------- */

/* Update data to be displayed after a change of the 
   triangulation in the simplification process. */

void Change()
{
  algo->ConvertData(&num_vert, &num_tri, &num_segm,
                    (float**)&vert, (int**)&tri, (int**)&segm);
  if (just_created)
  {
    num_vert0 = num_vert;
    num_tri0 = num_tri;
    num_segm0 = num_segm;
    just_created = 0;
  }
}

/* ----------------------------------------------------------------------- */

/* Generate help strings (more help strings depend on the specific process
   of refinement or decimation, see files i_ref.c and i_dec.c) */

static char aux_s[255];

char * CurrentStageString(void)
{
  sprintf(aux_s,"Current stage: %d", current_stage);
  return aux_s;
}

/* ----------------------------------------------------------------------- */

/* Input / output operations */

int CheckFileExtension(char * name, char * ext)
{
  int length = strlen(name);
  if (length<5) return 0;
  if (name[length-4] != '.') return 0;
  if (name[length-3] != ext[0]) return 0;
  if (name[length-2] != ext[1]) return 0;
  if (name[length-1] != ext[2]) return 0;
  return 1;
}

/* ----------------------------------------------------------------------- */

/* Simplification algorithm */

void NextStep()
{
  if ( !algo->NoMoreUpdates())
  {
    algo->UpdateStep();
    step_count++;
    Change();
  }
  /* else don't call algo->EndTriangulation();
     yet because after it is not possible to call
     ConvertData or WriteData, call it only
     before quitting the program */
  else GoToNextStage();
}

void AllSteps()
{
  while( !algo->NoMoreUpdates() )
  {
    algo->UpdateStep();
    step_count++;
  }
  Change();
  /* don't call algo->EndTriangulation();
     yet because after it is not possible to call
     ConvertData or WriteData, call it only
     before quitting the program */
  GoToNextStage();
}

void EndOfProgram(void)
{
   /* call this now because we are sure that we
      don't need to write or display the data */
   if (algo) algo->EndTriangulation();
   exit(0);
}

/* ----------------------------------------------------------------------- */

/* Menu management: each stage has its menu */

void GoToNextStage(void)
{
  glutDetachMenu(GLUT_RIGHT_BUTTON);
  current_stage++;
  glutSetMenu(menus[current_stage-1]);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

/* ----------------------------------------------------------------------- */

/* MAIN */

int main(int argc, char **argv)
{
  int i;
  
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(window_w, window_h);
  glutCreateWindow(ProgramTitle());

  for (i=0; i<NumberOfStages(); i++)
  {
    menus[i] = MakeMenuForStage(i+1);
  }
  current_stage = 1;
  glutSetMenu(menus[0]);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
    
  glutDisplayFunc(DisplayCallBack);
  glutKeyboardFunc(KeyboardCallBack);
  glutSpecialFunc(SpecialCallBack);
  glutReshapeFunc(ReshapeCallBack);
  glutMouseFunc(MouseCallBack);
  glutMotionFunc(MotionCallBack);
  glutMainLoop();
  
  return 1;
}

/* ----------------------------------------------------------------------- */
