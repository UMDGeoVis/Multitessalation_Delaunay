/* ----------------------------------------------------------------------- */
/*                               i_global.h                                */
/* ----------------------------------------------------------------------- */

/*
User interface for terrain decimation or refinement. Global declarations.
*/

/* ----------------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <GL/glut.h>

#include "mttracer.h"
#include "decrnddel.h"
#include "decerrdel.h"
#include "refdel.h"
#include "referrdel.h"
#include "refrnddel.h"
#include "decCDT.h"
#include "decrndCDT.h"
#include "decrnddbCDT.h"
#include "decerrCDT.h"
#include "refCDT.h"
#include "decrnddb.h"
#include "decerrdb.h"
#include "stepbystep.h"

/* ----------------------------------------------------------------------- */

/**** From file i_dec.c or i_ref.c ****/

extern PTTriangulation tt;
extern MTTracer mtt;
extern StepByStep  *algo;

extern int just_created;
extern char input_file[];

extern int current_stage; /* stage of execution of the interface */
extern int step_count;    /* number of simplification steps done */
extern int menus[]; /* every stage has its menu */

extern char input_description[];
extern int writing_input;    /* 1 iff user is writing input */
extern char written_input[];

extern char * ProgramTitle(void);
extern int DataLoaded(void);
extern int AlgoIsRunning(void);

extern int NumberOfStages(void);
extern int NumberHelpForStage(int stage);
extern char * HelpForStage(int stage, int i);
extern int NumberLinesForStage(int stage);
extern char * LineForStage(int stage, int i);

extern char * CurrentStageString(void);
extern int MakeMenuForStage(int i);

extern void KeyboardCallBack(unsigned char c, int x, int y);
   
/* Get initial data from file (build initial triangulation if needed) */
extern void LoadScene(char * filename);

/* Get vertices, triangles, segments from current triangulation */
extern void Change(void);

/* ----------------------------------------------------------------------- */

/**** From file i_common.c ***/

extern int CheckFileExtension(char * name, char * ext);
extern void NextStep();
extern void AllSteps();
extern void EndOfProgram(void);
extern void GoToNextStage(void);

/* ----------------------------------------------------------------------- */

/**** From file i_render.c ****/

extern int window_w, window_h; /* window width and height */

extern void SetLimits();
extern void ChangeStyle(); /* switch solid / wireframe */

/* Data to be rendered */

extern unsigned int num_pts0;   /* initial number of points */
extern float * pts;             /* array of initial points */

/*
Used just in decimation to contain the initial triangles,
when the triangulation has not yet been created.
Then, the array itri is deallocated.
*/
extern unsigned int num_itri0;  /* initial number of triangles */
extern int * itri;              /* array of initial triangles */

extern unsigned int num_cst0;   /* initial number of constraints */
extern int * cts;               /* array of constraints */

extern unsigned int num_vert0;  /* initial number of vertices */
extern unsigned int num_vert;   /* number of vertices */
extern float * vert;            /* array of vertices */

extern unsigned int num_tri0;   /* initial number of triangles */
extern unsigned int num_tri;    /* number of triangles */
extern int * tri;               /* array of triangles */

extern unsigned int num_segm0;  /* initial number of constraint segments */
extern unsigned int num_segm;   /* number of constraint segments */
extern int * segm;              /* array of constraint segments */

/* Get initial point set from file */
extern int LoadPoints(char * filename, int and_triangles, int and_edges);

/* Callbacks (more are callbacks are in files i_dec.c and i_ref.c) */

extern void DisplayCallBack(void);
extern void ReshapeCallBack(int w, int h);
extern void DefaultKeyboardCallBack(unsigned char c, int x, int y);
extern void SpecialCallBack(int c, int x, int y);
extern void MouseCallBack(int button, int state, int x, int y);
extern void MotionCallBack(int x, int y);

/* ----------------------------------------------------------------------- */
