/* ----------------------------------------------------------------------- */
/*                                i_dec.c                                  */
/* ----------------------------------------------------------------------- */

/*
User interface for terrain decimation. Main program.
*/

/* ----------------------------------------------------------------------- */

#include "i_global.h"

/* ----------------------------------------------------------------------- */

PTTriangulation tt;
MTTracer mtt;
StepByStep  *algo = NULL;

char output_file[255];
char input_file[255] = "no file";
int just_created;

/*
Parameters for mesh construction algorithm.
Default parameters are: 
triangulation is not constrained,
termination criterion is based on number of updates,
zero updates are to be made,
error level is meaningless (here set to zero),
norm for error calculation is maximum error over all triangles,
random choice of next vertex to be processed,
decimation is not simultaneous,
no bound on max degree of removable vertices,
do not extend optimitazion (CDT),
do not remove vertices with one incident constraint,
do not remove vertices with tree or more incident constraints.
*/
struct 
{
  int constrained;            /* Boolean, 0 for Delaunay, 1 for CDT */
  int termination_criterion;  /* TERM_NUPD or TERM_ERR */
  int num_upd;                /* number of updates to be made (TERM_NUPD) */
                              /* if -1 with TERM_NUPD means all updates */
  float error_level;          /* error level to be reached (TERM_ERR) */
  int norm;                   /* NORM_MAX or NORM_MED or NORM_SQM */
  int random;                 /* Boolean, random or error-driven */
  int simultaneous;           /* Boolean, simultaneous decimation */
  int degree;                 /* max degree of removable vertex */
                              /* if 0 means any degree */
  int extend_optimization;    /* allow extension of optimization (CDT) */
  int remove_one;             /* allow removing vertices with 1 constraint */
  int remove_three;           /* allow removing vertices with 3 or more */
} 
algo_params = {0, TERM_NUPD, 0, 0.0, NORM_MAX, 1, 0, 0, 0, 0};

int step_count = 0;

/* ----------------------------------------------------------------------- */

int PrepareToLoadScene(char * filename);

/* ----------------------------------------------------------------------- */

#define STAGE_TYPE 1
#define STAGE_DEG 2
#define STAGE_LOAD 3
#define STAGE_TERM 4
#define STAGE_RUN 5
#define STAGE_END 6
int current_stage = 1;  /* current stage of execution */

int DataLoaded(void) {  return (current_stage>STAGE_LOAD);  }
int AlgoIsRunning(void) {  return (current_stage>=STAGE_RUN);  }

int NumberOfStages(void) {  return STAGE_END;  }

int NumberHelpForStage(int stage)
{
  switch(stage)
  {
    case STAGE_TYPE: return 3;
    case STAGE_DEG:  return 1;
    case STAGE_LOAD: return 1;
    case STAGE_TERM: return 1;
    case STAGE_RUN:  return 1;
    case STAGE_END:  return 1;  
  }
}

/* helps are all on one line, but the function parameters allow         
   for multiple lines */
char * HelpForStage(int stage, int i)
{
  switch(stage)
  {
    case STAGE_TYPE:
      if (i==0) return "1) Set decimation type";
      else if (i==1) return "   and norm for evaluating the error";
      else return "   of the triangulation";
    case STAGE_DEG:
      return "2) Set additional simplification options";
    case STAGE_LOAD:
      return "3) Load data from file";
    case STAGE_TERM:
      return "4) Set termination conditions";
    case STAGE_RUN:
      return "5) Run simplification algorithm and write results";
    case STAGE_END:
      return "6) End";
  }
  return "";
}

static char aux_s[255];

char * CurrentSimplifString(void)
{
  if (algo_params.constrained)
  {
    if (algo_params.random)
      return "Mesh decimation type: Random CDT";
    else
      return "Mesh decimation type: Error-driven  CDT";
  }
  else
  {
    if (algo_params.simultaneous)
    {
       if (algo_params.random) 
         return "Mesh decimation type: simultaneous Random Delaunay";
       else
         return "Mesh decimation type: simultaneous Error-driven Delaunay";
    }
    else
    {
       if (algo_params.random)  
         return "Mesh decimation type: Random Delaunay (default)";
       else
         return "Mesh decimation type: Error-driven Delaunay";
    }
  }
}

char * CurrentTerminationString(void)
{
  if (algo_params.termination_criterion==TERM_NUPD)
  {
     if (algo_params.num_upd<0)
       sprintf(aux_s,"%s", "Perform all updates");
     else if (algo_params.num_upd==0)
       sprintf(aux_s, "Perform no update (default)");
     else
       sprintf(aux_s, "Perform %d updates", algo_params.num_upd);
  }
  else
  {
    sprintf(aux_s, "Reach error level %g", algo_params.error_level);
  }
  return aux_s;
}

char * CurrentNormalString(void)
{
  switch (algo_params.norm)
  {
    case NORM_MAX: return "Error norm: Maximum error (default)"; break;
    case NORM_MED: return "Error norm: Mean error"; break;
    case NORM_SQM: return "Error norm: Square mean error"; break;
  }
  return "****";
}

int NumberLinesForStage(int stage)
{
  switch(stage)
  {
    case STAGE_TYPE: return 2;
    case STAGE_DEG: if (algo_params.constrained) return 4; else return 1;
    case STAGE_LOAD: return 0;
    case STAGE_TERM: return 3;
    case STAGE_RUN:  return 0;
    case STAGE_END:  return 0;
  }
}

char * LineForStage(int stage, int i)
{
  switch(stage)
  {
    case STAGE_TYPE:
      if (i==0) return CurrentSimplifString();
      if (i==1) return CurrentNormalString();
      break;
    case STAGE_DEG:
      if (i==0) {  if (algo_params.degree>0)
                    sprintf(aux_s,"Max degree of removable vertices: %d",
                            algo_params.degree);
                   else 
                    sprintf(aux_s,"No limitation for the degree of removable vertices");
                   return aux_s;
                }
      if (i==1) {  if (algo_params.extend_optimization)
                       return "Extend optimization";
                   else return "Do not extend optimization (default)";
                }
      if (i==2) {  if (algo_params.remove_one)
                       return "Shorten open constraint chains";
                   else
                       return "Do not shorten open constraint chains (default)";
                }
      if (i==3) {  if (algo_params.remove_three)
                       return "Break closed constraint chains";
                   else
                       return "Do not break closed constraint chains (default)";
                }
      break;                              
    case STAGE_TERM:
      if (i==0)
      {  sprintf(aux_s,"Data loaded from %s", input_file); return aux_s;  }
      if (i==1) return "Termination condition:";
      if (i==2) return CurrentTerminationString();
      break;
  }
  return "";
}


/* ----------------------------------------------------------------------- */

char input_description[255];
int writing_input = 0; /* 1 iff user is writing input */
char written_input[255] = "";
int written_chars = 0;

void EchoKeyboardCallBack(unsigned char c, int x, int y)
{
  if (c==27) /* ESC */
  {
    writing_input = written_chars = 0;
    written_input[0] = '\0';
    if (current_stage==STAGE_TERM) /* restore default */
    {  algo_params.termination_criterion==TERM_NUPD;
       algo_params.num_upd = -1;
    }
  } 
  else if ((c==10)||(c==13)) /* input string is complete */
  {
    writing_input = written_chars = 0;
    switch (current_stage)
    {
      case STAGE_DEG: /* written string is vertex degree */
          sscanf(written_input, "%d", &algo_params.degree);
          if (algo_params.degree<3) 
          {
            strcpy(written_input,"Max degree must be >= 3, press any key");
            writing_input = 1;
            written_chars = 0;
          }
          break;
      case STAGE_LOAD: /* written string is input file name */
          strcpy(input_file,written_input);
          if (PrepareToLoadScene(input_file))
          {
            GoToNextStage();
            written_input[0] = '\0';
          } 
          else
          {
            strcpy(written_input,"SORRY: EMPTY FILE, press any key");
            writing_input = 1;
            written_chars = 0;
          }
          break;
      case STAGE_TERM: /* termination parameters */
          if (algo_params.termination_criterion==TERM_NUPD)
            sscanf(written_input, "%d", &algo_params.num_upd);
          else
            sscanf(written_input, "%f", &algo_params.error_level);
          break;
      case STAGE_RUN:
      case STAGE_END: /* written string is output file name */
          strcpy(output_file,written_input);
          algo->WriteData(output_file);
          break;
    }
  }
  else if (c==8) /* back space */
  {
    if (written_chars) written_input[--written_chars] = '\0';
  }
  else
  {
    written_input[written_chars++] = c;
    written_input[written_chars] = '\0';
  }
}

void KeyboardCallBack(unsigned char c, int x, int y)
{
  if (writing_input) EchoKeyboardCallBack(c,x,y);
  else DefaultKeyboardCallBack(c,x,y);
  glutPostRedisplay();
}

/* ----------------------------------------------------------------------- */

void ApplyAlgoParams()
{

   if (algo_params.termination_criterion == TERM_NUPD)
   {
     if (algo_params.num_upd == -1) /* all updates */
       mtt.SetTerminateCondition(TERM_NUPD, num_pts0);
     else
       mtt.SetTerminateCondition(TERM_NUPD, algo_params.num_upd);
   }
   else
   {
       mtt.SetTerminateCondition(TERM_ERR, 
                                 algo_params.norm,
                                 algo_params.error_level);
   }
}

/* ----------------------------------------------------------------------- */

/* Input / output operations */

int PrepareToLoadScene(char * filename)
{
   if (algo_params.constrained)
   {  if (!CheckFileExtension(filename,"cdt")) 
      {  fprintf(stderr,"Warning: Input file extension is not .cdt\n");
         fprintf(stderr,"Are you sure that it contains a constrained triangulation?\n");
      }         
   }
   else
   {  if (!CheckFileExtension(filename,"tri"))
      {  fprintf(stderr,"Warning: Input file extension is not .tri\n");
         fprintf(stderr,"Are you sure that it contains a triangulation?\n");
      }
   }
   if (!LoadPoints(filename, 1, algo_params.constrained)) return 0;
   SetLimits();
   return 1;
}

void LoadScene(char * filename)
{
   ApplyAlgoParams();
   if (algo_params.constrained)
   {
     if (algo_params.random)
       tt = new TDecRndCDT( algo_params.degree, &mtt,
                            (boolean)algo_params.extend_optimization,
                            (boolean)algo_params.remove_one,
                            (boolean)algo_params.remove_three);
     else
       tt = new TDecErrCDT( algo_params.degree, &mtt,
                            (boolean)algo_params.extend_optimization,
                            (boolean)algo_params.remove_one,
                            (boolean)algo_params.remove_three);
   }
   else
   {
     if (algo_params.simultaneous)
     {
       if (algo_params.random)
         tt = new TDecRndDeBerg( algo_params.degree, &mtt );
       else
         tt = new TDecErrDeBerg( algo_params.degree, &mtt );
     }
     else
     {
       if (algo_params.random)
         tt = new TDecRndDelaunay( algo_params.degree, &mtt );
       else
         tt = new TDecErrDelaunay( algo_params.degree, &mtt );
     }
   }
   if (tt == NULL) cout << "Error1 out of memory\n";
     algo = new StepByStep(tt);
   if (algo == NULL) cout << "Error2 out of memory\n";

   algo->ReadData(filename);
   num_itri0 = 0; if (itri) free(itri);
   algo->InitialTriangulation();
   just_created = 1;
   Change();
}

/* ----------------------------------------------------------------------- */

/* Menu management: each stage has its menu */

int menus[STAGE_END];
       
void StageOptionsTYPE(int v)
{
  switch (v)
  {
    case 1: /* Random Delaunay */
      algo_params.constrained = algo_params.simultaneous = 0;
      algo_params.random = 1; 
      break;
    case 2: /* Error-driven Delaunay */
      algo_params.constrained = algo_params.random =
                                algo_params.simultaneous = 0;
      break;
    case 3: /* Random simultaneous Delaunay */
      algo_params.constrained = 0;
      algo_params.random = algo_params.simultaneous = 1; 
      break;
    case 4: /* Error-driven simultaneous Delaunay */
      algo_params.constrained = algo_params.random = 0;
      algo_params.simultaneous = 1;
      break;
    case 5: /* Random CDT */
      algo_params.constrained = algo_params.random = 1;
      algo_params.simultaneous = 0;
      break;
    case 6: /* Error-driven CDT */
      algo_params.constrained = 1;
      algo_params.random = algo_params.simultaneous = 0;
      break;
    case 11: /* max error */ algo_params.norm = NORM_MAX; break;
    case 12: /* mean error */ algo_params.norm = NORM_MED; break;
    case 13: /* square mean error */ algo_params.norm = NORM_SQM; break;
    case 88: GoToNextStage(); break;
    case 99: EndOfProgram(); break;
  }
}

void StageOptionsDEG(int v)
{
  switch (v)
  {
    case 1: /* Any degree */
      algo_params.degree = 0; break;
    case 2: /* Bounded degree */
      /* version from command line:
      printf("Bound on vertex degree:");
      scanf("%d", &algo_params.degree); */
      /* version from window: keyboard callback will take care
         of collecting the input and performing the next tasks */
      writing_input = 1;
      strcpy(input_description, "Bound on vertex degree:");
      break;
    case 11: /* Extend optimization */
      algo_params.extend_optimization = !algo_params.extend_optimization;
      break;
    case 12: /* Remove vertices with one constraint */
      algo_params.remove_one = !algo_params.remove_one; break;
    case 13: /* Remove vertices with three or more constraint */
      algo_params.remove_three = !algo_params.remove_three; break;
    case 88: GoToNextStage(); break;
    case 99: EndOfProgram(); break;
  }
}

void StageOptionsLOAD(int v)
{
  switch (v)
  {
    case 1: /* load */ 
      /* Version from command line:
      printf("Name of file to be open:");
      scanf("%s", input_file);
      PrepareToLoadScene(input_file);
      don't break, after loading we go to next stage */
      /* version from window: keyboard callback will take care
         of collecting the input and performing the next tasks */
      writing_input = 1;
      strcpy(input_description, "Name of file to be open:");
      break;
    case 88: GoToNextStage(); break;
    case 99: EndOfProgram(); break;
  }
  glutPostRedisplay();
}

void StageOptionsTERM(int v)
{
  switch (v)
  {
    case 1: /* all updates */
      /* mark this case with a negative number, then it will
         be set to the number of vertices */
      algo_params.num_upd = -1;
      algo_params.termination_criterion = TERM_NUPD;
      break;
    case 2: /* number of updates */ 
      algo_params.termination_criterion = TERM_NUPD;
      /* version from command line:
      printf("Number of updates to be made:");
      scanf("%d", &algo_params.num_upd); */
      /* version from window: keyboard callback will take care
         of collecting the input and performing the next tasks */
      writing_input = 1;
      strcpy(input_description, "Number of updates to be made:");
      break;
    case 3: /* error level */
      algo_params.termination_criterion = TERM_ERR;
      /* version from command line:
      printf("Error level to be reached:");
      scanf("%f", &algo_params.error_level);*/
      /* version from window: keyboard callback will take care
         of collecting the input and performing the next tasks */
      writing_input = 1;
      strcpy(input_description, "Error level to be reached:");      
      break;
    case 88: 
      LoadScene(input_file);
      GoToNextStage();
      break;
    case 99: EndOfProgram(); break;
  }
  glutPostRedisplay();
}

void StageOptionsRUN(int v)
{
  switch (v)
  {
    case 1: /* next step */ NextStep(); break;
    case 2: /* all steps */ AllSteps(); break;
    case 3: /* write */
      /* version from command line:
      printf("Name of file to be written:");
      scanf("%s", output_file);
      algo->WriteData(output_file);*/
      /* version from window: keyboard callback will take care
         of collecting the input and performing the next tasks */
      writing_input = 1;
      strcpy(input_description, "Name of file to be written:");
      break;
    case 88: GoToNextStage(); break;
    case 99: EndOfProgram(); break;
  }
  glutPostRedisplay();
}

int MakeMenuForStage(int i)
{
  int m, m1, m2;
  switch (i)
  {
    case STAGE_TYPE: 
      m1 = glutCreateMenu(StageOptionsTYPE);
      glutAddMenuEntry("Random Delaunay",1);
      glutAddMenuEntry("Error-driven Delaunay",2);
      glutAddMenuEntry("Random Delaunay (simultaneously)",3);
      glutAddMenuEntry("Error-driven Delaunay (simultaneously)",4);
      glutAddMenuEntry("Random CDT",5);
      glutAddMenuEntry("Error-driven CDT",6);
      m2 = glutCreateMenu(StageOptionsTYPE);
      glutAddMenuEntry("maximum error",11);
      glutAddMenuEntry("mean error",12);
      glutAddMenuEntry("square mean error",13);
      m = glutCreateMenu(StageOptionsTYPE);
      glutAddSubMenu("mesh decimation type",m1);
      glutAddSubMenu("error norm",m2);
      glutAddMenuEntry("go to next stage",88);
      glutAddMenuEntry("exit",99);
      break;
    case STAGE_DEG:
      m1 = glutCreateMenu(StageOptionsDEG);
      glutAddMenuEntry("Any degree",1);
      glutAddMenuEntry("Set max degree",2);
      m2 = glutCreateMenu(StageOptionsDEG);
      glutAddMenuEntry("Allow extension of optimization(y/n)",11);
      glutAddMenuEntry("Shorten open constraint chains(y/n)",12);
      glutAddMenuEntry("Break closed constraint chains(y/n)",13);
      m = glutCreateMenu(StageOptionsDEG);
      glutAddSubMenu("Bound on vertex degree",m1);
      glutAddSubMenu("Change default CDT decimation options (CDT only)",m2);
      glutAddMenuEntry("go to next stage",88);
      glutAddMenuEntry("exit",99);
      break;
    case STAGE_LOAD:
      m = glutCreateMenu(StageOptionsLOAD);
      glutAddMenuEntry("load data from file",1);
      glutAddMenuEntry("exit",99);
      break;      
    case STAGE_TERM:
      m1 = glutCreateMenu(StageOptionsTERM);
      glutAddMenuEntry("all possible updates",1);
      glutAddMenuEntry("number of updates",2);
      glutAddMenuEntry("error level",3);
      m = glutCreateMenu(StageOptionsTERM);
      glutAddSubMenu("termination",m1);
      glutAddMenuEntry("go to next stage",88);
      glutAddMenuEntry("exit",99);
      break;      
    case STAGE_RUN:
      m = glutCreateMenu(StageOptionsRUN);
      glutAddMenuEntry("next step",1);
      glutAddMenuEntry("all steps",2);
      glutAddMenuEntry("write triangulation to file",3);
      glutAddMenuEntry("go to next stage",88);
      glutAddMenuEntry("exit",99);
      break;      
    case STAGE_END:
      m = glutCreateMenu(StageOptionsRUN);
      glutAddMenuEntry("write triangulation to file",3);
      glutAddMenuEntry("exit",99);
      break;      
  }
  return m;
}

/* ----------------------------------------------------------------------- */

/* Title of program */

char * ProgramTitle(void) {  return "Mesh decimation...";  }

/* ----------------------------------------------------------------------- */
