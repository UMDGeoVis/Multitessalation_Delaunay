/* ----------------------------------------------------------------------- */
/*                                i_ref.c                                  */
/* ----------------------------------------------------------------------- */

/*
User interface for terrain refinement. Main program.
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
random choice of next vertex to be processed.
*/
struct 
{
  int constrained;            /* Boolean, 0 for Delaunay, 1 for CDT */
  int termination_criterion;  /* TERM_NUPD or TERM_ERR */
  int num_upd;                /* number of updates to be made (TERM_NUPD) */
  float error_level;          /* error level to be reached (TERM_ERR) */
  int norm;                   /* NORM_MAX or NORM_MED or NORM_SQM */
  int random;                 /* Boolean, random or error-driven */
} 
algo_params = {0, TERM_NUPD, 0, 0.0, NORM_MAX, 1};

int step_count = 0;

/* ----------------------------------------------------------------------- */

int PrepareToLoadScene(char * filename);

/* ----------------------------------------------------------------------- */

#define STAGE_TYPE 1
#define STAGE_LOAD 2
#define STAGE_TERM 3
#define STAGE_RUN 4
#define STAGE_END 5
int current_stage = 1;  /* current stage of execution */

int DataLoaded(void) {  return (current_stage>STAGE_LOAD);  }
int AlgoIsRunning(void) {  return (current_stage>=STAGE_RUN);  }

int NumberOfStages(void) {  return STAGE_END;  }

int NumberHelpForStage(int stage)
{
  switch(stage)
  {
    case STAGE_TYPE: return 3;
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
      if (i==0) return "1) Set refinement type";
      else if (i==1) return "   and norm for evaluating the error";
      else return "   of the triangulation";
    case STAGE_LOAD:
      return "2) Load data from file";
    case STAGE_TERM:
      return "3) Set termination conditions";
    case STAGE_RUN:
      return "4) Run simplification algorithm and write results";
    case STAGE_END:
      return "5) End";
  }
  return "";
}

static char aux_s[255];

char * CurrentSimplifString(void)
{
  if (algo_params.constrained)
    return "Mesh refinement type: CDT";
  else
  {
    if (algo_params.random) 
      return "Mesh refinement type: Random Delaunay (default)";
    else 
      return "Mesh refinement type: Error-driven Delaunay";
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
   {  if (!CheckFileExtension(filename,"seg")) 
      {  fprintf(stderr,"Warning: Input file extension is not .seg\n");
         fprintf(stderr,"Are you sure that it contains points and segments?\n");
      }
   }
   else
   {  if (!CheckFileExtension(filename,"pts"))
      {  fprintf(stderr,"Warning: Input file extension is not .pts\n");
         fprintf(stderr,"Are you sure that it contains points?\n");
      }
   }
   if (!LoadPoints(filename, 0, algo_params.constrained)) return 0;
   SetLimits();
   return 1;
}

void LoadScene(char * filename)
{
   ApplyAlgoParams();
   if (algo_params.constrained)
   {
       tt = new TRefCDT( &mtt );
   }
   else
   {
     if (algo_params.random)
       tt = new TRefRndDelaunay( &mtt );
     else
       tt = new TRefErrDelaunay( &mtt );
   }
   if (tt == NULL) cout << "Error1 out of memory\n";
     algo = new StepByStep(tt);
   if (algo == NULL) cout << "Error2 out of memory\n";

   algo->ReadData(filename);
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
      algo_params.constrained = 0;
      algo_params.random = 1; 
      break;
    case 2: /* Error-driven Delaunay */
      algo_params.constrained = algo_params.random = 0;
      break;
    case 3: /* CDT */
      algo_params.constrained = 1;
      break;
    case 11: /* max error */ algo_params.norm = NORM_MAX; break;
    case 12: /* mean error */ algo_params.norm = NORM_MED; break;
    case 13: /* square mean error */ algo_params.norm = NORM_SQM; break;
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
      glutAddMenuEntry("CDT",3);
      m2 = glutCreateMenu(StageOptionsTYPE);
      glutAddMenuEntry("maximum error",11);
      glutAddMenuEntry("mean error",12);
      glutAddMenuEntry("square mean error",13);
      m = glutCreateMenu(StageOptionsTYPE);
      glutAddSubMenu("Mesh refinement type",m1);
      glutAddSubMenu("error norm",m2);
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
