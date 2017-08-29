/* ----------------------------------------------------------------------- */
/*                               i_render.c                                */
/* ----------------------------------------------------------------------- */

/*
User interface for terrain decimation or refinement. Rendering.
*/

/* ----------------------------------------------------------------------- */

#include "i_global.h"

/* ----------------------------------------------------------------------- */

GLfloat black[3] = {0.0,0.0,0.0};
GLfloat red[3] = {1.0,0.0,0.0};
GLfloat yellow[3] = {1.0,1.0,0.0};
GLfloat green[3] = {0.0,1.0,0.0};
GLfloat blue[3] = {0.0,0.0,1.0};
GLfloat magenta[3] = {1.0,0.0,1.0};
GLfloat white[3] = {1.0,1.0,1.0};
GLfloat cyan[3] = {0.0,1.0,1.0};
GLfloat red1[3] = {1.0,0.7,0.2}; /* orange */
GLfloat yellow1[3] = {1.0, 1.0, 0.2}; /* gold yellow */
GLfloat green1[3] = {0.8, 1.0, 0.8}; /* light green */

int window_w = 500, window_h = 500; /* window width and height */
int is_help = 1;

/*
There are two movements: translation and rotation, and two control 
devices: arrow keys and mouse. This boolean variable tells which movement
is controlled by arrow keys, the other one is controlled by mouse.
*/
int is_rotation = 0; /* 1 iff rotation controlled by arrow keys */

/* ----------------------------------------------------------------------- */

unsigned int num_pts0 = 0;     /* initial number of points */
float * pts;                   /* array of initial points */

/* (used just in decimation) */
unsigned int num_itri0 = 0;    /* initial number of triangles */
int * itri;                    /* array of initial triangles */

unsigned int num_cst0 = 0;     /* initial number of constraints */
int * cst;                     /* array of constraints */

unsigned int num_vert0 = 0;    /* initial number of vertices */
unsigned int num_vert = 0;     /* number of vertices */
float * vert;                  /* array of vertices */

unsigned int num_tri0 = 0;     /* initial number of triangles */
unsigned int num_tri = 0;      /* number of triangles */
int * tri;                     /* array of triangles */

unsigned int num_segm0 = 0;    /* initial number of constraint segments */
unsigned int num_segm = 0;     /* number of constraint segments */
int * segm;                    /* array of constraint segments */

float minX, minY, minZ;
float maxX, maxY, maxZ;

#define PT_X(v) (pts[v*3])
#define PT_Y(v) (pts[v*3+1])
#define PT_Z(v) (pts[v*3+2])

#define CT_V1(c) (cst[2*c])
#define CT_V2(c) (cst[2*c+1])

#define IT_V1(t) (itri[3*t])
#define IT_V2(t) (itri[3*t+1])
#define IT_V3(t) (itri[3*t+2])

#define VERT_X(v) (vert[v*3])
#define VERT_Y(v) (vert[v*3+1])
#define VERT_Z(v) (vert[v*3+2])

#define TRIANG_V1(t) (tri[t*3])
#define TRIANG_V2(t) (tri[t*3+1])
#define TRIANG_V3(t) (tri[t*3+2])

#define SEGM_V1(s) (segm[s*2])
#define SEGM_V2(s) (segm[s*2+1])

int keep_ratio = 1; /* maintain proportions among x,y,z */

/* ----------------------------------------------------------------------- */

int points_on = 1;        /* Boolean, draw input points */
int constraints_on  = 1;  /* Boolean, draw input constraints */
int triangles_on = 1;     /* Boolean, draw triangles */
int segments_on = 1;      /* Boolean, draw constraint segments */
int vertices_on = 1;      /* Boolean, draw vertices */

float zoom_factor = 1.0;
float zoom_z = 0.3; //1.0; /* scale factor for rendering z over x,y */

/* rotation angles */
float x_angle = 0.0;
float y_angle = 0.0;

#define ANGLE_STEP 5.0

/* translation movements */
float x_move = 0.0;
float y_move = 0.0;
float z_move = 0.0;

/* drawing style */
int style = GL_LINE;

/* 2D or 3D view: 2 or 3 for 2d and 3d rendering */
int dimensions = 2;

int lastX, lastY;

/* ----------------------------------------------------------------------- */

void ChangeStyle()
{  if (style==GL_LINE) style = GL_FILL; else style = GL_LINE;  }

/* ----------------------------------------------------------------------- */

void SetLimits()
{
  int i;

  minX = minY = minZ = 0.0;
  maxX = maxY = maxZ = 1.0;
  for (i=0;i<num_pts0;i++)
  {
    if ((i==0) || (PT_X(i)>maxX)) maxX = PT_X(i);
    if ((i==0) || (PT_Y(i)>maxY)) maxY = PT_Y(i);
    if ((i==0) || (PT_Z(i)>maxZ)) maxZ = PT_Z(i);
    if ((i==0) || (PT_X(i)<minX)) minX = PT_X(i);
    if ((i==0) || (PT_Y(i)<minY)) minY = PT_Y(i);
    if ((i==0) || (PT_Z(i)<minZ)) minZ = PT_Z(i);
  }
  if (minX==maxZ) maxX += 1.0;
  if (minY==maxZ) maxY += 1.0;
  if (minZ==maxZ) maxZ += 1.0;
}

int LoadPoints(char * filename, int and_triangles, int and_edges)
{
  FILE * fd = fopen(filename,"r");
  int n, i;
  float x, y, z;
  int triangles_found = 0;
    
  if (!fd) printf("Sorry, cannot open input file %s\n", filename);
  if (!fd) return 0;
  if ( fscanf(fd,"%d",&n) == EOF ) return 0; /* nothing */
  num_pts0 = n;
  pts = (float*) malloc (3*num_pts0*sizeof(float));
  for (i=0;i<num_pts0;i++)
  {
    fscanf(fd,"%f %f %f",&PT_X(i), &PT_Y(i), &PT_Z(i));
  }  
  if (and_triangles)
  {
    if ( fscanf(fd,"%d",&n) == EOF ) return 1; /* just points */
    num_itri0 = n;
    itri = (int*) malloc (3*num_itri0*sizeof(int));
    for (i=0;i<num_itri0;i++)
    {
      fscanf(fd,"%d %d %d",&IT_V1(i), &IT_V2(i), &IT_V3(i));
    }
    if (!and_edges) {  fclose(fd); return 2;  } /* points and triangles */
    else triangles_found = 1;
  }
  if (and_edges)
  {
    if ( fscanf(fd,"%d",&n) == EOF ) return 1; /* just points */
    num_cst0 = n;
    cst = (int*) malloc (2*num_cst0*sizeof(int));
    for (i=0;i<num_cst0;i++)
    {
      fscanf(fd,"%d %d",&CT_V1(i), &CT_V2(i));
    }
    fclose(fd);
    return 2+triangles_found; /* points and constraints (and triangles) */
  }
  fclose(fd);
  return 1+triangles_found; /* just points (and triangles) */
}  

/* ----------------------------------------------------------------------- */

void SetLight(void)
{
  static GLfloat pos0[4] = {0.0, 1.0, 5.0, 0.0 };
  static GLfloat ambient0[4] = { 1.0, 1.0, 0.0, 1.0 };
  static GLfloat diffuse0[4] = { 0.6, 0.6, 0.6, 1.0 };
  static GLfloat pos1[4] = {2.0, 4.0, -5.0, 0.0 };
  static GLfloat ambient1[4] = { 1.0, 0.0, 0.0, 1.0 };
  static GLfloat diffuse1[4] = { 0.6, 0.6, 0.6, 1.0 };
  glLightfv( GL_LIGHT0, GL_POSITION, pos0 );
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
  glLightfv( GL_LIGHT1, GL_POSITION, pos1 );
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse1);
}

/* ----------------------------------------------------------------------- */

/* normal vector to a triangle (note: it is not normalized) */
void TriangleNormal9f(float x1, float y1, float z1,
                      float x2, float y2, float z2,
                      float x3, float y3, float z3,
                      float *x, float *y, float *z)
{
  double  a[3];
  double  b[3];
  a[0] = x1-x2; a[1] = y1-y2; a[2] = z1-z2;
  b[0] = x1-x3; b[1] = y1-y3; b[2] = z1-z3;
  (*x) = a[1]*b[2] - a[2]*b[1];
  (*y) = a[2]*b[0] - a[0]*b[2];
  (*z) = a[0]*b[1] - a[1]*b[0];
}

void SetTriangle9f(float x1, float y1, float z1,
                   float x2, float y2, float z2,	
                   float x3, float y3, float z3)
{
  if (dimensions==3)
  {
    float nx, ny, nz;
    TriangleNormal9f(x1,y1,z1*zoom_z,x2,y2,z2*zoom_z,x3,y3,z3*zoom_z,
                    &nx,&ny,&nz);
    glNormal3f(nx,ny,nz);
    glVertex3f(x1,y1,z1*zoom_z);
    glVertex3f(x2,y2,z2*zoom_z);
    glVertex3f(x3,y3,z3*zoom_z);
  }
  else /* ==2 */
  {
    glVertex2f(x1,y1);
    glVertex2f(x2,y2);
    glVertex2f(x3,y3);
  }
}

void SetSegment6f(float x1, float y1, float z1,
                  float x2, float y2, float z2)
{
  if (dimensions==3)
  {
    glVertex3f(x1,y1,z1*zoom_z);
    glVertex3f(x2,y2,z2*zoom_z);
  }
  else /* ==2 */
  {
    glVertex2f(x1,y1);
    glVertex2f(x2,y2);
  }
}

void SetScene()
{
  int i;

  if (triangles_on)
  {
    if (dimensions==3)
    {  glEnable(GL_LIGHTING);
       glMaterialfv( GL_BACK, GL_AMBIENT_AND_DIFFUSE, green1);
       glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, yellow1);
    }
    glColor3fv(green1);         
    glBegin(GL_TRIANGLES);
    if (num_tri)
      for (i=0;i<num_tri;i++)
      {
        SetTriangle9f
         (VERT_X(TRIANG_V1(i)),VERT_Y(TRIANG_V1(i)),VERT_Z(TRIANG_V1(i)),
          VERT_X(TRIANG_V2(i)),VERT_Y(TRIANG_V2(i)),VERT_Z(TRIANG_V2(i)),
          VERT_X(TRIANG_V3(i)),VERT_Y(TRIANG_V3(i)),VERT_Z(TRIANG_V3(i)));
      }
    else 
      for (i=0;i<num_itri0;i++)
      {
        SetTriangle9f
         (PT_X(IT_V1(i)),PT_Y(IT_V1(i)),PT_Z(IT_V1(i)),
          PT_X(IT_V2(i)),PT_Y(IT_V2(i)),PT_Z(IT_V2(i)),
          PT_X(IT_V3(i)),PT_Y(IT_V3(i)),PT_Z(IT_V3(i)));             
      }
    glEnd();
  }
  glDisable(GL_LIGHTING);
  if (constraints_on)
  {
    glColor3fv(red1);
    glBegin(GL_LINES);
    for (i=0;i<num_cst0;i++)
    {
      SetSegment6f(PT_X(CT_V1(i)), PT_Y(CT_V1(i)), PT_Z(CT_V1(i)),
                   PT_X(CT_V2(i)), PT_Y(CT_V2(i)), PT_Z(CT_V2(i)));
    }
    glEnd();
  }
  if (segments_on)
  {
    glColor3fv(red);
    glBegin(GL_LINES);
    for (i=0;i<num_segm;i++)
    {
      SetSegment6f(VERT_X(SEGM_V1(i)),VERT_Y(SEGM_V1(i)),VERT_Z(SEGM_V1(i)),
                   VERT_X(SEGM_V2(i)),VERT_Y(SEGM_V2(i)),VERT_Z(SEGM_V2(i)));
    }
    glEnd();
  }
  if (points_on)
  {
    glColor3fv(white);
    glBegin(GL_POINTS);
    for (i=0;i<num_pts0;i++)
    {
      if (dimensions==3)
        glVertex3f(PT_X(i),PT_Y(i),PT_Z(i)*zoom_z);
      else
        glVertex2f(PT_X(i),PT_Y(i));
    }    
    glEnd();
  }
  if (vertices_on)
  {
    glColor3fv(yellow);
    glBegin(GL_POINTS);
    for (i=0;i<num_vert;i++)
    {
      if (dimensions==3)
        glVertex3f(VERT_X(i),VERT_Y(i),VERT_Z(i)*zoom_z);
      else
        glVertex2f(VERT_X(i),VERT_Y(i));
    }    
    glEnd();
  }
}

void ViewScene(void)
{
  /* projection transformations */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  /* modeling trnasformations */
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if (dimensions==3) SetLight();

  glTranslatef(x_move,y_move,z_move);
  if (dimensions==3)
    glScalef(zoom_factor*0.8,zoom_factor*0.8,zoom_factor*0.8*zoom_z);
  else 
    glScalef(zoom_factor*0.8,zoom_factor*0.8,1.0);
  if (dimensions==3)
  {
    glRotatef(x_angle,1.0,0.0,0.0);
    glRotatef(y_angle,0.0,1.0,0.0);
    /* first put terrain in vertical position */
    glRotatef(-90.0,1.0,0.0,0.0);
  }

  /* scale to fit into window */
  if (keep_ratio)
  {
     float max_dim = (maxX-minX);
     if ((maxY-minY) > max_dim) max_dim = (maxY-minY);
     if (dimensions==3)
     {
       if ((maxZ-minZ)*zoom_z > max_dim) max_dim = (maxZ-minZ)*zoom_z;   
       glScalef(2.0/max_dim, 2.0/max_dim, 2.0/max_dim);
     }
     else 
       glScalef(2.0/max_dim, 2.0/max_dim, 1.0);
  }
  else
  {  
     if (dimensions==3)
       glScalef(2.0/(maxX-minX), 2.0/(maxY-minY), 2.0/(maxZ-minZ)*zoom_z);
     else
       glScalef(2.0/(maxX-minX), 2.0/(maxY-minY), 1.0);
  }
  /* translate to move center of scene to origin */
  if (dimensions==3)
    glTranslatef(-0.5*(maxX+minX), -0.5*(maxY+minY), -0.5*(maxZ+minZ)*zoom_z);
  else
    glTranslatef(-0.5*(maxX+minX), -0.5*(maxY+minY), 0.0);
}

/* ----------------------------------------------------------------------- */

void DeviceCoords(void)
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0,window_w,0,window_h);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

/* ----------------------------------------------------------------------- */

static int row;

void displayString(int x, int y, char * s)
{
  int i;
  glRasterPos2f((float)x,(float)y);
  for (i=0; s[i]!='\0'; i++)
  {
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, s[i]);
  }
}

void displayStringNextRow(int x, char * s)
{
  int i;
  glRasterPos2f((float)x,(float)row);
  for (i=0; s[i]!='\0'; i++)
  {
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, s[i]);
  }
  row -= 15;
}

void Footnote(void)
{
  displayString( 40,5, "Press * to switch help/graphic mode");
}

/* ----------------------------------------------------------------------- */

/* Used as display callback when help screen is active */
void HelpCallBack(void)
{
  int i, j;
  char aux[255];    

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glClearColor(black[0], black[1], black[2], 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  DeviceCoords();
  glColor3fv(white);
  row = window_h-20;
  displayStringNextRow(0, "Stages:");
  for (i=1; i<=NumberOfStages(); i++)
  {
    for (j=0; j<NumberHelpForStage(i); j++)
       displayStringNextRow(10, HelpForStage(i,j));
  }
  displayStringNextRow(0, CurrentStageString());
  for (i=1; i<=current_stage; i++)
  {
    for (j=0; j<NumberLinesForStage(i); j++)
    {  if (j==0) row -= 15;
       displayStringNextRow(0, LineForStage(i,j));
    }
  }
  row -= 15;
  if (DataLoaded())
  {
    sprintf(aux,"Initial %d points, z ranges from %g to %g", 
            num_pts0, minZ, maxZ);
    displayStringNextRow( 10, aux);
    row -= 15;
  }
  if (AlgoIsRunning())
  {
    sprintf(aux,"Performed %d updates",step_count);
    displayStringNextRow( 0, aux);
    sprintf(aux,"Initial %d vertices, now %d",num_vert0, num_vert);
    displayStringNextRow( 0, aux);
    sprintf(aux,"Initial %d triangles, now %d",num_tri0, num_tri);
    displayStringNextRow( 0, aux);
    sprintf(aux,"Initial %d segments, now %d",num_segm0, num_segm);
    displayStringNextRow( 0, aux);
  }
  Footnote();
}

void EchoCallBack(char * message)
{
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  DeviceCoords();
  glColor3fv(white);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glBegin(GL_QUADS);
    glVertex2i(window_w/5, 2*window_h/5);
    glVertex2i(4*window_w/5, 2*window_h/5);
    glVertex2i(4*window_w/5, 4*window_h/5);
    glVertex2i(window_w/5, 4*window_h/5);
  glEnd();
  glColor3fv(blue);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
    glVertex2i(window_w/5, 2*window_h/5);
    glVertex2i(4*window_w/5, 2*window_h/5);
    glVertex2i(4*window_w/5, 4*window_h/5);
    glVertex2i(window_w/5, 4*window_h/5);
  glEnd();
  glColor3fv(blue);
  displayString( window_w/5+10,4*window_h/5-20, message);
  displayString( window_w/5+10,4*window_h/5-35, written_input);
}

/* ----------------------------------------------------------------------- */

/* Used as display callback when rendering screen is active */
void GraphicCallBack(void)
{
  int i;
  
  if (dimensions==3)
  {
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_NORMALIZE);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
  }
  else
  {
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
  }

  glClearColor(blue[0], blue[1], blue[2], 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ViewScene();
    glPointSize(4.0);
    glLineWidth(1.0);
    if (dimensions==3) glPolygonMode(GL_FRONT_AND_BACK,style);
    else  glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
  SetScene();

  DeviceCoords();
  glColor3fv(white);
  if (dimensions==3)
  {
    if (is_rotation)
      displayString( 10,20, "Arrows rotate, mouse translates. Press F1 to switch.");
    else
      displayString( 10,20, "Arrows translate, mouse rotates. Press F1 to switch.");
    displayString( 10,35, "Press F2 to switch solid/wireframe");      
  }
  Footnote();
}

/* ----------------------------------------------------------------------- */

void DisplayCallBack(void)
{
  if (is_help) HelpCallBack(); else GraphicCallBack();
  if (writing_input) EchoCallBack(input_description);
  glFlush();
  glutSwapBuffers();
}

/* ----------------------------------------------------------------------- */

void NewRotation(float * angle, float delta)
{
  (*angle) += delta;
  if ((*angle)<-180.0) (*angle) += 360.0;
  if ((*angle)>180.0) (*angle) -= 360.0;
}

void ChangeRotation(int axis, int sign)
{
  if (axis==1) NewRotation(&x_angle, sign*ANGLE_STEP);
  else NewRotation(&y_angle, sign*ANGLE_STEP);
}

void RestoreRotation(void)
{  x_angle = y_angle = 0.0;  }

/* ----------------------------------------------------------------------- */

void NewTranslation(float * move, float delta) {  (*move) += delta;  }

void ChangeTranslation(int axis, int sign)
{
  if (axis==1) NewTranslation(&x_move, 0.1*sign);
  else NewTranslation(&y_move, 0.1*sign);
}

void RestoreTranslation(void)
{  x_move = y_move = 0.0;  }

/* ----------------------------------------------------------------------- */

void DefaultKeyboardCallBack(unsigned char c, int x, int y)
{
  switch (c)
  {
    case '+': zoom_factor += 0.1; break;
    case '-': zoom_factor -= 0.1; break;
    case '=': zoom_factor = 1.0; break;
    case 'p': 
    case 'P': points_on = !points_on; break;
    case 'c':
    case 'C': constraints_on = !constraints_on; break;
    case 'v': 
    case 'V': vertices_on = !vertices_on; break;
    case 'T':
    case 't': triangles_on = !triangles_on; break;
    case 'S':
    case 's': segments_on = !segments_on; break;
    case 27: /* ESC */
    case 'Q':
    case 'q': exit(0);
    case '*': is_help = !is_help; break;
    case '2': dimensions = 2; break;
    case '3': dimensions = 3; break;
    case '<': zoom_z *= 2.0; printf("Now zoom_z=%g\n", zoom_z); break;
    case '>': zoom_z *= 0.5; printf("Now zoom_z=%g\n", zoom_z); break;
  }
}

/* ----------------------------------------------------------------------- */

void SpecialCallBack(int c, int x, int y)
{
  switch (c)
  {
    case 1: /* F1 */ is_rotation = !is_rotation; break;
    case 2: /* F2 */ ChangeStyle(); break;
    case 100: /* left arrow */
      if (is_rotation) ChangeRotation(2,-1);
      else  ChangeTranslation(1,-1);
      break;
    case 101: /* up arrow */
      if (is_rotation) ChangeRotation(1,1);
      else ChangeTranslation(2,1);
      break;
    case 102: /* right arrow */
      if (is_rotation) ChangeRotation(2,1);
      else ChangeTranslation(1,1);
      break;
    case 103: /* down arrow */
      if (is_rotation) ChangeRotation(1,-1);
      else ChangeTranslation(2,-1); 
      break;
    case 106: /* home */
      RestoreRotation(); RestoreTranslation(); break;
  }
  glutPostRedisplay();
}

/* ----------------------------------------------------------------------- */

void MouseCallBack(int button, int state, int x, int y)
{
  if ( (button==GLUT_LEFT_BUTTON) && (state == GLUT_DOWN) )
  {  lastX = x; lastY = y;  }
}

void MotionCallBack(int x, int y)
{
  if (is_rotation)
  /* rotation controlled by arrow keys, translation by mouse */
  {
    NewTranslation(&x_move, 0.05*(x-lastX));
    NewTranslation(&y_move, -0.05*(y-lastY));
  }
  else /* rotation controlled by mouse, translation by arrow keys */
  {
    NewRotation(&y_angle, (0.05*(float)(x-lastX))*ANGLE_STEP);
    NewRotation(&x_angle, (0.05*(float)(y-lastY))*ANGLE_STEP);
  }
  lastX = x; lastY = y;
  glutPostRedisplay();
}

/* ----------------------------------------------------------------------- */

void ReshapeCallBack(int w, int h)
{
    int d = ( (w>h) ? h : w );
    glViewport(0, 0, d, d);
}

/* ----------------------------------------------------------------------- */
