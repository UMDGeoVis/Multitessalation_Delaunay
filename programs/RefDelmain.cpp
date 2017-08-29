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

/* ------------------------------------------------------------------------ */
/*           ITERATIVE REFINEMENT OF A DELAUNAY TRIANGULATION               */
/* ------------------------------------------------------------------------ */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "error.h"
#include "geom.h"
#include "utils.h"

#include "mttracer.h"

#include "builddel.h"
#include "refdel.h"
#include "refrnddel.h"
#include "referrdel.h"

#define RANDOM 0
#define ERROR  1

int main( int argc, char **argv )
{
   char c = '\0';
   int nextpt;
   int interactive_mode = 0;
   int a; /* current argument scanned on command line, if not interactive */

   PTTriangulation T;
   MTTracer MT;

   cerr.precision(16);
   cerr << endl;

   if ( argc < 3 )
   {
      cerr << "usage: " << argv[0] << " infile outfile [format_string]" << endl; 
      exit(-1);
   }

   if ( argc == 3 ) interactive_mode = 1;
   else a = 3; /* for command line mode */

   //
   // Choice of the next point to be inserted: error-driven or random
   //

   if (interactive_mode)
      c = AskLetter("Choose points Error-based or Random", "eErR");
   else
      c = ParseLetter(argc, argv, a++, "eErR");
   cerr << "Next point option = " << c << endl;
   switch( c )
   {
      case 'E': case 'e': nextpt = ERROR;  break;
      case 'R': case 'r': nextpt = RANDOM; break;
   }

   //
   // Termination of the update sequence
   //

   if (interactive_mode)
      c = AskLetter("Termination based on Update number, Error level, All updates",
                    "uUeEaA");
   else
      c = ParseLetter(argc, argv, a++, "uUeEaA");
   cerr << "Termination option = " << c << endl;
   switch( c )
   {
      case 'U': case 'u':
      {
         int nupd;
         if (interactive_mode)
            nupd = AskIntegerAtLeast("Maximum number of updates",0);
         else
            nupd = ParseIntegerAtLeast(argc, argv, a++, 0);
         MT.SetTerminateCondition( TERM_NUPD, nupd );
      }
      break;
      case 'E': case 'e':
      {
         int norm;
         double errlev;
         char c1;
         if (interactive_mode)
            c1 = AskLetter("MaX, mean of Sums, mean of sQuares", "xXsSqQ");
         else
            c1 = ParseLetter(argc, argv, a++, "xXsSqQ");
         cerr << "Norm option = " << c1 << endl;
         switch(c1)
         {
            case 'X': case 'x': norm = NORM_MAX; break;
            case 'S': case 's': norm = NORM_MED; break;
            case 'Q': case 'q': norm = NORM_SQM; break;
         }
         if (interactive_mode)
            errlev = AskPositiveFloat("Minimum error level");
         else
            errlev = ParsePositiveFloat(argc, argv, a++);
         cerr << "Error threshold = " << errlev << endl;
         MT.SetTerminateCondition( TERM_ERR, norm, (double)errlev );
      }
      break;
   }  // end switch (c)

   //
   // Create triangulation based on the parameters
   //       

   if ( nextpt == RANDOM )
          T = new TRefRndDelaunay( &MT );
   else // nextpt == ERROR
          T = new TRefErrDelaunay( &MT );
   check( (T == NULL), "INSUFFICIENT MEMORY" );
   cerr << endl;

   //
   // Execute triangulation
   //

#ifdef MT_TRACER
       {
         static char first_line[] = "Multi-Tesselation build by calling:\n";
         char * cmd_line; /* commad line */
         int cmd_i; /* cursor on command line */
         int cmd_len; /* length of command line */
         cmd_len = argc;
         for (cmd_i=0; cmd_i<argc; cmd_i++)
         {
           cmd_len += strlen(argv[cmd_i]);
         }
         cmd_len += strlen(first_line);
         cmd_line = (char *) malloc((cmd_len+2)*sizeof(char));
         strcpy(cmd_line,first_line);
         for (cmd_i=0; cmd_i<argc; cmd_i++) 
         {
            strcat(cmd_line," ");
            strcat(cmd_line,argv[cmd_i]);
         }
         MT.MT_SetDescription(cmd_line);
         free (cmd_line);
       }
#endif

   T->BuildTriangulation( argv[1], argv[2] );
   cerr << "triangulation completed" << endl << endl;
   return 0;
}
