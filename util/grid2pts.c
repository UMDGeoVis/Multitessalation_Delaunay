/*****************************************************************************
File converter.
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
/*
Convert a grid file into a point file with the right
syntax for program Delaunaymt.

The format of the grid file must be the following:
  num_row  
  num_column
  z1 
  z2
  ...
  zK
where num_row and num_column are two integers specifying the number
of rows and columns in the grid, zI is the height of the I-th point,
and K = num_row*num_column.
*/

#include <stdio.h>

int main(int argc, char **argv)
{
  FILE * f_in, * f_out; /* input and output files */
  int n_row, n_col; /* number of rows and columns in the grid */
  int i,j;          /* indexes for scanning the grid */
  float delta_x;    /* length for horizontal grid edges */
  float delta_y;    /* length for vertical grid edges */
  float x,y,z;      /* coordinates of a point */

  /**** PROCESS COMMAND LINE OPTIONS ****/
  if (argc<5) 
  {
    fprintf(stderr,"Usage: %s delta_x delta_y input_file output_file\n",
            argv[0]);
    fprintf(stderr,"where delta_x, delta_y = lenght of grid edges\n");
    return 0;
  }
  if ( sscanf(argv[1], "%f",&delta_x) != 1)
  {
    fprintf(stderr,"Invalid first argument (delta_x)\n");
    return 0;
  }
  if ( sscanf(argv[2], "%f",&delta_y) != 1)
  {
    fprintf(stderr,"Invalid second argument (delta_y)\n");
    return 0;
  }
  fprintf(stderr,"delta_x = %f delta_y = %f input_file = %s output_file %s\n",
          delta_x, delta_y, argv[3], argv[4]);
  f_in = fopen(argv[3],"r");
  if (!f_in) 
  {
    fprintf(stderr,"Cannot open input file %s\n",argv[3]);
    return 0;
  }
  f_out = fopen(argv[4],"w");
  if (!f_out)
  {
    fprintf(stderr,"Cannot open output file %s\n",argv[4]);
    return 0;
  }

  /**** READ HEAD OF GRID FILE ****/
  if ( fscanf(f_in,"%d %d", &n_row, &n_col) != 2 )
  {
    fprintf(stderr,"Cannot read number of rows and columns\n");
    return 0;
  }
      
  /**** WRITE POINTS ****/
  fprintf(f_out, "%d\n", n_row*n_col);
  for (i=0; i<n_row; i++)
  for (j=0; j<n_col; j++)
  {
     if ( fscanf(f_in,"%f", &z) != 1 )
     {
        fprintf(stderr,"Cannot read height of point %d,%d\n",i,j);
        return 0;
     }
     x = i*delta_x; y = j*delta_y;
     fprintf(f_out, "%g %g %g\n",x,y,z);
  }
  return 1;
}

