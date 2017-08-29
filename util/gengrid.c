/*
Generate a grid of random elevations and write them to a file.
The format of the grid file is:
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
#include <stdlib.h> /* per rand */

int main(void)
{
  int DimX, DimY; /* dimensioni griglia */
  float MinZ, MaxZ; /* min e max quota */
  float z;
  int i,j;
  int randvalue;
  
  fprintf(stderr,"Dimensioni griglia DimX DimY:");
  fscanf(stdin,"%d %d", &DimX, &DimY);
  fprintf(stderr,"Min e Max quota:");
  fscanf(stdin,"%f %f", &MinZ, &MaxZ);

  /* scrittura numero vertici */
  printf("%d\n%d\n", DimX,DimY);
  
  /* scrittura quote */
  for (j=0;j<DimY;j++)
    for (i=0;i<DimX;i++)
    {
      z = MinZ + ( (MaxZ-MinZ)*rand()/RAND_MAX );
      printf("%f\n", z);
    }

  return 1;
}
