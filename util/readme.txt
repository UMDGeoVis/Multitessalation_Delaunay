========================
PROGRAMS
========================

gengrid.c
	Generate a grid of random elevations and write them to a file.

grid2tri.c
	Convert a grid file into a triangulation file with the right
	syntax for program Delaunaymt.
	Divide the grid into rectangles and divide each rectangle into
	two triangles by drawing one of its diagonals.

grid2pts.c
	Convert a grid file into a point file with the right
	syntax for program Delaunaymt.

All programs are written in C and compiled with: 
gcc PROG.c -o PROG

========================
FORMAT FOR GRID FILE
========================

Text file containing:

  num_row  
  num_column
  z_1 
  z_2
  ...
  z_K

where:
	num_row (int) is the number of rows in the grid,
	num_column (int) is the number of columns in the grid,
	z_I (float) is the height of the I-th point, 
	K (int) is equal to num_row*num_column.
