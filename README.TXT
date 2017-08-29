===============
doc
===============
Documentation. Entry file is index.html

===============
sample_data
===============
Some sample data in the formats needed and produced by programs.
For more details, see file readme.txt inside.

===============
programs
===============
Sources and makefiles for Delaunay refinement and decimation,
constrained Delaunay refinement and decimation.
They are compiled with:
  cd programs
  make -f simpmake RefDel
  make -f simpmake DecDel
  make -f simpmake SiDecDel
  make -f simpmake_cdt RefCDT
  make -f simpmake_cdt DecCDT
  make -f simpmake_cdt SiDecCDT 2> logg
for version which does not contruct a multi-tesselation, and with same
commands and makefiles mtmake and mtmake_cdt (instead of simpmake and
simpmake_cdt) for version which contructs a multi-tesselation.
More instructions within directory programs.

Note: compilation gives many warnings for deprecated conversion from
type string to char*, but this does not affect the correctness of the
programs.

===============
gui
===============
Version with graphical user interface, using libraries OpenGL and Glut.

===============
util
===============
Some useful utility programs. See file readme.txt inside.


