===============================
Some sample terrain data:
===============================

ustica.pts   Points, accepted input for RefDel

ustica.seg   Points + constraint edges, accepted input for RefCTD

ustica.tri   Triangulation, accepted input for DecDel, SiDecDel

ustica.cdt   Constrained triangulazion, accepted input for DecCDT, SiDecCDT


===============================
Exaples:
===============================

Output is an (unconstrained) Delaunary triangulation:

programs/RefDel ustica.pts  output.tri

programs/DecDel ustica.tri output.tri
programs/SiDecDel ustica.tri output.tri

Output is a constrained Delaunary triangulation:

programs/RefDel ustica.seg output.cdt

programs/RefCTD ustica.cdt  output.cdt
programs/SiRefCTD ustica.cdt  output.cdt

===============================
Data format:
===============================


POINTS (extension .pts):

N
x1 y1 z1
...
xN yN zN

where N=number of points, (xI yI,zI)=coordinates of the I-th point.

POINTS AND SEGMENTS (extension .seg):

N
x1 y1 z1
...
xN yN zN
M 
a1 b1
...
aM bM

where N=number of points, (xI yI,zI)=coordinates of the I-th point,
M=number of segments, (aI,bI)=indexes of the endpoints of the I-th 
segment, where indexes range from 0 to N-1.

TRIANGULATION (extension .tri):

N
x1 y1 z1
...
xN yN zN
M 
a1 b1 c1
...
aM bM cM

where N=number of points, (xI yI,zI)=coordinates of the I-th point,
M=number of triangles, (aI,bI,cI)=indexes of the three vertices of
the I-th triangle, where indexes range from 0 to N-1.

CONSTRAINED TRIANGULATION (extension .cdt):

N
x1 y1 z1
...
xN yN zN
M 
a1 b1 c1
...
aM bM cM
K
d1 e1
...
dK eK

where N=number of points, (xI yI,zI)=coordinates of the I-th point,
M=number of triangles, (aI,bI,cI)=indexes of the three vertices of
the I-th triangle, K=number of constraint edges, (dI,eI)=indexes of
the endpoints of the I-th edge, where indexes range from 0 to N-1.
