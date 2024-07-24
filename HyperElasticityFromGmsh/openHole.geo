SetFactory("OpenCASCADE");                                                                               

lc = 1e-1;

xC = 0.5;
yC = 0.5;
r =  0.3;


x0 = 0;
y0 = 0;
L  = 1;

// Define square with hole

Point(1) = {x0, y0, 0, lc};
Point(2) = {x0+L, y0, 0, lc};
Point(3) = {x0+L, y0+L, 0, lc};
Point(4) = {x0, y0+L, 0, lc};

Point(5) = {x0, yC+r, 0, lc};
Point(6) = {x0+r, yC, 0, lc};
Point(7) = {x0, yC-r, 0, lc};
Point(8) = {x0, yC, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};

//Define circle


Point(9) = {xC, yC, 0, lc};
Point(10) = {xC+r, yC, 0, lc};
Point(11) = {xC, yC+r, 0, lc};
Point(12) = {xC-r, yC, 0, lc};
Point(13) = {xC, yC-r, 0, lc};

Circle(8) = {10, 9, 11};
Circle(9) = {11, 9, 12};
Circle(10) = {12, 9, 13};
Circle(11) = {13, 9, 10};

Curve Loop(2) = {8, 9, 10, 11};


Plane Surface(1) = {1};
Plane Surface(2) = {2};

BooleanDifference{Surface{1}; Delete;}{ Surface{2}; Delete;}

Delete { Curve{:}; }
Delete { Point{:}; }

Physical Curve("BotEdges", 12) = {12};
Physical Curve("TopEdges", 15) = {15};
Physical Surface(1) = {1, 2, 3, 4};


Mesh.MeshSizeMin = 0.01;
Mesh.MeshSizeMax = 0.02;

Coherence;

//Save metaMaterial.msh
