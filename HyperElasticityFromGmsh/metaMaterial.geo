SetFactory("OpenCASCADE");                                                                               

lc = 1e-1;

xC = 0.65;
yC = 0.65;
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
Line(4) = {4, 5};
Circle(5) = {5, 8, 6};
Circle(6) = {6, 8, 7};
Line(7) = {7,1};
//Line(5) = {5, 8};
//Line(6) = {8, 7};
//Line(7) = {7, 1};


//Curve Loop(1) = {1, 2, 3, 4,5,6,7};
//Curve Loop(1) = {1, 2, 3, 4,5,6,7};
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7};

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

Symmetry {1, 0, 0, -1} {
  Duplicata { Surface{1}; }
}

Symmetry {0, 1, 0, -1} {
  Duplicata { Surface{1}; }
}

Symmetry {0, 1, 0, -1} {
  Duplicata { Surface{2}; }
}
Physical Curve("HorEdges", 11) = {1, 12, 22, 32};
//Physical Curve("VerEdges", 12) = {23,27,7,2};//13,17,36,33};
Physical Curve("VerEdges", 12) = {23,27,7,2,13,17,36,33};
Physical Surface(1) = {1, 2, 3, 4};

Delete { Curve{:}; }
Delete { Point{:}; }
//Delete { Surface{:}; }
//Delete { Volume{:}; Surface{dels()}; Curve{:}; Point{:}; }

Mesh.MeshSizeMin = 0.01;
Mesh.MeshSizeMax = 0.02;

Coherence;

//Save metaMaterial.msh
