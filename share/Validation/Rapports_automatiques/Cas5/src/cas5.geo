lc = 0.02;
Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {2,0.0,0.0,lc};
Point(3) = {2,1,0.0,lc};
Point(4) = {0.0,1,0.0,lc};
Point(5) = {0.0,-0.2,0.0,lc};
Point(6) = {2.0,-0.2,0.0,lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line(5) = {1,5};
Line(6) = {5,6};
Line(7) = {6,2};

Line Loop(8) = {2,3,4,1};
Line(9) = {2,1};
Line Loop(10) = { 5,6,7,9 };
Plane Surface(6) = {8};
Plane Surface(7) = {10};

Coherence;
Physical Line("Bas") = {6};
Physical Line("Droit") = {2};
Physical Line("Haut") = {3};
Physical Line("Gauche") = {4};
Physical Line("Cote") = {5,7};
Physical Surface("dom") = {6,7};

