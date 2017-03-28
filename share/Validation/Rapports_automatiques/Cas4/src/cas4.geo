lc = 0.02;
Point(1) = {0.0,0.0,0.0,lc};
Point(2) = {2,0.0,0.0,lc};
Point(3) = {2,1,0.0,lc};
Point(4) = {0.0,1,0.0,lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {2,3,4,1};
Plane Surface(6) = {5};
Coherence;
Physical Line("Bas") = {1};
Physical Line("Droit") = {2};
Physical Line("Haut") = {3};
Physical Line("Gauche") = {4};
Physical Surface("dom") = {6};

