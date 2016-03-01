// Gmsh project created on Mon Nov 11 16:25:22 2013


r1 = 100.0;

h1 = 0.8;

Point(1) = {0, 0, 0, h1};

Point(2) = {r1,0, 0, h1};

Point(3) = {-r1,0, 0, h1};


Circle(1) = {2, 1, 3};

Circle(2) = {3,1,2};


Line Loop(1) = {1,2};

Plane Surface(1) = {1};


