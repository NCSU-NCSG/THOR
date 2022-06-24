// -----------------------------------------------------------------------------
//
//  Godiva Quarter Sphere
//
//  Constructive Solid Geometry, OpenCASCADE geometry kernel
//
// -----------------------------------------------------------------------------

// Instead of constructing a model in a bottom-up fashion with Gmsh's built-in
// geometry kernel, starting with version 3 Gmsh allows you to directly use
// alternative geometry kernels. Here we use the OpenCASCADE kernel:

SetFactory("OpenCASCADE");

// Godiva is just a sphere of radius 8.7407, also make a boxes to remove left half, bottom, and back
Sphere(1) = {0,0,0,8.7407};
Box(2) = {-10,-10,-10, 10,20,20};
Box(3) = {0,-10,-10, 10,20,10};
Box(4) = {0,-10,0, 10,10,10};

// this creates a new volume removing the box from the sphere.
// this also deletes the previous volumes so we don't get overlaying meshes.
BooleanDifference(5) = { Volume{1}; Delete; }{ Volume{2:4}; Delete; };

// Assign a mesh size to all the points of all the volumes:
MeshSize{ PointsOf{ Volume{:}; } } = 0.1;