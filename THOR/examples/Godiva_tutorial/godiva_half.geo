// -----------------------------------------------------------------------------
//
//  Godiva Half Sphere
//
//  Constructive Solid Geometry, OpenCASCADE geometry kernel
//
// -----------------------------------------------------------------------------

// Instead of constructing a model in a bottom-up fashion with Gmsh's built-in
// geometry kernel, starting with version 3 Gmsh allows you to directly use
// alternative geometry kernels. Here we use the OpenCASCADE kernel:

SetFactory("OpenCASCADE");

// Godiva is just a sphere of diameter 8.7407, also make a box to remove the left half
Sphere(1) = {0,0,0,4.37035};
Box(2) = {-5,-5,-5, 5,10,10};

// this creates a new volume removing the box from the sphere.
// this also deletes the previous volumes so we don't get overlaying meshes.
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

// Assign a mesh size to all the points of all the volumes:
MeshSize{ PointsOf{ Volume{:}; } } = 2;
