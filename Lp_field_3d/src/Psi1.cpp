#include <cassert>
#include <iostream>
#include <vector>

#include "GeometryUtils.h"
#include "IO.h"

#include "Distance.h"


using namespace std;


void printHelp() {
  cout << "psi1: " << endl;
  cout << "\tCompute a first order approximate smooth distance field to a surface" << endl;
  cout << "\tapproximated by a triangle mesh." << endl;
  cout << "\tThe approach follows: Mean Value coordinates by Ju, Schaeffer and Warren." << endl;

  cout << endl;
  
  cout << "Syntax:" << endl;
  cout << "\t./psi1 inputmesh outputvtk xm ym zm xM yM zM nx ny nz normalized?" << endl;
  cout << "where:" << endl;
  cout << "\t inputmesh: input triangle mesh in OFF format," << endl;
  cout << "\t outputvtk: output scalar field in VTK format," << endl;
  cout << "\t xm ym zm: coordinates of the lower corner of the object's bounding box," << endl;
  cout << "\t xM yM zM: coordinates of the upper corner of the object's bounding box," << endl;
  cout << "\t nx ny nz: grid resolution," << endl;
  cout << "\t normalized?: normalize the field (0: no, 1: yes)." << endl;
}


// 
// Compute psi1 for a given mesh on a grid
// 
int main(int argc, char** argv) {
  if (argc != 13) {
    printHelp();
    return -1;
  }

  // read mesh file
  TriangleMesh triMesh;
  string fileName = argv[1];
  readTriMesh(fileName, triMesh);

  // Create a structured grid for evaluating psi
  vector<Point3> structuredGrid;
  double lcx = atof(argv[3]);
  double lcy = atof(argv[4]);
  double lcz = atof(argv[5]);

  double ucx = atof(argv[6]);
  double ucy = atof(argv[7]);
  double ucz = atof(argv[8]);
  
  Point3 leftCorner(lcx,lcy,lcz);
  Point3 rightCorner(ucx,ucy,ucz);

  unsigned subx = atoi(argv[9]);
  unsigned suby = atoi(argv[10]);
  unsigned subz = atoi(argv[11]);
  Subdivision3 subdivisions(subx, suby, subz);
  
  createGrid(leftCorner, rightCorner, subdivisions, structuredGrid);

  bool isNormalized = false;
  int normalized = atoi(argv[12]);
  if (normalized == 1) {
    isNormalized = true;
  }
  
  
  // evaluate psi1 on the grid
  vector<double> results;
  if (isNormalized) {
    cout << "Compute normalized psi1" << endl;
    results = computeNormalizedPsi1(structuredGrid, triMesh);
  } else {
    cout << "Compute (non-normalized) psi1" << endl;
    results = computePsi1(structuredGrid, triMesh);
  }

  string outFileName = argv[2];
  createVTKFile(outFileName, subdivisions, structuredGrid, results);
  
  return 0;
}
