#include <cassert>
#include <iostream>
#include <vector>

#include "GeometryUtils.h"
#include "IO.h"

#include "Distance.h"
#include "TriangleIntegral.h"


using namespace std;


void printHelp() {
  cout << "smoothdisthybrid: " << endl;
  cout << "\tCompute an approximate smooth distance field to a surface" << endl;
  cout << "\tapproximated by a triangle mesh. The potential integral is" << endl;
  cout << "\tcomputed using an hybrid method: an analytical expression" << endl;
  cout << "\tclose to the surface and a quadrature rule far." << endl; 

  cout << endl;

  cout << "Syntax:" << endl;
  cout << "\t./smoothdisthybrid inputmesh outputvtk xm ym zm xM yM zM nx ny nz order normalized? switcheps" << endl;
  cout << "where:" << endl;
  cout << "\t inputmesh: input triangle mesh in OFF format," << endl;
  cout << "\t outputvtk: output scalar field in VTK format," << endl;
  cout << "\t xm ym zm: coordinates of the lower corner of the object's bounding box," << endl;
  cout << "\t xM yM zM: coordinates of the upper corner of the object's bounding box," << endl;
  cout << "\t nx ny nz: grid resolution," << endl;
  cout << "\t order: approximation order (currently available: 1, 2, 3, 4, 5), " << endl;
  cout << "\t normalized?: normalize the field (0: no, 1: yes)." << endl;
  cout << "\t switcheps: parameter used for switching between numerical and analytical" << endl;
}


// 
// Compute the smooth distance for a given mesh on a grid
// 
int main(int argc, char** argv) {
  if (!(argc == 14 || argc == 15)) {
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

  int order = atoi(argv[12]);
  cout << "order: " << order << endl;

  bool isNormalized = false;
  int normalized = atoi(argv[13]);
  if (normalized == 1) {
    isNormalized = true;
  }

  double switcheps = 0.01;
  if (argc == 15) {
    switcheps = atof(argv[14]);
  }


  // evaluate psi on the grid
  vector<double> results;
  if (isNormalized) {
    cout << "Use normalized psi" << endl;
    results =
        computeNormalizedPsiHybrid(structuredGrid, order, triMesh, switcheps);
  } else {
    cout << "Use non-normalized psi" << endl;
    results = computePsiHybrid(structuredGrid, order, triMesh, switcheps);
  }

  string outFileName = argv[2];
  createVTKFile(outFileName, subdivisions, structuredGrid, results);
  
  return 0;
}
