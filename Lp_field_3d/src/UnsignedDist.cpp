#include <cassert>
#include <iostream>
#include <vector>

#include "GeometryUtils.h"
#include "IO.h"

#include "Distance.h"


using namespace std;


void printHelp() {
  cout << "UnsignedDist: " << endl;
  cout << "\tCompute the unsigned distance to a surface" << endl;
  cout << "\tapproximated by a triangle mesh." << endl;

  cout << endl;

  cout << "Syntax:" << endl;
  cout << "\t./UnsignedDist inputmesh outputvtk xm ym zm xM yM zM nx ny nz squared?" << endl;
  cout << "where:" << endl;
  cout << "\t inputmesh: input triangle mesh in OFF format," << endl;
  cout << "\t outputvtk: output scalar field in VTK format," << endl;
  cout << "\t xm ym zm: coordinates of the lower corner of the object's bounding box," << endl;
  cout << "\t xM yM zM: coordinates of the upper corner of the object's bounding box," << endl;
  cout << "\t nx ny nz: grid resolution," << endl;
	cout << "\t squared?: 0/1 whether the distance is squared or not" << endl;	
}


// 
// Compute the exact distance for a given mesh on a grid
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
  std::vector<Point3> structuredGrid;
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

	bool is_squared = true;
	if (atoi(argv[12]) == 0) is_squared = false;
	
  // evaluate exact distance
  std::vector<double> results = computeUnsignedDist(structuredGrid, triMesh, is_squared);

  string outFileName = argv[2];
  createVTKFile(outFileName, subdivisions, structuredGrid, results);
  
  return 0;
}

