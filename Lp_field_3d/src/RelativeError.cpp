#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include <algorithm>
#include <functional>


#include "GeometryUtils.h"
#include "IO.h"

#include "Distance.h"


using namespace std;


// TODO
void printHelp() {
  cerr << "At least 11 arguments required" << endl;
}


const double RELATIVE_ERROR_EPSILON = 1e-9;


//
// Compute the relative error between the exact distance field
// and the smooth approximate distance field
//
std::vector<double> computeRelativeError(
    const std::vector<double>& exactField,
    const std::vector<double>& smoothField) {

  assert(exactField.size() == smoothField.size());

  std::vector<double> relativeErrorField;
  for (unsigned idx = 0; idx < exactField.size(); ++idx) {
    double exact = exactField[idx];
    double smooth = smoothField[idx];

    double denominator = fabs(exact);
    if (denominator < RELATIVE_ERROR_EPSILON) {
      denominator += RELATIVE_ERROR_EPSILON;
    }

    double result = fabs(exact - smooth) / denominator;
    relativeErrorField.push_back(result);
  }

  return relativeErrorField;
}


// 
// Compute the smooth distance for a given mesh on a grid
// 
int main(int argc, char** argv) {
  if (argc != 12) {
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

  // evaluate exact distance
  std::vector<double> exact = computeUnsignedDist(structuredGrid, triMesh, false);

  // evaluate smooth distance
  std::vector<double> smooth = computePsi(structuredGrid, 2, triMesh);
  std::transform(smooth.begin(), smooth.end(), smooth.begin(),
                 std::ptr_fun<double,double>(std::fabs));
  
  std::vector<double> smoothNormalized =
      computeNormalizedPsi(structuredGrid, 2, triMesh);
  std::transform(smoothNormalized.begin(), smoothNormalized.end(),
                 smoothNormalized.begin(),
                 std::ptr_fun<double,double>(std::fabs));

  // compute relative error
  std::vector<double> relativeError = computeRelativeError(exact, smooth);
  
  string outFileName = argv[2];
  createVTKFile(outFileName, subdivisions, structuredGrid, relativeError);

  std::vector<double> relativeErrorNormalized = computeRelativeError(exact, smoothNormalized);
  string normalizedFileName = outFileName + "normalized";
  createVTKFile(normalizedFileName, subdivisions, structuredGrid, relativeErrorNormalized);
  
  return 0;
}

