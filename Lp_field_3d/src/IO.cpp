#include <string>
#include <iostream>
#include <fstream>
#include "GeometryUtils.h"

using namespace std;

void readOff(const string& fileName, TriangleMesh& triMesh) {
  ifstream input(fileName.c_str());

  // read header
  string header;
  input >> header;

  unsigned int numVertices;
  unsigned int numFaces;
  unsigned int numEdges;
  input >> numVertices >> numFaces >> numEdges;


  // read vertex coordinates
  for (unsigned int vertCount = 0; vertCount < numVertices; vertCount++) {
    double x, y, z;
    input >> x >> y >> z;
    triMesh.addVertex(x, y, z);
  }

  // read face
  for (unsigned int faceCount = 0; faceCount < numFaces; faceCount++) {
    unsigned int numVert;
    unsigned int v0, v1, v2;
    input >> numVert >> v0 >> v1 >> v2;
    triMesh.addFace(v0, v1, v2);
  }

  /* 
     if (!(input.eof())) {
     cerr << "Error in readOff: EOF not reached" << endl;
     }
  */

  if (triMesh.getNumVertices() != numVertices) {
    cerr << "Error in readOff: wrong number of vertices" << endl;
  }

  if (triMesh.getNumTriangles() != numFaces) {
    cerr << "Error in readOff: wrong number of faces" << endl;
  }

  input.close();
}


void readTriMesh(const string& fileName, TriangleMesh& triMesh) {
  readOff(fileName, triMesh);
}


void createVTKFile(
    const string& outFileName,
    const Subdivision3& sub, const std::vector<Point3>& grid,
    const std::vector<double>& data) {

  ofstream out(outFileName.c_str());

  // header
  out << "# vtk DataFile Version 3.0" << endl;
  out << "vtk output" << endl;
  out << "ASCII" << endl;
  out << "DATASET STRUCTURED_GRID" << endl;
  out << "DIMENSIONS " <<
      sub.x << " " <<
      sub.y << " " <<
      sub.z << endl;
  out << "POINTS " << sub.x*sub.y*sub.z << " double" << endl;
  
  // structured grid
  std::vector<Point3>::const_iterator it;
  for (it = grid.begin(); it != grid.end(); ++it) {
    Point3 curr = *it;
    out << curr.x << " " << curr.y << " " << curr.z << endl;
  }
  out << endl;
  
  // data
  // header
  out << endl;
  out << "POINT_DATA " << sub.x*sub.y*sub.z << endl;
  out << "SCALARS Density double" << endl;
  out << "LOOKUP_TABLE default" << endl;

  // data
  std::vector<double>::const_iterator datait;
  for (datait = data.begin(); datait != data.end(); ++datait) {
    out << *datait << endl;
  }

  out << endl;

  out.close();
}
