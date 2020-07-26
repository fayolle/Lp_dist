#include <string>
#include "GeometryUtils.h"

#ifndef IO_H
#define IO_H

void readTriMesh(const std::string& fileName, TriangleMesh& triMesh);
void createVTKFile(
    const std::string& outFileName,
    const Subdivision3& sub, const std::vector<Point3>& grid,
    const std::vector<double>& data);

#endif
