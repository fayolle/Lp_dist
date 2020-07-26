#include <vector>
#include "GeometryUtils.h"
#include "TriangleIntegral.h"


#ifndef DISTANCE_H
#define DISTANCE_H


// Currently only supports k = 1,2,3,4,5
std::vector<double> computePsi(const std::vector<Point3>& pts, int k, const TriangleMesh& tri_mesh);
std::vector<double> computeNormalizedPsi(const std::vector<Point3>& pts, int k, const TriangleMesh& tri_mesh);


// For experiments only (?)
// We can control the number of samples for the numerical integration
// at the expense of runtime speed
std::vector<double> computePsiNumerical(const std::vector<Point3>& pts, int k, const NumSamples& num_samples, const TriangleMesh& tri_mesh);
std::vector<double> computeNormalizedPsiNumerical(const std::vector<Point3>& pts, int k, const NumSamples& num_samples, const TriangleMesh& tri_mesh);


// Faster but number of integration samples is fixed
std::vector<double> computePsiNumericalFast(const std::vector<Point3>& pts, int k,  const TriangleMesh& tri_mesh);
std::vector<double> computeNormalizedPsiNumericalFast(const std::vector<Point3>& pts, int k, const TriangleMesh& tri_mesh);


// Combination of analytical and numerical integration
std::vector<double> computePsiHybrid(const std::vector<Point3>& pts, int k, const TriangleMesh& tri_mesh, double switcheps=0.01);
std::vector<double> computeNormalizedPsiHybrid(const std::vector<Point3>& pts, int k, const TriangleMesh& tri_mesh, double switcheps=0.01);


// Compute Psi_1 using Ju et al approach
std::vector<double> computePsi1(const std::vector<Point3>& pts, const TriangleMesh& tri_mesh);
std::vector<double> computeNormalizedPsi1(const std::vector<Point3>& pts, const TriangleMesh& tri_mesh);


// Exact (unsigned) dist
std::vector<double>
computeUnsignedDist(const std::vector<Point3>& pts, const TriangleMesh& tri_mesh, bool is_squared);


#endif

