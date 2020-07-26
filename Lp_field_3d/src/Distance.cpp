#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <list>

#include <omp.h>

#include "GeometryUtils.h"
#include "LinearAlgebra.h"
#include "TriangleIntegral.h"


const double EPSILON = 1e-7;
const double PI = 3.14159265358979;


// functions for robust predicates
double exactinit();
double orient2d(double *pa, double *pb, double *pc);
double orient3d(double* pa, double* pb, double* pc, double* pd);


// alternative function for computing the potential phi (follows code from M. Carley)
double potential_phi(const Point3& pt, int k, const Triangle& tri);


// 
// Smooth distance computations
//

// Compute the rotation matrix and translation along z that transforms
// the coordinate system such that the input triangle "tri" after
// transformation is in the plane z=0.
// This is an alternative approach that does not use the SVD.
static void getTransformation(const Triangle& tri, Matrix33& rot, double& ztrans) {
  Point3 v0 = tri.v0;
  Point3 v1 = tri.v1;
  Point3 v2 = tri.v2;

  Vector3 u0 = v1 - v0;
  Vector3 e0 = normalize(u0); 
  Vector3 u1 = v2 - v0;
  Vector3 e1 = u1 - projection(u1, u0);
  e1 = normalize(e1);

  Vector3 n = tri.getUnitNormal();

  Matrix33 rotationMatrix;
  // rotationMatrix = [vec column e0, vec column e1, vec column n]
  rotationMatrix(0,0) = e0.x;
  rotationMatrix(0,1) = e1.x;
  rotationMatrix(0,2) = n.x;

  rotationMatrix(1,0) = e0.y;
  rotationMatrix(1,1) = e1.y;
  rotationMatrix(1,2) = n.y;

  rotationMatrix(2,0) = e0.z;
  rotationMatrix(2,1) = e1.z;
  rotationMatrix(2,2) = n.z;

  // The matrix that we got is the transposed of the sought matrix
  rot = transpose(rotationMatrix);

  // get the translation amount along ez
  Point3 transformedV0 = rot * v0;
  ztrans = transformedV0.z;
}


// Given the rotation matrix and translation along ez, compute the transformed
// coordinates of the input point inputPt
static Point3 applyTransformation(const Point3& inputPt, const Matrix33& rot, double ztrans) {
  Point3 result = rot * inputPt;
  result.z = result.z - ztrans;
  return result;
}


// eval I^{p}_{m,n}[alpha,theta] with m=0 and n=-p
// (see Eq. 14 in Carley's paper)
static double evalI(int p, double alpha, double theta) {
  assert(alpha >= 0);
  assert(alpha <= 1.0);
  
  /*
    if (alpha < 0) { 
    std::cerr << "alpha is negative, ";
    std::cerr << alpha << std::endl; 
    }
  
    if (alpha > 1.0) {
    std::cerr << "alpha is greater than 1, ";
    std::cerr << alpha << std::endl;
    }
  */

  assert(p == -2 || p == -3 || p == -4 || p == -5 || p == -6);
	
  if (p == -2) {
    // used in phi_1
    double alphaSq = alpha * alpha;
    double alphaPrime = sqrt(1.0 - alphaSq);
    return (theta - alphaPrime*atan(alphaPrime*tan(theta))) / alphaSq;

  } else if (p == -3) {
    // used in phi_2
    // I_{0,3}^{-3} in Table 3 of Carley's paper
    double alphaSq = alpha * alpha;
    double alphaPrime = sqrt(1.0 - alphaSq);
    double alphaPrimeSq = alphaPrime * alphaPrime;
    double delta = sqrt(1.0 - alphaSq * sin(theta) * sin(theta));
    double alphaCube = alphaSq * alpha;
    
    /* return (-alphaPrimeSq / alphaSq)*(sin(theta) / delta) +
       1.0 / alphaCube * asin(alpha * sin(theta)); */
    
    // I_{0,3}^{-3} from my own calculation with Mathematica:
    return (-alphaPrimeSq/alphaSq)*(sin(theta)/delta) + 
        atan(alpha*sin(theta)/delta) / alphaCube;

  } else if (p == -4) {
    // used in phi_3
    double alphaSq = alpha * alpha;
    double alphaPrime = sqrt(1.0 - alphaSq);
    double delta = sqrt(1.0 - alphaSq * sin(theta) * sin(theta));
    double alphaFour = alphaSq * alphaSq;

    /* return (theta
       - (1.0 + 0.5*alphaSq)*alphaPrime*atan(alphaPrime*tan(theta))
       + (alphaSq*(-0.5 + 0.5*alphaSq)*sin(2.0*theta))/(2.0 - alphaSq + alphaSq*cos(2.0*theta))
       ) / (alphaSq*alphaSq); */

    return (2.0*theta + ((-2 + alphaSq + alphaFour)*(atan(alphaPrime*tan(theta))))/(alphaPrime) + (alphaSq * (alphaSq - 1) * sin(2.0*theta))/(2.0 * delta * delta)) / (2.0 * alphaFour);
    
  } else if (p == -5) {
    // used in phi_4
    double alphaSq = alpha*alpha;
    //double alphaFive = alphaSq * alphaSq * alpha;
    double delta = sqrt(1.0 - alphaSq * sin(theta) * sin(theta));
    double alphaFour = alphaSq * alphaSq;
    double alphaFive = alphaFour * alpha;
    double alphaPrimeSq = 1.0 - alphaSq;
    double deltaSq = delta * delta;
    double deltaCube = deltaSq * delta;
    
    /* double a = atan((sqrt(2.0)*alpha*sin(theta))/sqrt(alphaSq*cos(2.0*theta)-alphaSq+2.0)) / alphaFive;
    
       double b = (2.0*sqrt(2.0)*(alphaSq-1.0)*sin(theta)*(-alphaSq*alphaSq+(alphaSq+2.0)*alphaSq*cos(2.0*theta)+alphaSq+3.0))
       / (3.0*alphaSq*alphaSq*pow((alphaSq*cos(2.0*theta)-alphaSq+2),1.5));
       return a + b; */

    double a = atan(alpha*sin(theta)/delta) / alphaFive;

    double b = ((-alphaPrimeSq)*
                (2.0*(2.0+alphaSq)*deltaSq - alphaPrimeSq)*sin(theta))
        / (3.0 * alphaFour * deltaCube);
    
    return a + b;

  } else if (p == -6) {
    // used in phi_5
    double alphaSq = alpha*alpha;
    double a = (4.0*(alphaSq-1)*(alphaSq-1)*alphaSq*sin(2.0*theta))
        / pow((alphaSq*cos(2.0*theta)-alphaSq+2.0),2.0);
    double alphaFour = alphaSq * alphaSq;
    double b = (3.0*(alphaFour+alphaSq-2.0)*alphaSq*sin(2.0*theta))
        / (alphaSq*cos(2.0*theta)-alphaSq+2.0);
    double alphaSix = alphaFour * alphaSq;
    double c = ((3.0*alphaSix+alphaFour+4.0*alphaSq-8.0)*atan(sqrt(1.0-alphaSq)*tan(theta)))
        / (sqrt(1.0-alphaSq));

    return (1.0/(8.0*alphaSix))*(a + b + c + 8.0*theta);

  } else {
    std::cerr << "only defined for the following values of p:"
              << " -2, -3, -4, -5, -6"
              << std::endl;
    return 0.0;
  }
}


// Follow Carley, Section 2
// Using robust predicates from J. Shewchuck
static double
subTriangleContribution(double z, int k, const Triangle& tri) {
  // compute the orientation (signed area) of the triangle
  double pa[2] = {tri.v0.x, tri.v0.y};
  double pb[2] = {tri.v1.x, tri.v1.y};
  double pc[2] = {tri.v2.x, tri.v2.y};
  double orientation = orient2d(pa, pb, pc);

  // if the orientation is not zero, calculate a, r1, r2, theta, phi taking into account the orientation (using equations 4, 8 and 33)
  if (fabs(orientation) <= EPSILON) {
    return 0.0;
  }

  // r1 = || p1 - p0 || (See Fig. 1)
  Vector3 v0v1 = tri.v1 - tri.v0;
  double r1 = v0v1.norm();
	
  // r2 = || p2 - p0 ||
  Vector3 v0v2 = tri.v2 - tri.v0;
  double r2 = v0v2.norm();
	
	
  if (r1 <= EPSILON || r2 <= EPSILON) {
    return 0.0;
  }
  
  // Eq. 33
  double theta = atan2(orientation, dot(v0v1,v0v2));

  
  // Eq. 4
  double costheta = cos(theta);
  double sintheta = sin(theta);
  if (fabs(sintheta) <= EPSILON) {
    return 0.0;
  }

  double a = (r2 * costheta - r1) / (r2 * sintheta);

  
  // Eq. 8 (right side)
  double phi = atan(a);

  
  // Eq. 12
  double zSq = z * z;
  double betaSq = (r1*r1)/(1.0+a*a) + zSq;
  double beta = std::sqrt(betaSq);
  double alpha = std::fabs(z) / beta;
  
  // evaluate the required integral using the appropriate formula
  int gamma = -(k + 3);
  double Ib = evalI(gamma+2, alpha, theta+phi);
  double Ia = evalI(gamma+2, alpha, phi);

  double I = (1.0 / (gamma+2)) *
      (std::pow(beta, gamma+2) * (Ib - Ia) -
       std::pow(fabs(z), gamma+2) * theta);
  
  return I;
}


// Field phik at point pt from one triangle
double computePhi(const Point3& pt, int k, const Triangle& tri) {
  // Transform the coordinate system such that the triangle lies in the plane z = 0
  // Compute the transformation (rotation and z-translation)
  Matrix33 rot;
  double ztrans = 0.0;
  getTransformation(tri, rot, ztrans);

  // Apply the transformation to pt and each vertex of tri
  Point3 transformedPt = applyTransformation(pt, rot, ztrans);

  Point3 v0;
  tri.getVertex0(v0);
  Point3 transformedV0 = applyTransformation(v0, rot, ztrans);
  assert(fabs(transformedV0.z) <= EPSILON);

  Point3 v1;
  tri.getVertex1(v1);
  Point3 transformedV1 = applyTransformation(v1, rot, ztrans);
  assert(fabs(transformedV1.z) <= EPSILON);

  Point3 v2;
  tri.getVertex2(v2);
  Point3 transformedV2 = applyTransformation(v2, rot, ztrans);
  assert(fabs(transformedV2.z) <= EPSILON);

  // if the field point lie in the support plane of the triangle
  // then the contribution from this triangle is 0
  double pa[3] = {pt.x, pt.y, pt.z};
  double pb[3] = {v0.x, v0.y, v0.z};
  double pc[3] = {v1.x, v1.y, v1.z};
  double pd[3] = {v2.x, v2.y, v2.z};
  double orient = orient3d(pa, pb, pc, pd);
  bool inPlane = (fabs(orient) <= EPSILON);
  if (inPlane) return 0.0;
	
  
  // projection of transformedPt on the plane supporting tri
  Point3 projected(transformedPt.x, transformedPt.y, 0.0);

  // Decompose the triangle into 3 subtriangles
  // (use notation from Carley's paper)
  Point3 p0 = projected;
  Point3 p1 = transformedV0;
  Point3 p2 = transformedV1;
  Point3 p3 = transformedV2;

  Triangle t1(p0, p1, p2);
  Triangle t2(p0, p2, p3);
  Triangle t3(p0, p3, p1);

  // Sum the contribution of each sub-triangle
  double result = 0.0;
  result += subTriangleContribution(transformedPt.z, k, t1);
  result += subTriangleContribution(transformedPt.z, k, t2);
  result += subTriangleContribution(transformedPt.z, k, t3);
  
  
  return -transformedPt.z * result;
}


static double computePhi(const Point3& pt, int k, const TriangleMesh& triMesh) {
  double sum = 0.0;
  // compute phi_k for each triangle
  std::size_t numTriangles = triMesh.getNumTriangles();
  for(std::size_t triIdx = 0; triIdx < numTriangles; ++triIdx) {
    Triangle tri;
    triMesh.getTriangle(triIdx, tri);

    sum += computePhi(pt, k, tri);
    //sum += potential_phi(pt, k, tri);  
  }     

  return sum;
}


// Determine if a point is inside a triangle; no dependencies (CGAL or other)
static bool pointInTriangle(const Point3& pt, const Triangle& tri) {
  double pa[3] = {pt.x, pt.y, pt.z};
  double pb[3] = {tri.v0.x, tri.v0.y, tri.v0.z};
  double pc[3] = {tri.v1.x, tri.v1.y, tri.v1.z};
  double pd[3] = {tri.v2.x, tri.v2.y, tri.v2.z};

  double orient = orient3d(pa, pb, pc, pd);
  bool inPlane = (fabs(orient) <= EPSILON);
  if (!inPlane) return false;

  // pt is in the support plane of triangle
  // Check if it is inside the triangle using
  // its barycentric coordinates
  Vector3 u = tri.v1 - tri.v0;
  Vector3 v = tri.v2 - tri.v0;
  Vector3 w = pt - tri.v0;
  Vector3 vCrossW = cross(v, w);
  Vector3 vCrossU = cross(v, u);

  if (dot(vCrossW, vCrossU) < 0.0) {
    return false;
  }

  Vector3 uCrossW = cross(u, w);
  Vector3 uCrossV = cross(u, v);

  if (dot(uCrossW, uCrossV) < 0) {
    return false;
  }

  double denom = uCrossV.norm();
  double r = vCrossW.norm() / denom;
  double t = uCrossW.norm() / denom;

  return (r <= 1 &&
          t <= 1 &&
          r + t <= 1);
}


static double
computePsi(const Point3& pt, int k, const TriangleMesh& tri_mesh) {

  // compute phi_k
  double phik = 0.0;
  std::size_t numTriangles = tri_mesh.getNumTriangles();
  for (std::size_t triIdx = 0; triIdx < numTriangles; ++triIdx) {
    Triangle tri;
    tri_mesh.getTriangle(triIdx, tri);

    // point located on the surface: distance is 0
    if (pointInTriangle(pt, tri)) {
      return 0.0; 
    }

    phik += computePhi(pt, k, tri);

    // This is an alternative computation for phi inspired from Carley.
    // The results are not better. 
    // phik += potential_phi(pt, k, tri);
  }

  // compute (1/phi_k)^(1/k)
  double sign = 1.0;
  if (phik < 0.0) sign = -1.0; 
  double inverse = 1.0 / fabs(phik);
  double exponent = 1.0 / k;
  double psik = pow(inverse, exponent);

  return sign * psik;
}


std::vector<double>
computePsi(const std::vector<Point3>& pts, int k, const TriangleMesh& triMesh) {
  // init needed by robust predicates
  exactinit();
  std::size_t gridsize = pts.size();
  std::vector<double> results(gridsize);
#pragma omp parallel for 
  for (int i = 0; i < static_cast<int>(gridsize); ++i) {
    results[i] = computePsi(pts[i], k, triMesh);
  }
  return results;
}


//
// Compute normalization weights for psi
// to have a unit normal derivative.
// k is the order (in psi_k)
// Assume n (dimension of space) = 3
//
static double computeNormalizationWeight(int k) {
  if (k == 0) {
    return 2.0 * PI;
  }

  if (k == 1) {
    return 2.0;
  }

  return ((k-2)+1.0)/((k-2)+3.0) * computeNormalizationWeight(k-2);
}


//
// Normalization of psi_k
//
static std::vector<double> normalize(const std::vector<double>& psi, int k) {
  // compute the weight to get a unit normal derivative
  double normal_weight_k = computeNormalizationWeight(k);
  double normal_weight = pow(normal_weight_k, 1.0/double(k));
  
  // apply it to each value psi in the input vector
  std::size_t gridsize = psi.size();
  std::vector<double> normalized(gridsize);
#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(gridsize); ++i) {
    normalized[i] = (normal_weight * psi[i]);
  }
  return normalized;
}


static void normalize(std::vector<double>& psi, int k) {
  // Same as above but modify the input vector psi
  double normal_weight_k = computeNormalizationWeight(k);
  double normal_weight = pow(normal_weight_k, 1.0/double(k));
  std::size_t gridsize = psi.size();
#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(gridsize); ++i) {
    psi[i] = normal_weight * psi[i];
  }
}


// Normalized version of psi_k
std::vector<double>
computeNormalizedPsi(const std::vector<Point3>& pts,
                     int k, const TriangleMesh& triMesh) {
  std::vector<double> psi = computePsi(pts, k, triMesh);
  normalize(psi, k);
  return psi;
}


//
// Compute psi_1 by using the approach of Ju et al.
// see pseudocode in Fig. 4 of:
// Mean value coordinates for closed triangular meshes
// Ju, Schaefer, Warren
//
static double computePsi1(const Point3& pt, const TriangleMesh& triMesh) {
  double phi1 = 0.0;
  std::size_t numTriangles = triMesh.getNumTriangles();
  for (std::size_t triIdx = 0; triIdx < numTriangles; ++triIdx) {
    Triangle tri;
    triMesh.getTriangle(triIdx, tri);

    // point located on the surface: distance is 0
    if (pointInTriangle(pt, tri)) {
      return 0.0; 
    }

    Vector3 uv0 = normalize(tri.v0 - pt);
    Vector3 uv1 = normalize(tri.v1 - pt);
    Vector3 uv2 = normalize(tri.v2 - pt);

    double l[3];
    l[0] = norm(uv1 - uv2);
    l[1] = norm(uv2 - uv0);
    l[2] = norm(uv0 - uv1);

    double theta[3];
    theta[0] = 2.0 * asin(l[0]/2.0);
    theta[1] = 2.0 * asin(l[1]/2.0);
    theta[2] = 2.0 * asin(l[2]/2.0);

    double h = (theta[0] + theta[1] + theta[2]) / 2.0;

    double c[3];
    c[0] = 2.0*sin(h)*sin(h-theta[0]) / (sin(theta[1])*sin(theta[2])) - 1.0;
    c[1] = 2.0*sin(h)*sin(h-theta[1]) / (sin(theta[2])*sin(theta[0])) - 1.0;
    c[2] = 2.0*sin(h)*sin(h-theta[2]) / (sin(theta[0])*sin(theta[1])) - 1.0;

    double ss = signDet(uv0, uv1, uv2);

    double s[3];
    s[0] = ss * sqrt(1.0 - std::min(c[0]*c[0],1.0));
    s[1] = ss * sqrt(1.0 - std::min(c[1]*c[1],1.0));
    s[2] = ss * sqrt(1.0 - std::min(c[2]*c[2],1.0));

    double psi1eps = 0.0001;
    if (fabs(s[0]) < psi1eps || fabs(s[1]) < psi1eps || fabs(s[2]) < psi1eps)
      continue;

    double dv0 = dist(pt, tri.v0);
    double dv1 = dist(pt, tri.v1);
    double dv2 = dist(pt, tri.v2);

    // Note:
    // in the pseudocode Fig. 4, a factor 1/2 is missing compared to Eq. 10
    double w[3];
    w[0] = (theta[0]-c[1]*theta[2]-c[2]*theta[1])
                                             / (2.0*dv0*sin(theta[1])*s[2]);
    w[1] = (theta[1]-c[0]*theta[2]-c[2]*theta[0])
                                             / (2.0*dv1*sin(theta[2])*s[0]);
    w[2] = (theta[2]-c[1]*theta[0]-c[0]*theta[1])
                                             / (2.0*dv2*sin(theta[0])*s[1]);

    phi1 += w[0] + w[1] + w[2];
  }

  return 1.0 / phi1;
}


std::vector<double> 
computePsi1(const std::vector<Point3>& pts, 
	    const TriangleMesh& triMesh) {
  exactinit();
  std::size_t gridsize = pts.size();
  std::vector<double> results(gridsize);
#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(gridsize); ++i) {
    results[i] = computePsi1(pts[i], triMesh);
  }
  return results;
}


std::vector<double>
computeNormalizedPsi1(const std::vector<Point3>& pts,
                      const TriangleMesh& triMesh) {
  std::vector<double> psi = computePsi1(pts, triMesh);
  normalize(psi, 1);
  return psi;
}


//
// Compute psi_k using numerical integration
//
static double computePhiNumerical(const Point3& pt, int order,
                                  const NumSamples& num_samples,
                                  const Triangle& tri) {
  Vector3 triangle_normal = tri.getUnitNormal();
  
  TriangleIntegral integral(num_samples, tri.v0, tri.v1, tri.v2);
  std::vector<double> fi;
	
  for (unsigned int s = 0; s < integral.get_num_samples(); ++s) {
    Point3 y = integral.get_sample(s);
    Vector3 xy = y - pt;
    double xydotn = dot(xy, triangle_normal);
    double r = norm(xy);
    fi.push_back(xydotn / pow(r, double(order)+3.0));
  }
  
  double phi = integral.integrate(fi);
  return phi;
}


static double computePhiNumericalFast27(const Point3& pt, int k,
                                        const Triangle& tri) {
  Vector3 triangleNormal = tri.getUnitNormal();
  Point3 samples[27];
  double weights[27];
  initArrays27(tri.v0, tri.v1, tri.v2, samples, weights );

  // the scalar product between xy and n should not change on
  // a flat triangle
  Point3 y = samples[0];
  Vector3 xy = y - pt;
  double xydotn = dot(xy, triangleNormal);
  double r = norm(xy);

  double phi = 0.0;
  phi += weights[0]*(xydotn / pow(r, double(k)+3.0));

  // other samples
  for (unsigned int s = 1; s < 27; ++s) {
    y = samples[s];
    xy = y - pt;
    r = norm(xy);
    phi += weights[s] * (xydotn / pow(r, double(k)+3.0));
  }

  return phi;
}


static double computePhiNumericalFast(const Point3& pt, int k,
                                      const Triangle& tri) {
  Vector3 triangleNormal = tri.getUnitNormal();
  Point3 samples[32];
  double weights[32];
  initArrays32(tri.v0, tri.v1, tri.v2, samples, weights );
  
  // the scalar product between xy and n should not change on
  // a flat triangle
  Point3 y = samples[0];
  Vector3 xy = y - pt;
  double xydotn = dot(xy, triangleNormal);
  double r = norm(xy);
  
  double phi = 0.0;
  phi += weights[0]*(xydotn / pow(r, double(k)+3.0));
  
  // other samples
  for (unsigned int s = 1; s < 32; ++s) {
    y = samples[s];
    xy = y - pt;
    r = norm(xy);
    phi += weights[s] * (xydotn / pow(r, double(k)+3.0));
  }
  
  return phi;
}


static double computePsiNumerical(const Point3& pt, int k,
                                  const NumSamples& num_samples,
                                  const TriangleMesh& triMesh) {
  // compute phi_k
  double phik = 0.0;
  std::size_t numTriangles = triMesh.getNumTriangles();
  for (std::size_t triIdx = 0; triIdx < numTriangles; ++triIdx) {
    Triangle tri;
    triMesh.getTriangle(triIdx, tri);
    
    // point located on the surface: distance is 0
    if (pointInTriangle(pt, tri)) {
      return 0.0; 
    }
    
    phik += computePhiNumerical(pt, k, num_samples, tri);
  }
  
  // compute (1/phi_k)^(1/k)
  double sign = 1.0;
  if (phik < 0.0) sign = -1.0; 
  double inverse = 1.0 / fabs(phik);
  double exponent = 1.0 / k;
  double psik = pow(inverse, exponent);
  
  return sign * psik;
}


static double computePsiNumericalFast(const Point3& pt, int k,
                                      const TriangleMesh& triMesh) {
  // compute phi_k
  double phik = 0.0;
  std::size_t numTriangles = triMesh.getNumTriangles();
  for (std::size_t triIdx = 0; triIdx < numTriangles; ++triIdx) {
    Triangle tri;
    triMesh.getTriangle(triIdx, tri);
    
    // point located on the surface: distance is 0
    if (pointInTriangle(pt, tri)) {
      return 0.0; 
    }
    
    phik += computePhiNumericalFast(pt, k, tri);
  }
  
  // compute (1/phi_k)^(1/k)
  double sign = 1.0;
  if (phik < 0.0) sign = -1.0; 
  double inverse = 1.0 / fabs(phik);
  double exponent = 1.0 / k;
  double psik = pow(inverse, exponent);
  
  return sign * psik;
}


std::vector<double>
computePsiNumerical(const std::vector<Point3>& pts, int order,
                    const NumSamples& num_samples, const TriangleMesh& triMesh) {
  exactinit();
  std::size_t gridsize = pts.size();
  std::vector<double> results(gridsize);
#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(gridsize); ++i) {
    results[i] = computePsiNumerical(pts[i], order, num_samples, triMesh);
  }
  
  return results;
}


std::vector<double>
computePsiNumericalFast(const std::vector<Point3>& pts, int order,
                        const TriangleMesh& triMesh) {
  exactinit();
  std::size_t gridsize = pts.size();
  std::vector<double> results(gridsize);
#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(gridsize); ++i) {
    results[i] = computePsiNumericalFast(pts[i], order, triMesh);
  }
  
  return results;
}


std::vector<double>
computeNormalizedPsiNumericalFast(const std::vector<Point3>& pts, int order,
                                  const TriangleMesh& triMesh) { 
  std::vector<double> psi = computePsiNumericalFast(pts, order, triMesh);
  normalize(psi, order);
  return psi;
}


std::vector<double>
computeNormalizedPsiNumerical(const std::vector<Point3>& pts, int order,
                              const NumSamples& num_samples,
                              const TriangleMesh& triMesh) { 
  std::vector<double> psi = computePsiNumerical(pts, order, num_samples, triMesh);
  normalize(psi, order);
  return psi;
}


//
// Hybrid method
//

// switcheps: used to decide whether to use the analytical
//  expression (when close to the surface,
//  i.e. if dist to surface <= switcheps * meshsize )
//  or the numerical expression (away)
static double
computePsiHybrid(const Point3& pt, int k,
		 const TriangleMesh& triMesh, double switcheps) {
  double trimeshlen = triMesh.getDiagonalLength();

  double psik = computePsiNumericalFast(pt, k, triMesh);
  if (fabs(psik) <= switcheps*trimeshlen) {
    psik = computePsi(pt, k, triMesh);
    return psik;
  }

  return psik;
}


std::vector<double>
computePsiHybrid(const std::vector<Point3>& pts, int order,
		 const TriangleMesh& triMesh, double switcheps) {
  exactinit();
  std::size_t gridsize = pts.size();
  std::vector<double> results(gridsize);
#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(gridsize); ++i) {
    results[i] = computePsiHybrid(pts[i], order, triMesh, switcheps);
  }
  return results;
}


std::vector<double>
computeNormalizedPsiHybrid(const std::vector<Point3>& pts, int order,
                           const TriangleMesh& triMesh, double switcheps) {
  std::vector<double> psi = computePsiHybrid(pts, order, triMesh, switcheps);
  normalize(psi, order);
  return psi;
}

