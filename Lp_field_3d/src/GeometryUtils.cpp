#include "GeometryUtils.h"
#include <iostream>
using namespace std;


Vector3 operator- (const Point3& a, const Point3& b) {
  return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}


Vector3 operator/ (const Vector3& v, double d) {
  return Vector3(v.x / d, v.y / d, v.z / d);
}


double dot(const Vector3& v0, const Vector3& v1) {
  return v0.x*v1.x + v0.y*v1.y + v0.z*v1.z;
}


Vector3 cross(const Vector3& v0, const Vector3& v1) {
  double crossx = v0.y * v1.z - v0.z * v1.y;
  double crossy = v0.z * v1.x - v0.x * v1.z;
  double crossz = v0.x * v1.y - v0.y * v1.x;

  return Vector3(crossx, crossy, crossz);
}


double dist(const Point3& a, const Point3& b) {
  Vector3 ab = b - a;
  return ab.norm();
}


Point3 operator* (const Point3& pt, double s) {
  Point3 result;
  result.x = pt.x * s;
  result.y = pt.y * s;
  result.z = pt.z * s;
  return result;
}


Point3 operator* (double s, const Point3& pt) {
  return pt * s;
}


Vector3 operator* (const Vector3& v, double s) {
  Vector3 result;
  result.x = v.x * s;
  result.y = v.y * s;
  result.z = v.z * s;
  return result;
}


Vector3 operator* (double s, const Vector3& v) {
  return v * s;
}


Point3 operator+ (const Point3& p, const Vector3& v) {
  Point3 result;
  result.x = p.x + v.x;
  result.y = p.y + v.y;
  result.z = p.z + v.z;
  return result;
}


Point3 operator+ (const Vector3& v, const Point3& p) {
  return p + v;
}


Vector3 operator+ (const Vector3& v1, const Vector3& v2) {
  Vector3 result;
  result.x = v1.x + v2.x;
  result.y = v1.y + v2.y;
  result.z = v1.z + v2.z;
  return result;
}


Vector3 operator- (const Vector3& v1, const Vector3& v2) {
  Vector3 result;
  result.x = v1.x - v2.x;
  result.y = v1.y - v2.y;
  result.z = v1.z - v2.z;
  return result;
}


double norm(const Vector3& v) {
  return v.norm();
}


Point3 operator+ (const Point3& p1, const Point3& p2) {
  return Point3(p1.x + p2.x,
                p1.y + p2.y,
                p1.z + p2.z);
}


// return a unit vector with same direction as v
Vector3 normalize(const Vector3& v) {
  double norm = v.norm();
  if (norm == 0.0) {
    norm = 1.0;
  }
  return Vector3(v.x/norm, v.y/norm, v.z/norm);
}


// compute the projection of v on u
Vector3 projection(const Vector3& v, const Vector3& u) {
  double num = dot(v, u);
  double den = dot(u, u);
  if (den == 0.0) {
    cerr << "Projection on a null vector" << endl;
    den = 1.0;
  }
  double alpha = num / den;
  return Vector3(alpha*u.x, alpha*u.y, alpha*u.z);
}


double signDet(const Vector3& v0, const Vector3& v1, const Vector3& v2) {
  double det = v0.x*v1.y*v2.z + v0.y*v1.z*v2.x +
    v0.z*v1.x*v2.y - v0.z*v1.y*v2.x - 
    v0.x*v1.z*v2.y - v0.y*v1.x*v2.z;
  return (det >= 0.0) ? 1.0 : -1.0;
}


// Create a regular grid
void createGrid(const Point3& leftCorner, const Point3& rightCorner,
                const Subdivision3& subdivisions,
                std::vector<Point3>& grid) {

  double dx = (rightCorner.x - leftCorner.x) / subdivisions.x;
  double dy = (rightCorner.y - leftCorner.y) / subdivisions.y;
  double dz = (rightCorner.z - leftCorner.z) / subdivisions.z;

  double currentX = leftCorner.x;
  double currentY = leftCorner.y;
  double currentZ = leftCorner.z;
      
  for (unsigned zsub = 0; zsub < subdivisions.z; ++zsub) {
    currentY = leftCorner.y;
    for (unsigned ysub = 0; ysub < subdivisions.y; ++ysub) {
      currentX = leftCorner.x;
      for (unsigned xsub = 0; xsub < subdivisions.x; ++xsub) {
        grid.push_back(Point3(currentX, currentY, currentZ));
        currentX = currentX + dx;
      }      
      currentY = currentY + dy;
    }
    currentZ = currentZ + dz;
  }
}


Point3 operator* (const Matrix33& m, const Point3& p) {
  Point3 result(0.0,0.0,0.0);
  result.x = m(0,0) * p.x + m(0,1) * p.y + m(0,2) * p.z;
  result.y = m(1,0) * p.x + m(1,1) * p.y + m(1,2) * p.z;
  result.z = m(2,0) * p.x + m(2,1) * p.y + m(2,2) * p.z;
  return result;
}
