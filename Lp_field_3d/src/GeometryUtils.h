#ifndef GEOMETRY_UTILS
#define GEOMETRY_UTILS


#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "LinearAlgebra.h"


// Utility type 
class Triple {
public:
  Triple(unsigned int v0, unsigned int v1, unsigned int v2) : 
    v0(v0), v1(v1), v2(v2) {}

  unsigned int getIndex0() { return v0; }
  unsigned int getIndex1() { return v1; }
  unsigned int getIndex2() { return v2; }


  unsigned int v0;
  unsigned int v1;
  unsigned int v2;
};


// Geometric type
class Point3 {
public:
  Point3() : x(0.0), y(0.0), z(0.0) {}
  Point3(double x, double y, double z) : x(x), y(y), z(z) {}

  double getX() { return x; }
  double getY() { return y; }
  double getZ() { return z; }


  double x;
  double y;
  double z;
};


class Vector3 {
public:
  Vector3() : x(0.0), y(0.0), z(0.0) {}
  Vector3(double x, double y, double z) : x(x), y(y), z(z) {}

  double getX() { return x; }
  double getY() { return y; }
  double getZ() { return z; }
  void setX(double v) { x = v; }
  void setY(double v) { y = v; }
  void setZ(double v) { z = v; }

  double norm() const {
    return sqrt(x*x + y*y + z*z);
  }

  double norm2() const {
    return x*x + y*y + z*z;
  }
  
  double x;
  double y;
  double z;
};


// Operations
Vector3 operator- (const Point3& a, const Point3& b);
Vector3 operator/ (const Vector3& v, double d);
double dot(const Vector3& v0, const Vector3& v1);
Vector3 cross(const Vector3& v0, const Vector3& v1);
double dist(const Point3& a, const Point3& b);
Point3 operator* (const Point3& pt, double scalar);
Point3 operator* (double scalar, const Point3& pt);
Vector3 operator* (const Vector3& v, double scalar);
Vector3 operator* (double scalar, const Vector3& v);
Point3 operator+ (const Point3& p, const Vector3& v);
Point3 operator+ (const Vector3& v, const Point3& p);
Vector3 operator+ (const Vector3& v1, const Vector3& v2);
Vector3 operator- (const Vector3& v1, const Vector3& v2);
double norm(const Vector3& v);

Point3 operator+ (const Point3& p1, const Point3& p2);

// return a unit vector with same direction as v
Vector3 normalize(const Vector3& v);

// compute the projection of v on u
Vector3 projection(const Vector3& v, const Vector3& u);

// return the sign of the det of the matrix made of v0,v1,v2
double signDet(const Vector3& v0, const Vector3& v1, const Vector3& v2);


// Geometric types
class Triangle {
public:
  Triangle() : v0(0.0,0.0,0.0), v1(0.0,0.0,0.0), v2(0.0,0.0,0.0) {}

  Triangle(const Point3& v0, const Point3& v1, const Point3& v2) :
    v0(v0), v1(v1), v2(v2) {}
    
  void setVertices(const Point3& v0, const Point3& v1, const Point3& v2) {
    this->v0 = v0;
    this->v1 = v1;
    this->v2 = v2;
  }

  void getVertex0(Point3& v) const {
    v = v0;
  }

  void getVertex1(Point3& v) const {
    v = v1;
  }

  void getVertex2(Point3& v) const {
    v = v2;
  }

  Vector3 getUnitNormal() const {
    Vector3 v01 = v1 - v0;
    Vector3 v02 = v2 - v0;

    Vector3 normal = cross(v01, v02);
    double norm = normal.norm();
    if (norm != 0.0) normal = normal / norm;

    return normal;
  }

  Point3 v0;
  Point3 v1;
  Point3 v2;
};


class TriangleMesh {
public:
  void addVertex(double x, double y, double z) {
    vertices.push_back(Point3(x,y,z));
  }

  void addFace(unsigned int v0, unsigned int v1, unsigned int v2) {
    triangles.push_back(Triple(v0,v1,v2));
  }

  std::size_t getNumVertices() const {
    return vertices.size();
  }

  std::size_t getNumTriangles() const {
    return triangles.size();
  }

  void getTriangle(std::size_t i, Point3& v0, Point3& v1, Point3& v2) const {
    Triple tri = triangles[i];
    unsigned int idx0 = tri.getIndex0();
    v0 = vertices[idx0];
    unsigned int idx1 = tri.getIndex1();
    v1 = vertices[idx1];
    unsigned int idx2 = tri.getIndex2();
    v2 = vertices[idx2];
  }

  void getTriangle(std::size_t i, Triangle& tri) const {
    Point3 v0(0.0,0.0,0.0);
    Point3 v1(0.0,0.0,0.0);
    Point3 v2(0.0,0.0,0.0);
    getTriangle(i, v0, v1, v2);
    tri.setVertices(v0, v1, v2);
  }

  struct CmpX {
    bool operator() (const Point3& a, const Point3& b) {
      return a.x < b.x;
    }
  };
  
  struct CmpY {
    bool operator() (const Point3& a, const Point3& b) {
      return a.y < b.y;
    }
  };
  
  struct CmpZ {
    bool operator() (const Point3& a, const Point3& b) {
      return a.z < b.z;
    }
  };
  
  void getBoundingBox(Point3& leftCorner, Point3& rightCorner) const {
    std::vector<Point3>::iterator minx =
        min_element(vertices.begin(), vertices.end(), CmpX());
    std::vector<Point3>::iterator miny =
        min_element(vertices.begin(), vertices.end(), CmpY());
    std::vector<Point3>::iterator minz =
        min_element(vertices.begin(), vertices.end(), CmpZ());

    std::vector<Point3>::iterator maxx =
        max_element(vertices.begin(), vertices.end(), CmpX());
    std::vector<Point3>::iterator maxy =
        max_element(vertices.begin(), vertices.end(), CmpY());
    std::vector<Point3>::iterator maxz =
        max_element(vertices.begin(), vertices.end(), CmpZ());

    leftCorner.x = (*minx).x;
    leftCorner.y = (*minx).y;
    leftCorner.z = (*minx).z;

    rightCorner.x = (*maxx).x;
    rightCorner.y = (*maxx).y;
    rightCorner.z = (*maxx).z;  
  }

  double getDiagonalLength() const {
    Point3 lowerleft;
    Point3 upperright;
    getBoundingBox(lowerleft, upperright);
    return dist(lowerleft, upperright);
  }

  
private:
  mutable std::vector<Triple> triangles;
  mutable std::vector<Point3> vertices;
};


// Create a regular grid
class Subdivision3 {
 public:
  Subdivision3() : x(0), y(0), z(0) {}
  Subdivision3(unsigned x, unsigned y, unsigned z) :
      x(x), y(y), z(z) {}
  
  unsigned x;
  unsigned y;
  unsigned z;
};

void createGrid(const Point3& leftCorner, const Point3& rightCorner,
                const Subdivision3& subdivisions,
                std::vector<Point3>& grid);


// Operations on Point3
Point3 operator* (const Matrix33& m, const Point3& p);



#endif
