#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <limits>


using namespace std;


// Constants
const double EPSILON = 1e-5;


// Types
struct Point2 {
  Point2() : x(0), y(0) {}
  Point2(double x, double y) : x(x), y(y) {}
  double x;
  double y;
  bool operator==(const Point2& p) const {
    return fabs(x-p.x)<=EPSILON && fabs(y-p.y)<=EPSILON;
  }
};

struct Vector2 {
  Vector2() : x(0), y(0) {}
  Vector2(double x, double y) : x(x), y(y) {}
  double x;
  double y;
};

typedef vector<Point2> VertexList;
typedef vector<VertexList> PolyLine;

typedef vector<vector<Point2> > Grid;
typedef vector<vector<double> > ValueGrid;


Point2 operator* (const Point2& lhs, double rhs) {
  return Point2(lhs.x * rhs, lhs.y * rhs);
}


Vector2 operator* (double s, const Vector2& rhs) {
  return Vector2(s*rhs.x, s*rhs.y);
}


Point2 operator+ (const Point2& p, const Vector2& v) {
  return Point2(p.x+v.x, p.y+v.y);
}


Vector2 operator- (const Point2& lhs, const Point2& rhs) {
  return Vector2(lhs.x - rhs.x, lhs.y - rhs.y);
}


Vector2 operator/ (const Vector2& lhs, double& rhs) {
  return Vector2(lhs.x/rhs, lhs.y/rhs);
}


void MakeZeroValueGrid(std::size_t m, std::size_t n, ValueGrid& u);
ValueGrid operator* (const ValueGrid& lhs, double rhs) {
  ValueGrid newgrid;
  MakeZeroValueGrid(lhs.size(), lhs[0].size(), newgrid);
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    for (std::size_t j = 0; j < lhs[i].size(); ++j) {
      newgrid[i][j] = lhs[i][j] * rhs;
    }
  }

  return newgrid;
}


ValueGrid operator* (double lhs, const ValueGrid& rhs) {
  return rhs * lhs;
}


ValueGrid operator/ (double lhs, const ValueGrid& rhs) {
  ValueGrid newgrid;
  MakeZeroValueGrid(rhs.size(), rhs[0].size(), newgrid);
  for (std::size_t i = 0; i < rhs.size(); ++i) {
    for (std::size_t j = 0; j < rhs[i].size(); ++j) {
      newgrid[i][j] = lhs / rhs[i][j];
    }
  }

  return newgrid;
}


/*
  Format of mesh is:
  number of loops
  number of vertices in first loop
  x1 y1
  x2 y2
  ...
  xn yn
  number of vertices in second loop
  x1 y1
  x2 y2
  ...
  xn yn
  ...
  It is assumed that:
  * each loop is closed
  * for each loop: vn = v1
  * the vertices are given in counter clockwise order
 */
void ReadMesh(const string& filename, PolyLine& mesh) {
  mesh.clear();
  
  ifstream f(filename.c_str());
  int number_loops;
  f >> number_loops;

  for (int loop = 0; loop < number_loops; ++loop) {
    int number_vertices;
    f >> number_vertices;

    vector<Point2> vertices;
    
    for (int v = 0; v < number_vertices; ++v) {
      double x, y;
      f >> x >> y;
      vertices.push_back(Point2(x,y));
    }
    
    mesh.push_back(vertices);
  }
}


// Make a regular grid of point coordinates.
// In grid[x][y], the first index corresponds to the x coordinate
// being fixed and the second index to the y coordinate being fixed.
// Example:
// grid[x][y0], grid[x][y1], ..., correspond to points on the grid
// with fixed x value and increasing y values
void MakeGrid(const Point2& pmin, const Point2& pmax,
              int number_steps_x, int number_steps_y,
              Grid& grid) {

  grid.clear();
  
  double delta_x, delta_y;
  delta_x = (pmax.x-pmin.x) / number_steps_x;
  delta_y = (pmax.y-pmin.y) / number_steps_y;

  for (int step_x = 0; step_x < number_steps_x; ++step_x) {
    double x = pmin.x + step_x * delta_x;
    vector<Point2> yraw;
    for (int step_y = 0; step_y < number_steps_y; ++step_y) {
      double y = pmin.y + step_y * delta_y;
      yraw.push_back(Point2(x,y));
    }
    grid.push_back(yraw);
  }
}


// Evaluate psi_1 on the grid p
double psi_1(const Point2& p, const PolyLine& v);
void psi_1(const Grid& p, const PolyLine& v, ValueGrid& psi) {
  psi.clear();
  
  for (unsigned int i = 0; i < p.size(); ++i) {
    vector<double> temp;
    for (unsigned int j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(psi_1(current, v));
    }
    psi.push_back(temp);
  }
}

double psi_1_normalized(const Point2&p, const PolyLine& v);
void psi_1_normalized(const Grid& p, const PolyLine& v, ValueGrid& psi) {
  psi.clear();
  
  for (unsigned int i = 0; i < p.size(); ++i) {
    vector<double> temp;
    for (unsigned int j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(psi_1_normalized(current, v));
    }
    psi.push_back(temp);
  }
}



// Evaluate psi_1 at a point
double phi_1(const Point2& p, const PolyLine& v);
double phi_1_HF(const Point2& p, const PolyLine& v);
double psi_1(const Point2& p, const PolyLine& v) {
  return 1.0 / phi_1(p, v);
}

double psi_1_normalized(const Point2& p, const PolyLine& v) {
  double psi1 = 1.0 / phi_1(p, v);
  double dpsi1dn = 1.0/2.0;
  double d = psi1/dpsi1dn;
  return d;
}


// Evaluate phi_1 on a grid.
// Vectorized version of phi_1.
void phi_1(const Grid& p, const PolyLine& v, ValueGrid& phi) {
  phi.clear();
  
  for (unsigned int i = 0; i < p.size(); ++i) {
    vector<double> temp;
    for (unsigned int j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(phi_1(current, v));
    }
    phi.push_back(temp);
  }
}


double ComputeNorm(const Vector2& v);
double ComputeSignedAngle(const Vector2& v1, const Vector2& v2);
double phi_1(const Point2& p, const PolyLine& v) {
  std::size_t number_loops = v.size();
  double sum = 0.0;

  for (std::size_t loop = 0; loop < number_loops; ++loop) {
    std::size_t number_vertices = v[loop].size();
    for (std::size_t vert = 0; vert < number_vertices; ++vert) {
      std::size_t vert_plus = (vert+1) % (number_vertices);
      if (p == v[loop][vert]) return std::numeric_limits<double>::infinity();
      if (p == v[loop][vert_plus]) return std::numeric_limits<double>::infinity();

      double a = ComputeNorm(v[loop][vert] - p);
      double b = ComputeNorm(v[loop][vert_plus] - p);
      
      double alpha = ComputeSignedAngle((v[loop][vert] - p)/a,
                                       (v[loop][vert_plus] - p)/b);
      
      double t = tan(alpha / 2.0);
      sum += t * (1.0/a + 1.0/b);
    }
  }

  return sum;
}


// Compute (a - 1) % n.
// Where a is an unsigned int \in [0, n-1].
std::size_t SubtractOneModuloN(std::size_t a, std::size_t n) {
  if (a == 0) {
    return n-1;
  }

  return (a - 1) % n;
}


double ComputeDet(const Vector2& v1, const Vector2& v2) {
  double det = 0.0;
  det = v1.x * v2.y - v2.x * v1.y;
  return det;
}


double ComputeScalarProduct(const Vector2& v1, const Vector2& v2) {
  double dot = 0.0;
  dot = v1.x * v2.x + v1.y * v2.y;
  return dot;
}


// Alternative approach for computing phi_1.
// Follow the pseudo-code given in the paper:
// K. Hormann and M. S. Floater, Mean value coordinates for arbitrary planar
// polygons, ACM Trans. on Graphics 25 (2006), 1424-1441
double phi_1_HF(const Point2& p, const PolyLine& v) {
  std::size_t number_loops = v.size();
  double sum = 0.0;

  for (std::size_t loop = 0; loop < number_loops; ++loop) {
    std::size_t number_vertices = v[loop].size();

    vector<Vector2> s(number_vertices);
    for (std::size_t vert = 0; vert < number_vertices; ++vert) {
      s[vert] = v[loop][vert] - p;
    }
    
    vector<double> r(number_vertices);
    vector<double> A(number_vertices);
    vector<double> D(number_vertices);
    for (std::size_t vert = 0; vert < number_vertices; ++vert) {
      std::size_t vert_plus = (vert + 1) % (number_vertices); 
      r[vert] = ComputeNorm(s[vert]);
      A[vert] = ComputeDet(s[vert], s[vert_plus]) / 2.0;
      D[vert] = ComputeScalarProduct(s[vert], s[vert_plus]);
    }
        
    for (std::size_t vert = 0; vert < number_vertices; ++vert) {
      std::size_t vert_plus = (vert+1) % (number_vertices);

      // unsigned int vert_minus = (vert-1) % (number_vertices);
      std::size_t vert_minus = SubtractOneModuloN(vert, number_vertices);

      // This looks strange.
      // I think that (A[vert]==0.0 || A[vert_minus]==0.0) should be
      // replaced by (A[vert]==0.0 && A[vert_minus]==0.0)
      if (std::fabs(r[vert])<=EPSILON ||
          std::fabs(A[vert])<=EPSILON ||
          std::fabs(A[vert_minus])<=EPSILON) {
        return numeric_limits<double>::max();
      } else {
        sum += (r[vert_minus] - (D[vert_minus]/r[vert])) / A[vert_minus];
        sum += (r[vert_plus] - D[vert]/r[vert]) / A[vert];
      }
    }
  }

  return sum/2.0;
}


double psi_3(const Point2& p, const PolyLine& v);
void psi_3(const Grid& p, const PolyLine& v, ValueGrid& psi) {
  psi.clear();

  for(std::size_t i = 0; i < p.size(); ++i) {
    vector<double> temp;
    for (std::size_t j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(psi_3(current, v));
    }
    psi.push_back(temp);
  }
}


double psi_3_normalized(const Point2& p, const PolyLine& v);
void psi_3_normalized(const Grid& p, const PolyLine& v, ValueGrid& psi) {
  psi.clear();

  for(std::size_t i = 0; i < p.size(); ++i) {
    vector<double> temp;
    for (std::size_t j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(psi_3_normalized(current, v));
    }
    psi.push_back(temp);
  }
}


double phi_3(const Point2& p, const PolyLine& v);
double sign(double arg);
double psi_3(const Point2& p, const PolyLine& v) {
  double phi3 = phi_3(p, v);
  double psi3 = sign(phi3) * std::pow(std::fabs(phi3), (1.0/3.0));
  psi3 = 1.0 / psi3;
  return psi3;
}


double psi_3_normalized(const Point2& p, const PolyLine& v) {
  double psi3 = psi_3(p, v);
  double dpsi3dn = (1.0 / (4.0/3.0));//^(1.0/3.0);
  dpsi3dn = std::pow(dpsi3dn, (1.0/3.0));
  return psi3/dpsi3dn;
}


double phi_3(const Point2& p, const PolyLine& v) {
  std::size_t number_loops = v.size();
  double sum = 0.0;

  for (std::size_t loop = 0; loop < number_loops; ++loop) {
    std::size_t number_vertices = v[loop].size();
    for (std::size_t vert = 0; vert < number_vertices; ++vert) {
      std::size_t vert_plus = (vert+1) % (number_vertices);

      if (p == v[loop][vert]) return std::numeric_limits<double>::infinity();
      if (p == v[loop][vert_plus]) return std::numeric_limits<double>::infinity();

      double a = ComputeNorm(v[loop][vert] - p);
      double b = ComputeNorm(v[loop][vert_plus] - p);
      double alpha = ComputeSignedAngle((v[loop][vert] - p)/a,
                                       (v[loop][vert_plus] - p)/b);
      double t = tan(alpha / 2.0);
      sum += 1.0 / 3.0 * t * (1.0/std::pow(a,3) + 1.0/std::pow(b,3))
          + 1.0/6.0 * t * (1.0 + std::pow(t,2)) * std::pow((1.0/a + 1.0/b),3);
    }
  }

  return sum;
}


// Evaluate phi_3 on a grid.
// Vectorized version of phi_3.
void phi_3(const Grid& p, const PolyLine& v, ValueGrid& phi) {
  phi.clear();
  
  for (std::size_t i = 0; i < p.size(); ++i) {
    vector<double> temp;
    for (std::size_t j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(phi_3(current, v));
    }
    phi.push_back(temp);
  }
}


double phi_5(const Point2& p, const PolyLine& v) {
  std::size_t number_loops = v.size();
  double sum = 0.0;

  for (std::size_t loop = 0; loop < number_loops; ++loop) {
    std::size_t number_vertices = v[loop].size();
    for (std::size_t vert = 0; vert < number_vertices; ++vert) {
      std::size_t vert_plus = (vert+1) % (number_vertices);

      if (p == v[loop][vert]) return std::numeric_limits<double>::infinity();
      if (p == v[loop][vert_plus]) return std::numeric_limits<double>::infinity();

      double a = ComputeNorm(v[loop][vert] - p);
      double b = ComputeNorm(v[loop][vert_plus] - p);
      double alpha = ComputeSignedAngle((v[loop][vert] - p)/a,
                                       (v[loop][vert_plus] - p)/b);
      double t = tan(alpha / 2.0);

      double a5 = std::pow(a, 5);
      double b5 = std::pow(b, 5);
      double a4 = std::pow(a, 4);
      double b4 = std::pow(b, 4);
      
      sum += 1.0/5.0 * (1.0/a5 + 1.0/b5) * t +
          ((3.0*a5 + 5.0*a4*b + 5.0*a*b4 + 3.0*b5) *
           t * (1.0 + std::pow(t,2))) / (30.0*a5*b5) +
          1.0/30.0 * (std::pow((1.0/a + 1.0/b),5))*t*std::pow((1.0+std::pow(t,2)),2);
    }
  }

  return sum;
}


// Evaluate phi_5 on a grid.
void phi_5(const Grid& p, const PolyLine& v, ValueGrid& phi) {
  phi.clear();
  
  for (std::size_t i = 0; i < p.size(); ++i) {
    vector<double> temp;
    for (std::size_t j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(phi_5(current, v));
    }
    phi.push_back(temp);
  }
}


double psi_5(const Point2& p, const PolyLine& v) {
  double phi5 = phi_5(p, v);
  double psi5 = sign(phi5) * std::pow(std::fabs(phi5), (1.0/5.0));
  psi5 = 1.0 / psi5;
  return psi5;
}


void psi_5(const Grid& p, const PolyLine& v, ValueGrid& psi) {
  psi.clear();

  for(std::size_t i = 0; i < p.size(); ++i) {
    vector<double> temp;
    for (std::size_t j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(psi_5(current, v));
    }
    psi.push_back(temp);
  }
}


double psi_5_normalized(const Point2& p, const PolyLine& v) {
  double psi5 = psi_5(p, v);
  double dpsi5dn = (1.0 / (16.0/15.0));
  dpsi5dn = std::pow(dpsi5dn, (1.0/5.0));
  return psi5/dpsi5dn;
}


void psi_5_normalized(const Grid& p, const PolyLine& v, ValueGrid& psi) {
  psi.clear();

  for(std::size_t i = 0; i < p.size(); ++i) {
    vector<double> temp;
    for (std::size_t j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(psi_5_normalized(current, v));
    }
    psi.push_back(temp);
  }
}


double squared_norm(const Vector2& v) {
  return (v.x*v.x) + (v.y*v.y);
}


double euclidean_dist(const Point2& p1, const Point2& p2) {
  return std::sqrt((p2.x-p1.x)*(p2.x-p1.x) +
                   (p2.y-p1.y)*(p2.y-p1.y));
}


double dist_to_segment(const Point2& p, const Point2& v1, const Point2& v2) {
  double param = ComputeScalarProduct(p-v1, v2-v1);
  param = param / squared_norm(v2-v1);
  double dist = 0.0;
  if (param < 0.0) {
    dist = euclidean_dist(v1, p);
  } else if (param > 1.0) {
    dist = euclidean_dist(v2, p);
  } else {
    dist = euclidean_dist(p, v1 + param*(v2 - v1));
  }
  return dist;
}


double compute_dist(const Point2& p, const PolyLine& v) {
  std::size_t number_loops = v.size();
  double dist = 99999.9;

  for (std::size_t loop = 0; loop < number_loops; ++loop) {
    std::size_t number_vertices = v[loop].size();
    for (std::size_t vert = 0; vert < number_vertices; ++vert) {
      std::size_t vert_plus = (vert+1) % (number_vertices);

      Point2 v1 = v[loop][vert];
      Point2 v2 = v[loop][vert_plus];

      double curr_dist = dist_to_segment(p, v1, v2);
      dist = std::min(dist, curr_dist);      
    }
  }

  return sign(psi_1(p, v)) * dist;
}


// Exact dist
void compute_dist(const Grid& p, const PolyLine& v, ValueGrid& dist) {
  dist.clear();

  for(std::size_t i = 0; i < p.size(); ++i) {
    vector<double> temp;
    for (std::size_t j = 0; j < p[i].size(); ++j) {
      Point2 current = p[i][j];
      temp.push_back(compute_dist(current, v));
    }
    dist.push_back(temp);
  }
}



double sign(double arg) {
  return arg > 0.0 ? 1.0 : -1.0;
}


double ComputeNorm(const Vector2& v) {
  return sqrt(v.x*v.x + v.y*v.y);
}


double ComputeSignedAngle(const Vector2& v1, const Vector2& v2) {
  double angle = atan2(v1.x*v2.y - v1.y*v2.x, v1.x*v2.x + v1.y*v2.y);
  return angle;
}


// Export the grid of values to a text file with the format:
// v11 v12 ... v1n
// v21 v22 ... v2n
// .....
// It can be read with Mathematica by the command ReadList[]
void ExportValueGrid(const string& filename, const ValueGrid& vgrid) {
  ofstream of(filename.c_str());

  for (std::size_t i = 0; i < vgrid.size(); ++i) {
    for (std::size_t j = 0; j < vgrid[i].size() - 1; ++j) {
      of << vgrid[i][j] << " ";
    }
    of << vgrid[i][vgrid.size()-1] << endl;
  }
}


void ExportGridX(const string& filename, const Grid& pgrid);
void ExportGridY(const string& filename, const Grid& pgrid);
 
void ExportGrid(const string& filename, const Grid& pgrid) {
  string gridx_filename = filename + "X";
  string gridy_filename = filename + "Y";

  ExportGridX(gridx_filename, pgrid);
  ExportGridY(gridy_filename, pgrid);
}


void ExportGridX(const string& filename, const Grid& pgrid) {
    ofstream of(filename.c_str());

    for (std::size_t i = 0; i < pgrid.size(); ++i) {
        for (std::size_t j = 0; j < pgrid[i].size()-1; ++j) {
            of << pgrid[i][j].x << " ";
        }
        of << pgrid[i][pgrid[i].size()-1].x << endl;
    }
}


void ExportGridY(const string& filename, const Grid& pgrid) {
    ofstream of(filename.c_str());

    for (std::size_t i = 0; i < pgrid.size(); ++i) {
        for (std::size_t j = 0; j < pgrid[i].size()-1; ++j) {
            of << pgrid[i][j].y << " ";
        }
        of << pgrid[i][pgrid[i].size()-1].y << endl;
    }
}


void ComputeBoundingBox(const PolyLine& v, Point2& pmin, Point2& pmax) {
  // for each loop
  //  for each vertex
  //    save the min and max coordinates in pmin and pmax

  double xmin = numeric_limits<double>::max();
  double xmax = -numeric_limits<double>::max();
  double ymin = numeric_limits<double>::max();
  double ymax = -numeric_limits<double>::max();
  
  for (std::size_t loop = 0; loop < v.size(); ++loop) {
    for (std::size_t vert = 0; vert < v[loop].size(); ++vert) {
      xmin = min(v[loop][vert].x, xmin);
      ymin = min(v[loop][vert].y, ymin);

      xmax = max(v[loop][vert].x, xmax);
      ymax = max(v[loop][vert].y, ymax);
    }
  }

  pmin.x = xmin;
  pmin.y = ymin;
  pmax.x = xmax;
  pmax.y = ymax;
}


// Enlarge the bounding box given by (pmin,pmax) by
// 2.0*scale*diagonal where diagonal is the bounding box
// diagonal.
void EnlargeBoundingBox(Point2& pmin, Point2& pmax, double scale) {
  double diagonal = ComputeNorm(pmax - pmin);
  // Enlarge the bounding box:
  pmin.x = pmin.x - scale*diagonal;
  pmax.x = pmax.x + scale*diagonal;
  pmin.y = pmin.y - scale*diagonal;
  pmax.y = pmax.y + scale*diagonal;
}


// Initialize the value grid u with 0
void MakeZeroValueGrid(std::size_t m, std::size_t n, ValueGrid& u) {
  u.clear();
  
  vector<double> temp(n);
  for (std::size_t i = 0; i < m; ++i) {
    u.push_back(temp);
  }
}


// Discrete 5-point Laplacian in Cartesian coordinates
void ComputeDiscreteLaplacian(const ValueGrid& u,
                              double hx, double hy,
                              ValueGrid& laplacian) {
  //Initialize grids u, uxx and uyy to 0
  ValueGrid uxx;
  ValueGrid uyy;
  MakeZeroValueGrid(u.size(),u[0].size(),uxx);
  MakeZeroValueGrid(u.size(),u[0].size(),uyy);
  MakeZeroValueGrid(u.size(),u[0].size(),laplacian);
  
  // General case
  for (std::size_t x = 1; x < u.size() - 1; ++x) {
    for (std::size_t y = 0; y < u[x].size(); ++y) {
      uxx[x][y] = (u[x+1][y] - 2.0*u[x][y] + u[x-1][y])
          / std::pow(hx,2);
    }
  }

  for (std::size_t x = 0; x < u.size(); ++x) {
    for (std::size_t y = 1; y < u[x].size() - 1; ++y) {
      uyy[x][y] = (u[x][y+1] - 2.0*u[x][y] + u[x][y-1])
          / std::pow(hy,2);
    }
  }

  // Left boundary
  for (std::size_t y = 0; y < u[0].size(); ++y) {
    uxx[0][y] = (2.0*u[0][y] - 4.0*u[1][y] + 2.0*u[2][y])
        / (2.0*std::pow(hx,2));
  }

  // Right boundary
  std::size_t end = u.size()-1;
  for (std::size_t y = 0; y < u[0].size(); ++y) {
    uxx[end][y] = (2.0*u[end][y] - 4.0*u[end-1][y] + 2.0*u[end-2][y])
        / (2.0*std::pow(hx,2));
  }

  // Lower boundary
  for (std::size_t x = 0; x < u.size(); ++x) {
    uyy[x][0] = (2.0*u[x][0] - 4.0*u[x][1] + 2.0*u[x][2])
        / (2.0*std::pow(hy,2));
  }

  // Upper boundary
  std::size_t yend = u[0].size()-1;
  for (std::size_t x = 0; x < u.size(); ++x) {
    uyy[x][yend] = (2.0*u[x][yend] - 4.0*u[x][yend-1] + 2.0*u[x][yend-2])
        / (2.0*std::pow(hy,2));
  }

  // Fill the laplacian:
  // laplacian u = uxx + uyy
  for (std::size_t i = 0; i < u.size(); ++i) {
    for (std::size_t j = 0; j < u[0].size(); ++j) {
      laplacian[i][j] = uxx[i][j] + uyy[i][j];
    }
  }
}


// Approximate \phi_{k+2} given: k, \phi_k.
// We have the relation:
// \Delta \phi_k = k (k+n) \phi_{k+2}
// where n is the dimension.
// \Delta is approximated by the discrete 5 point Laplacian.
void ApproximatePhiByLaplacian(
    const ValueGrid& phik,
    unsigned int k,
    const Grid& coordinates,
    const PolyLine& v,
    ValueGrid& approximate_phik2)
{
  ValueGrid laplacian_phik;
  
  double hx = coordinates[1][0].x - coordinates[0][0].x;
  double hy = coordinates[0][1].y - coordinates[0][0].y;

  ComputeDiscreteLaplacian(phik, hx, hy, laplacian_phik);
  approximate_phik2 = laplacian_phik * (1.0/(k*(k+2.0)));
}

    
void ApproximatePsi3ByLaplacian(
    const Grid& coordinates,
    const PolyLine& v,
    ValueGrid& approximate_psi)
{
  ValueGrid phi1;
  phi_1(coordinates, v, phi1);

  ValueGrid approximate_phi3;
  ApproximatePhiByLaplacian(phi1, 1, coordinates, v, approximate_phi3);

  for (std::size_t i=0; i<approximate_phi3.size(); ++i) {
    for (std::size_t j=0; j<approximate_phi3[i].size(); ++j) {
      double phi3 = approximate_phi3[i][j];
      approximate_phi3[i][j] = sign(phi3) * std::pow(std::fabs(phi3), (1.0/3.0));
    }
  }

  approximate_psi = 1.0 / approximate_phi3;  
}


void ApproximatePsi5ByLaplacian(
    const Grid& coordinates,
    const PolyLine& v,
    ValueGrid& approximate_psi)
{
  ValueGrid phi1;
  phi_1(coordinates, v, phi1);

  ValueGrid approximate_phi3;
  ApproximatePhiByLaplacian(phi1, 1, coordinates, v, approximate_phi3);

  ValueGrid approximate_phi5;
  ApproximatePhiByLaplacian(approximate_phi3,
                            3, coordinates, v,
                            approximate_phi5);

  // Compute \psi_5
  for (std::size_t i=0; i<approximate_phi5.size(); ++i) {
    for (std::size_t j=0; j<approximate_phi5[i].size(); ++j) {
      double phi5 = approximate_phi5[i][j];
      approximate_phi5[i][j] = sign(phi5) * std::pow(std::fabs(phi5), (1.0/5.0));
    }
  }

  approximate_psi = 1.0 / approximate_phi5;
}


// Additional IO routines
template<class T>
std::vector<T>
flatten(const std::vector<std::vector<T> >& grid) {
  std::vector<T> flat_grid;
  for (std::size_t i = 0; i < grid.size(); ++i) {
    for (std::size_t j = 0; j < grid[i].size(); ++j) {
        flat_grid.push_back(grid[i][j]);
    }
  }

  return flat_grid;
}


void create_vtk_file(
    const std::string& out_filename,
    const Grid& coordinates, 
    const ValueGrid& values,
    const std::string& field_title="Density")
{
  ofstream out(out_filename.c_str());

  std::size_t subx = values.size();
  std::size_t suby = values[0].size();
  std::size_t subz = 1;

  // flatten the grids:
  std::vector<Point2> grid = flatten(coordinates);
  std::vector<double> data = flatten(values); 
  
  // header
  out << "# vtk DataFile Version 3.0" << endl;
  out << "vtk output" << endl;
  out << "ASCII" << endl;
  out << "DATASET STRUCTURED_GRID" << endl;
  out << "DIMENSIONS " <<
      subx << " " <<
      suby << " " <<
      subz << endl;
  out << "POINTS " << subx*suby*subz << " double" << endl;
  
  // structured grid
  std::vector<Point2>::const_iterator it;
  for (it = grid.begin(); it != grid.end(); ++it) {
    Point2 curr = *it;
    out << curr.x << " " << curr.y << " 0.0" << endl;
  }
  out << endl;
  
  // data
  // header
  out << endl;
  out << "POINT_DATA " << subx*suby*subz << endl;
  out << "SCALARS " << field_title.c_str() << " double" << endl;
  out << "LOOKUP_TABLE default" << endl;

  // data
  std::vector<double>::const_iterator datait;
  for (datait = data.begin(); datait != data.end(); ++datait) {
    out << *datait << endl;
  }

  out << endl;

  out.close();
}


void DisplayUsage(const string& prog_name) {
  cout << "Usage: " << endl;
  cout << prog_name.c_str() << " resx resy input.polyline output.base.name" << endl;
  cout << "where: " << endl;
	cout << "\b resx and resy: grid resolution in x respectively y ";
	cout << "direction" << endl;
  cout << "\b input.polyline: 2d mesh defined by a list ";
  cout << "of vertices given in ccw order." << endl;
  cout << "\b output.base.name: base name to be used for ";
  cout << "the generated files" << endl;
  cout << endl;
}


#ifdef TEST
int main(void) {
  PolyLine v;
  ReadMesh("./riderr.vert", v);
  Point2 p(-3.6, 0.0);
  cout << psi_1(p, v) << endl;
  cout << psi_1_normalized(p, v) << endl;
  cout << psi_3(p, v) << endl;
  cout << psi_3_normalized(p, v) << endl;
  cout << psi_5(p, v) << endl;
  cout << psi_5_normalized(p, v) << endl;
  return 0;
}
#else
int main(int argc, char** argv) {
  if (argc!=5) {
    DisplayUsage(argv[0]);
    exit(-1);
  }

  string polyline_filename = argv[3];
  PolyLine v;
  ReadMesh(polyline_filename, v);

  Grid coordinates;
  Point2 pmin, pmax;
  ComputeBoundingBox(v, pmin, pmax);
  EnlargeBoundingBox(pmin, pmax, 0.1);

  // Number of steps in each direction for the grid
  int number_steps_x = atoi(argv[1]);
  int number_steps_y = atoi(argv[2]);
  MakeGrid(pmin, pmax,
           number_steps_x, number_steps_y, coordinates);

  string output_base_name = argv[4];


  // Save psi1, normalized psi1
  ValueGrid psi1;
  psi_1(coordinates, v, psi1);
  create_vtk_file(output_base_name+"_psi1.vtk", coordinates, psi1, "$\\psi_1$");

  ValueGrid psi1_normalized;
  psi_1_normalized(coordinates, v, psi1_normalized);
  create_vtk_file(output_base_name+"_psi1_normalized.vtk", coordinates,
                  psi1_normalized, "$\\Psi_1$");
  

  // Save psi3, normalized psi3
  ValueGrid psi3;
  psi_3(coordinates, v, psi3);
  create_vtk_file(output_base_name+"_psi3.vtk", coordinates, psi3, "$\\psi_3$");

  ValueGrid psi3_normalized;
  psi_3_normalized(coordinates, v, psi3_normalized);
  create_vtk_file(output_base_name+"_psi3_normalized.vtk", coordinates,
                  psi3_normalized, "$\\Psi_3$");


  // Save psi5, normalized psi5
  ValueGrid psi5;
  psi_5(coordinates, v, psi5);
  create_vtk_file(output_base_name+"_psi5.vtk", coordinates, psi5, "$\\psi_5$");

  ValueGrid psi5_normalized;
  psi_5_normalized(coordinates, v, psi5_normalized);
  create_vtk_file(output_base_name+"_psi5_normalized.vtk", coordinates,
                  psi5_normalized, "$\\Psi_5$");

    
  // Save exact dist
  ValueGrid dist;
  compute_dist(coordinates, v, dist);
  create_vtk_file(output_base_name+"_dist.vtk", coordinates, dist, "dist");
  

  return EXIT_SUCCESS;
}
#endif
