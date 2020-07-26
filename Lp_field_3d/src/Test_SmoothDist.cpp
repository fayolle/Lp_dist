#include <cassert>
#include <iostream>
#include <vector>

#include "GeometryUtils.h"
#include "IO.h"

#include "Distance.h"


using namespace std;

double computePhi(const Point3& pt, int k, const Triangle& tri);

// 
// Test smooth distance to some triangles at some points
// 
int main(int argc, char** argv) {
	Point3 p1(1,-1,1);
	Point3 p2(1,1,1);
	Point3 p3(-1,1,1);
	
	Triangle t(p1,p2,p3);
	
	Point3 eval(-2,-2,-2);
	
	computePhi(eval, 3, t);
	
  return 0;
}
