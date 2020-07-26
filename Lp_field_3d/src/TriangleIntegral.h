/*
  Compute integral of a function on a triangle by quadrature.
  References:
  - The finite element method: its basis and fundamentals,
  O.C. Zienkiewicz, R.L. Taylor & J.Z. Zhu, pp. 164-166
  - Nonlinear Finite Element Methods, Peter Wriggers,
  pp. 118
  - Guassian Quadrature Formulas for triangles, G.R. Cowper
  - Asymmetric cubature formulas for polynomial integration
  in the triangle and square, Journal of Computational and
  Applied Mathematics, Mark Taylor
 */

#ifndef TRIANGLEINTEGRAL_H
#define TRIANGLEINTEGRAL_H

#include <cmath>
#include <vector>
#include "GeometryUtils.h"


enum NumSamples {
  one_point = 1,
  three_points = 3,
  seven_points = 7,
  twenty_four_points = 24,
  twenty_seven_points = 27,
  thirty_two_points = 32
};


// convert from int to enum and enum to int
template<typename T>
T my_enum_convert(int);


class TriangleIntegral {
 public:
  TriangleIntegral(NumSamples num_samples, const Point3& p1, const Point3& p2, const Point3& p3);
  ~TriangleIntegral();
  Point3 get_sample(int i);
  double get_weight(int i);
  double integrate(double fi[]);
  double integrate(const std::vector<double>& fi);
  double get_area();
  unsigned int get_num_samples();
  
 private:
  double _area;
  double* _weights;
  Point3* _samples;
  unsigned int _num_samples;

  const static double _samples_weights_1[];
  const static double _samples_weights_3[];
  const static double _samples_weights_4[];
  const static double _samples_weights_7[];
  const static double _samples_weights_24[];
  const static double _samples_weights_27[];
  const static double _samples_weights_32[];
};


void initArrays27(const Point3& p1, const Point3& p2, const Point3& p3,
				  Point3 samples[], double weights[]);
void initArrays32(const Point3& p1, const Point3& p2, const Point3& p3,
				  Point3 samples[], double weights[]);

#endif
