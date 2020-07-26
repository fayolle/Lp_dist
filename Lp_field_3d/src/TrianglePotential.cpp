#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>

#include "GeometryUtils.h"
#include "LinearAlgebra.h"


const double EPSILON = 1e-7;
const double PI = 3.14159265358979;


// functions for robust predicates
double exactinit();
double orient2d(double *pa, double *pb, double *pc);
double orient3d(double* pa, double* pb, double* pc, double* pd);


// eval I^{p}_{m,n}[alpha,theta] with m=0 and n=-p
// (see Eq. 14 in Carley's paper)
static double evalI(int p, double alpha, double theta) {
  if (alpha < 0) { 
    std::cerr << "alpha is negative, ";
    std::cerr << alpha << std::endl; 
  }
  
  if (alpha > 1.0) {
    std::cerr << "alpha is greater than 1, ";
    std::cerr << alpha << std::endl;
  }

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




void subtriangle_polar_quad(double r1, double r2, double theta, 
                            double z, int l, 
                            double *I1)
{
  double a, b, b2, z2, k, kd, t1, t0, S0, S1, C0, C1;

  double phi;

  a = (r2*cos(theta) - r1)/r2/sin(theta) ; phi = atan(a) ;

  b2 = (r1*r1 + z*z*(1+a*a))/(1+a*a) ; b = sqrt(b2) ;
  k = fabs(z)/b ; kd = sqrt(1.0-k*k) ;

  assert(z != 0.0) ;
  z2 = z*z ;

  t1 = theta+(phi) ; t0 = (phi) ;
  S0 = sin(t0) ; S1 = sin(t1) ; C0 = cos(t0) ; C1 = cos(t1) ;


  int gamma = -(l + 3);
  double Ib = evalI(gamma+2, k, t1);
  double Ia = evalI(gamma+2, k, t0);
  
  *I1 = (1.0 / (gamma+2.0) * (pow(b, gamma+2) * (Ib - Ia) - pow(fabs(z), gamma + 2) * theta));
}


static void subtriangle_potential_quad(double *p,
				       double *x1, double *x2, 
				       double r1, double r2,
				       double ont,
				       int k, 
				       double *I)
{
  double J1, th;

  if ( r1 == 0.0 || r2 == 0.0 ) return;
  if ( fabs(th = ((x2[0]-p[0])*(x1[0]-p[0]) + 
		  (x2[1]-p[1])*(x1[1]-p[1]))/r1/r2) >= 1.0 ) return;
  if ( (th = ont*acos(th)) == 0.0 ) return;
  assert(!std::isnan(th));

  if ( p[2] != 0.0 ) 
    subtriangle_polar_quad(r1, r2, th, p[2], k, &J1);
	  
  *I += J1; 
  
  return;
}


void triangle_potential_quad(double *x1, double *x2, double *x3,
                             double *p, 
                             double o12, double o23, double o31,
                             int k, 
                             double *I)
{
  double r1, r2, r3 ;  
  *I = 0.0 ;
  
  r1 = sqrt((p[0]-x1[0])*(p[0]-x1[0]) +  (p[1]-x1[1])*(p[1]-x1[1])) ;
  r2 = sqrt((p[0]-x2[0])*(p[0]-x2[0]) +  (p[1]-x2[1])*(p[1]-x2[1])) ;
  r3 = sqrt((p[0]-x3[0])*(p[0]-x3[0]) +  (p[1]-x3[1])*(p[1]-x3[1])) ;

  subtriangle_potential_quad(p, x1, x2, r1, r2, o12, k, I) ;
  subtriangle_potential_quad(p, x2, x3, r2, r3, o23, k, I) ;
  subtriangle_potential_quad(p, x3, x1, r3, r1, o31, k, I) ;

}


double triangle_quad_shape0(double *x1, double *x2, double *x3,
                            double *p, 
                            double o12, double o23, double o31, int k)
{
  double I;

  triangle_potential_quad(x1, x2, x3, p, o12, o23, o31, k, &I) ;

  return -p[2] * I;
}



void triangle_orientations(double *x1, double *x2, double *x3,
			   double *p, double *o12, double *o23, 
			   double *o31)
{
  double pa1[] = {p[0], p[1]};
  double pb1[] = {x1[0], x1[1]};
  double pc1[] = {x2[0], x2[1]};  
  *o12 = orient2d(pa1, pb1, pc1); 

  double pa2[] = {p[0], p[1]};
  double pb2[] = {x2[0], x2[1]};
  double pc2[] = {x3[0], x3[1]};
  *o23 = orient2d(pa2, pb2, pc2);

  double pa3[] = {p[0], p[1]};
  double pb3[] = {x3[0], x3[1]};
  double pc3[] = {x1[0], x1[1]};
  *o31 = orient2d(pa3, pb3, pc3);
}



void transformation_rotation_matrix(double *x1, double *x2, double *x3,
                                    double *A) 
{
  A[0] = x2[0] - x1[0] ; A[1] = x2[1] - x1[1] ; A[2] = x2[2] - x1[2] ; 
  
  A[6] = A[1]*(x3[2] - x2[2]) - A[2]*(x3[1] - x2[1]) ;
  A[7] = A[2]*(x3[0] - x2[0]) - A[0]*(x3[2] - x2[2]) ;
  A[8] = A[0]*(x3[1] - x2[1]) - A[1]*(x3[0] - x2[0]) ;

  {
    double _len = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]) ;
    A[0] /= _len ; A[1] /= _len ; A[2] /= _len ; 
    _len = sqrt(A[6]*A[6] + A[7]*A[7] + A[8]*A[8]) ;
    A[6] /= _len ; A[7] /= _len ; A[8] /= _len ; 
  }

  A[3] = A[7]*A[2] - A[8]*A[1] ;
  A[4] = A[8]*A[0] - A[6]*A[2] ;
  A[5] = A[6]*A[1] - A[7]*A[0] ;
}



// Function to be called by:
// computePhi(const Point3& pt, int k, const TriangleMesh& triMesh)
double 
potential_phi(const Point3& pt, int k, const Triangle& tri) {
  Point3 v0;
  tri.getVertex0(v0);

  Point3 v1;
  tri.getVertex1(v1);

  Point3 v2;
  tri.getVertex2(v2);

  // if the field point lie in the support plane of the triangle
  // then the contribution from this triangle is 0
  double pa[3] = {pt.x, pt.y, pt.z};
  double pb[3] = {v0.x, v0.y, v0.z};
  double pc[3] = {v1.x, v1.y, v1.z};
  double pd[3] = {v2.x, v2.y, v2.z};
  double orient = orient3d(pa, pb, pc, pd);
  bool inPlane = (fabs(orient) <= EPSILON);
  if (inPlane) return 0.0;


  double x1[3] = {v0.x, v0.y, v0.z};
  double x2[3] = {v1.x, v1.y, v1.z};
  double x3[3] = {v2.x, v2.y, v2.z};

  double p[3] = {pt.x, pt.y, pt.z};


  // transformation matrix
  double A[9];
  transformation_rotation_matrix(x1, x2, x3, A);


  // transformed point and vertex coordinates
  double y[3], y1[3], y2[3], y3[3];
  y[0] = A[0]*(p[0]-x1[0]) + A[1]*(p[1]-x1[1]) + A[2]*(p[2]-x1[2]) ;
  y[1] = A[3]*(p[0]-x1[0]) + A[4]*(p[1]-x1[1]) + A[5]*(p[2]-x1[2]) ;
  y[2] = A[6]*(p[0]-x1[0]) + A[7]*(p[1]-x1[1]) + A[8]*(p[2]-x1[2]) ;
  
  y2[0] = A[0]*(x2[0]-x1[0]) + A[1]*(x2[1]-x1[1]) + A[2]*(x2[2]-x1[2]) ;
  y2[1] = A[3]*(x2[0]-x1[0]) + A[4]*(x2[1]-x1[1]) + A[5]*(x2[2]-x1[2]) ;
  y2[2] = A[6]*(x2[0]-x1[0]) + A[7]*(x2[1]-x1[1]) + A[8]*(x2[2]-x1[2]) ;
  
  y3[0] = A[0]*(x3[0]-x1[0]) + A[1]*(x3[1]-x1[1]) + A[2]*(x3[2]-x1[2]) ;
  y3[1] = A[3]*(x3[0]-x1[0]) + A[4]*(x3[1]-x1[1]) + A[5]*(x3[2]-x1[2]) ;
  y3[2] = A[6]*(x3[0]-x1[0]) + A[7]*(x3[1]-x1[1]) + A[8]*(x3[2]-x1[2]) ;
	
  y1[0] = y1[1] = y1[2] = 0.0 ;
     
     
  // triangles orientation
  double o12, o23, o31; 
  triangle_orientations(y1, y2, y3, y, &o12, &o23, &o31);


  double result = triangle_quad_shape0(y1, y2, y3, y, o12, o23, o31, k);

  return result;
}


