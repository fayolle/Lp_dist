#include "LinearAlgebra.h"


double det(const Matrix33& m) {
  return -(m(0,2)*m(1,1)*m(2,0)) + 
    (m(0,1)*m(1,2)*m(2,0)) + 
    (m(0,2)*m(1,0)*m(2,1)) - 
    (m(0,0)*m(1,2)*m(2,1)) - 
    (m(0,1)*m(1,0)*m(2,2)) + 
    (m(0,0)*m(1,1)*m(2,2));
}


Matrix33 transpose(const Matrix33& m) {
  Matrix33 trans;
  for (unsigned row = 0; row < 3; ++row) {
    for (unsigned col = 0; col < 3; ++col) {
      trans(row,col) = m(col,row);
    }
  }
  return trans;
}


// concat matrices a (3x2) and b (3x1) into a matrix of size 3x3
// (append the colum of b to those of a)
void concatMatrices(const Matrix32& a, const Matrix31& b, Matrix33& concat) {
  concat(0,0) = a(0,0);
  concat(1,0) = a(1,0);
  concat(2,0) = a(2,0);
  concat(0,1) = a(0,1);
  concat(1,1) = a(1,1);
  concat(2,1) = a(2,1);
  concat(0,2) = b(0,0);
  concat(1,2) = b(1,0);
  concat(2,2) = b(2,0);
}

