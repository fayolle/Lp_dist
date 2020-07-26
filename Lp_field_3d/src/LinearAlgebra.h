#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H


// Matrix type
template <unsigned ROW, unsigned COL>
class Matrix {
public:
  Matrix() {
    for (unsigned row = 0; row < ROW; ++row) {
      for (unsigned col = 0; col < COL; ++col) {
	m[row][col] = 0.0;
      }
    }
  }

  double& operator() (unsigned row, unsigned col) {
    return m[row][col];
  }

  const double& operator() (unsigned row, unsigned col) const {
    return m[row][col];
  }
  
  void permuteColumns(unsigned i, unsigned j) {
    for (unsigned row = 0; row < ROW; ++row) {
      double temp = m[row][i];
      m[row][i] = m[row][j];
      m[row][j] = temp;
    }
  }

    
private:
  double m[ROW][COL];
};


// Instantiate matrices with useful dimension
typedef Matrix<3, 3> Matrix33;
typedef Matrix<3, 2> Matrix32;
typedef Matrix<3, 1> Matrix31;


// Operations
double det(const Matrix33& m);
Matrix33 transpose(const Matrix33& m);
void concatMatrices(const Matrix32& a, const Matrix31& b, Matrix33& concat);



#endif
