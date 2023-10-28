#pragma once
#include <vector>

typedef std::vector<std::vector<double>> mtrx;
typedef std::vector<double> vc_dbl;

class Matrix {
 private:
  int size_;
  mtrx data_;

 public:
  Matrix() = default;
  Matrix(vc_dbl&, int);
  Matrix(int);
  Matrix(const Matrix&);
  Matrix(mtrx&, int);
  ~Matrix() = default;

  double norm();
  double det();
  std::pair<Matrix, Matrix> LU();
  std::pair<Matrix, Matrix> QR();

  Matrix operator-() const;
  Matrix operator+(const Matrix&);
  Matrix operator-(const Matrix&);
  Matrix operator*(const Matrix&);

  void show();
};