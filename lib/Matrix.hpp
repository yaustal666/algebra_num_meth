#pragma once
#include <iostream>
#include <vector>

typedef std::vector<std::vector<double>> mtrx;
typedef std::vector<double> vc_dbl;

struct LU {
  mtrx l_;
  mtrx u_;
};

class Matrix {
 public:
  int size_;
  double det_ = 0;

  mtrx data_;

  LU lu_;

 public:
  Matrix() = default;

  /// @brief creates a matrix from vector and size,
  /// vector split on n = size parts and result is n x n matrix
  Matrix(vc_dbl&, int);

  /// @brief create matrix from 2dim vector
  Matrix(mtrx&, int);

  Matrix(const Matrix&);

  /// @brief create zero matrix of n size
  Matrix(int);

  /// @brief create n x n diagonal matrix with value
  Matrix(int, double);

  ~Matrix() = default;

  double Norm();
  double Det();
  Matrix Transponate() const;

  std::pair<Matrix, Matrix> LUFactorization();
  Matrix KholetskiyFactorization();

  Matrix operator-() const;
  Matrix operator+(const Matrix&);
  Matrix operator-(const Matrix&);
  Matrix operator*(const Matrix&);
  vc_dbl operator*(const vc_dbl&);

  Matrix operator+(double);
  Matrix operator-(double);
  Matrix operator*(double);

  bool operator==(const Matrix&) const;
  bool operator!=(const Matrix&) const;

  vc_dbl SolveL(const vc_dbl&, const mtrx&);
  vc_dbl SolveU(const vc_dbl&, const mtrx&);
  vc_dbl SolveSystemViaLU(const vc_dbl&);
  vc_dbl SolveSystemViaKholetskiy(const vc_dbl&, const Matrix&);

  Matrix InvertibleMatrix();

  void show();
  bool IsDegenerate();

 public:
  bool CheckIFDegenerateAndAddLU();
};
