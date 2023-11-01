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
  void Transponate();

  std::pair<Matrix, Matrix> LUFactorization();

  Matrix operator-() const;
  Matrix operator+(const Matrix&);
  Matrix operator-(const Matrix&);
  Matrix operator*(const Matrix&);

  Matrix operator+(double);
  Matrix operator-(double);
  Matrix operator*(double);

  bool operator==(const Matrix&) const;
  bool operator!=(const Matrix&) const;

  vc_dbl SolveL(const vc_dbl&);
  vc_dbl SolveU(const vc_dbl&);
  vc_dbl SolveSystem(const vc_dbl&);

  Matrix InvertibleMatrix();

  void show();
  bool IsDegenerate();

 private:
  bool CheckIFDegenerateAndAddLU();
};
