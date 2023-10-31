#pragma once
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
  mtrx data_;
  double det_;
  LU lu_;

 public:
  Matrix() = default;
  Matrix(vc_dbl&, int);
  Matrix(int);
  Matrix(const Matrix&);
  Matrix(mtrx&, int);
  ~Matrix() = default;

  double Norm();
  double Det();
  Matrix InvertibleMatrix();
  std::pair<Matrix, Matrix> LU();
  bool AddLU();
  std::pair<Matrix, Matrix> QR();

  Matrix operator-() const;
  Matrix operator+(const Matrix&);
  Matrix operator-(const Matrix&);
  Matrix operator*(const Matrix&);

  bool operator==(const Matrix&) const;
  bool operator!=(const Matrix&) const;

  vc_dbl SolveL(const vc_dbl&);
  vc_dbl SolveU(const vc_dbl&);
  vc_dbl SolveSystem(const vc_dbl&);

  void show();
};
