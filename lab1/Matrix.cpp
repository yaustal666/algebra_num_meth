#include "Matrix.hpp"

#include <cmath>
#include <fstream>
#include <iostream>

Matrix::Matrix(vc_dbl& vc, int size) {
  size_ = size;
  for (int i = 0; i < size_; i++) {
    std::vector<double> tmp;
    for (int j = i * size_; j < size_ * (i + 1); j++) {
      tmp.push_back(vc[j]);
    }
    data_.push_back(tmp);
  }

  bool k = this->AddLU();
  if (k)
    det_ = this->Det();
  else
    det_ = 0;
}

Matrix::Matrix(mtrx& m, int size) {
  size_ = size;
  data_ = m;

  bool k = this->AddLU();
  if (k) det_ = this->Det();
}

Matrix::Matrix(int size) {
  size_ = size;
  data_.resize(size_);
  for (int i = 0; i < size_; i++) data_[i].resize(size_, 0);
}

Matrix::Matrix(const Matrix& other) {
  size_ = other.size_;
  data_ = other.data_;
  lu_ = other.lu_;
}
// --------------------------------

std::pair<Matrix, Matrix> Matrix::LU() {
  Matrix L(size_);
  Matrix U(size_);

  for (int j = 0; j < size_; j++) U.data_[0][j] = data_[0][j];
  for (int i = 0; i < size_; i++) L.data_[i][0] = data_[i][0] / U.data_[0][0];
  for (int i = 0; i < size_; i++) L.data_[i][i] = 1;

  for (int i = 1; i < size_; i++) {
    for (int j = 1; j < size_; j++) {
      if (i <= j) {
        double tmp1 = 0;

        for (int k = 0; k <= i - 1; k++) {
          tmp1 += L.data_[i][k] * U.data_[k][j];
        }

        U.data_[i][j] = data_[i][j] - tmp1;

        // zero-determinant check
        if (U.data_[i][j] != U.data_[i][j]) {
          return std::make_pair(Matrix(1), Matrix(1));
        }

      } else {
        double tmp2 = 0;

        for (int k = 0; k <= j - 1; k++) {
          tmp2 += L.data_[i][k] * U.data_[k][j];
        }

        L.data_[i][j] = (1 / U.data_[j][j]) * (data_[i][j] - tmp2);

        // zero-determinant check
        if (L.data_[i][j] != L.data_[i][j]) {
          return std::make_pair(Matrix(1), Matrix(1));
        }
      }
    }
  }

  return std::make_pair(L, U);
}

double Matrix::Norm() {
  double res = 0;

  for (auto i : data_) {
    for (auto j : i) {
      res += j * j;
    }
  }

  return sqrt(res);
}

double Matrix::Det() {
  double res = 0;

  int i = 0, j = 0;
  double tmp = 1;
  while (i < size_) {
    tmp *= lu_.l_[i][j];
    i++;
    j++;
  }

  res = tmp;

  i = 0;
  j = 0;
  tmp = 1;
  while (i < size_) {
    tmp *= lu_.u_[i][j];
    i++;
    j++;
  }

  res *= tmp;

  return res;
}

// Matrix Matrix::InvertibleMatrix() {}
// ------------operators------------

Matrix Matrix::operator-() const {
  Matrix cp(*this);

  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      cp.data_[i][j] = -cp.data_[i][j];
    }
  }

  return cp;
}

Matrix Matrix::operator+(const Matrix& other) {
  Matrix res(*this);

  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      res.data_[i][j] += other.data_[i][j];
    }
  }

  return res;
}

Matrix Matrix::operator-(const Matrix& other) {
  Matrix res(*this);

  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      res.data_[i][j] -= other.data_[i][j];
    }
  }

  return res;
}

Matrix Matrix::operator*(const Matrix& other) {
  Matrix res(*this);

  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      double tmp = 0;
      for (int k = 0; k < size_; k++) {
        tmp += data_[i][k] * other.data_[k][j];
      }
      res.data_[i][j] = tmp;
    }
  }

  return res;
}

bool Matrix::operator==(const Matrix& other) const {
  if (size_ != other.size_) return false;

  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      if (data_[i][j] != other.data_[i][j]) return false;
    }
  }

  return true;
}

bool Matrix::operator!=(const Matrix& other) const { return !(*this == other); }

// ---------solving-systems---------
vc_dbl Matrix::SolveL(const vc_dbl& free_vc) {
  vc_dbl res;
  res.push_back(free_vc[0] / lu_.l_[0][0]);

  for (int i = 1; i < size_; i++) {
    double tmp = free_vc[i];
    for (int j = 0; j < i; j++) {
      tmp -= res[j] * lu_.l_[i][j];
    }
    tmp /= lu_.l_[i][i];
    res.push_back(tmp);
  }

  return res;
}

vc_dbl Matrix::SolveU(const vc_dbl& free_vc) {
  vc_dbl res;
  int last = size_ - 1;
  res.push_back(free_vc[last] / lu_.u_[last][last]);

  for (int i = last - 1; i >= 0; i--) {
    double tmp = free_vc[i];
    for (int j = 0; j < last - i; j++) {
      tmp -= res[j] * lu_.u_[i][j];
    }
    tmp /= lu_.u_[i][i];
    res.push_back(tmp);
  }

  return res;
}
// vc_dbl Matrix::SolveSystem(const vc_dbl& free_vc) {}

// --------------util---------------
void Matrix::show() {
  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      std::cout << data_[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

bool Matrix::AddLU() {
  std::pair<Matrix, Matrix> lu = this->LU();

  Matrix l = lu.first;
  Matrix z(1);
  if (l == z) {
    this->show();
    std::cout << "^^^^^ эта матрица вырожденная!" << std::endl;
    return 0;
  }

  lu_.l_ = lu.first.data_;
  lu_.u_ = lu.second.data_;
  return 1;
}