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
}

Matrix::Matrix(mtrx& m, int size) {
  size_ = size;
  data_ = m;
}

Matrix::Matrix(int size) {
  size_ = size;
  data_.resize(size_);
  for (int i = 0; i < size_; i++) data_[i].resize(size_, 0);
}

Matrix::Matrix(const Matrix& other) {
  size_ = other.size_;
  data_ = other.data_;
}
// --------------------------------

double Matrix::norm() {
  double res = 0;

  for (auto i : data_) {
    for (auto j : i) {
      res += j * j;
    }
  }

  return sqrt(res);
}

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

        if (U.data_[i][j] != U.data_[i][j]) {
          std::cout << "Вырожденная матрица, упс!";
          return std::make_pair(Matrix(1), Matrix(1));
        }

      } else {
        double tmp2 = 0;

        for (int k = 0; k <= j - 1; k++) {
          tmp2 += L.data_[i][k] * U.data_[k][j];
        }

        L.data_[i][j] = (1 / U.data_[j][j]) * (data_[i][j] - tmp2);
      }
    }
  }

  return std::make_pair(L, U);
}

double Matrix::det() { double res; }

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

// --------------util---------------
void Matrix::show() {
  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      std::cout << data_[i][j] << " ";
    }
    std::cout << std::endl;
  }
}
