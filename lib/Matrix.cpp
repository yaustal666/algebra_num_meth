#include "Matrix.hpp"

#include <algorithm>
#include <cmath>

// ------------constructors-------------------------------
Matrix::Matrix(vc_dbl& vc, int size) {
  size_ = size;
  for (int i = 0; i < size_; i++) {
    vc_dbl tmp;
    for (int j = i * size_; j < size_ * (i + 1); j++) {
      tmp.push_back(vc[j]);
    }
    data_.push_back(tmp);
  }

  bool check = CheckIFDegenerateAndAddLU();
  if (!check) {
    det_ = this->Det();
  }
}

Matrix::Matrix(mtrx& m, int size) {
  size_ = size;
  data_ = m;

  bool check = CheckIFDegenerateAndAddLU();
  if (!check) {
    det_ = this->Det();
  }
}

Matrix::Matrix(int size) {
  size_ = size;
  data_.resize(size_);
  for (int i = 0; i < size_; i++) data_[i].resize(size_, 0);
}

Matrix::Matrix(int size, double diag) {
  size_ = size;
  data_.resize(size_);
  for (int i = 0; i < size_; i++) data_[i].resize(size_, 0);

  for (int i = 0; i < size_; i++) data_[i][i] = diag;
}

Matrix::Matrix(const Matrix& other) {
  size_ = other.size_;
  data_ = other.data_;
  lu_ = other.lu_;
  det_ = other.det_;
}

// --------------algebra_operations--------------------------

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
      if (std::abs(tmp) < 0.00000001) tmp = 0;
      res.data_[i][j] = tmp;
    }
  }

  return res;
}

vc_dbl Matrix::operator*(const vc_dbl& other) {
  vc_dbl res;

  for (int i = 0; i < size_; i++) {
    double tmp = 0;
    for (int j = 0; j < size_; j++) {
      tmp += data_[i][j] * other[j];
    }
    res.push_back(tmp);
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

Matrix Matrix::operator+(double num) {
  Matrix res(*this);

  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      res.data_[i][j] += num;
    }
  }

  return res;
}

Matrix Matrix::operator-(double num) {
  Matrix res(*this);

  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      res.data_[i][j] -= num;
    }
  }

  return res;
}

Matrix Matrix::operator*(double num) {
  Matrix res(*this);

  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      res.data_[i][j] *= num;
    }
  }

  return res;
}

// --------------matrix_properties-----------------------
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

Matrix Matrix::Transponate() const {
  Matrix copy(*this);
  for (int i = 0; i < size_; i++) {
    for (int j = i; j < size_; j++) {
      double t = copy.data_[i][j];
      copy.data_[i][j] = copy.data_[j][i];
      copy.data_[j][i] = t;
    }
  }

  return copy;
}

Matrix Matrix::InvertibleMatrix() {
  Matrix idm = Matrix(size_, 1);

  mtrx res_data;
  for (int i = 0; i < size_; i++) {
    vc_dbl tmp = this->SolveSystemViaLU(idm.data_[i]);
    res_data.push_back(tmp);
  }

  Matrix res(res_data, size_);
  res.Transponate();

  return res;
}
// --------------Factorizations--------------------------
std::pair<Matrix, Matrix> Matrix::LUFactorization() {
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

Matrix Matrix::KholetskiyFactorization() {
  Matrix U(size_);
  U.data_[0][0] = sqrt(data_[0][0]);

  for (int j = 1; j < size_; j++) {
    U.data_[0][j] = data_[0][j] / U.data_[0][0];
  }

  for (int i = 1; i < size_; i++) {
    for (int j = 1; j < size_; j++) {
      if (i == j) {
        double tmp1 = 0;
        for (int k = 0; k < i; k++) {
          tmp1 += U.data_[k][i] * U.data_[k][i];
        }

        U.data_[i][i] = sqrt(data_[i][i] - tmp1);
      }

      if (i < j) {
        double tmp2 = 0;
        for (int k = 0; k < i; k++) {
          tmp2 += U.data_[k][i] * U.data_[k][j];
        }

        U.data_[i][j] = (data_[i][j] - tmp2) / U.data_[i][i];
      }
    }
  }

  return U;
}
// --------------solve_equasion_with_LU------------------
vc_dbl Matrix::SolveL(const vc_dbl& free_vc, const mtrx& lower) {
  vc_dbl res;
  res.push_back(free_vc[0] / lower[0][0]);

  for (int i = 1; i < size_; i++) {
    double tmp = free_vc[i];
    for (int j = 0; j < i; j++) {
      tmp -= res[j] * lower[i][j];
    }

    tmp /= lower[i][i];
    res.push_back(tmp);
  }

  return res;
}

vc_dbl Matrix::SolveU(const vc_dbl& free_vc, const mtrx& upper) {
  vc_dbl res;
  int last = size_ - 1;
  res.push_back(free_vc[last] / upper[last][last]);

  for (int i = last - 1; i >= 0; i--) {
    double tmp = free_vc[i];
    for (int j = last; j > i; j--) {
      tmp -= res[last - j] * upper[i][j];
    }

    tmp /= upper[i][i];

    res.push_back(tmp);
  }

  std::reverse(res.begin(), res.end());
  return res;
}

vc_dbl Matrix::SolveSystemViaLU(const vc_dbl& free_vc) {
  vc_dbl tmp = this->SolveL(free_vc, lu_.l_);

  vc_dbl res = this->SolveU(tmp, lu_.u_);

  return res;
}

vc_dbl Matrix::SolveSystemViaKholetskiy(const vc_dbl& free_vc,
                                        const Matrix& u) {
  Matrix ut = u.Transponate();

  vc_dbl tmp = this->SolveL(free_vc, ut.data_);

  vc_dbl res = this->SolveU(tmp, u.data_);

  return res;
}

// --------------Util------------------------------------
void Matrix::show() {
  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_; j++) {
      std::cout << data_[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

bool Matrix::CheckIFDegenerateAndAddLU() {
  std::pair<Matrix, Matrix> lu = this->LUFactorization();

  Matrix l = lu.first;
  Matrix z(1);
  if (l == z) {
    return 1;
  }

  lu_.l_ = lu.first.data_;
  lu_.u_ = lu.second.data_;
  return 0;
}

bool Matrix::IsDegenerate() { return det_ == 0; }