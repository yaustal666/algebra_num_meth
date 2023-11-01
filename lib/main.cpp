#include <cmath>
#include <fstream>
#include <iostream>

#include "Matrix.hpp"

void space() { std::cout << std::endl << std::endl; }
void empty_line() { std::cout << std::endl; }

double vnorm(const vc_dbl& vc) {
  double res = 0;
  for (auto i : vc) {
    res += i * i;
  }

  res = sqrt(res);
  return res;
}

int main() {
  vc_dbl vc1;
  int melemcnt = 625;
  int msize = 25;

  // generate matrix
  system("python3 randmtrx.py");

  std::ifstream fin;
  fin.open("matrix.txt");

  for (int i = 0; i < melemcnt; i++) {
    double a;
    fin >> a;
    vc1.push_back(a);
  }

  Matrix m1(vc1, msize);

  // -----------show_mtrx----------------

  std::cout << "Our Matrix (Vandermond): ";
  empty_line();
  m1.show();
  space();

  std::cout << "Our Matrix Transp: ";
  m1.Transponate();
  m1.CheckIFDegenerateAndAddLU();
  empty_line();
  m1.show();
  space();

  std::cout << "Determinant: " << m1.det_;

  space();

  if (m1.det_ == 0) return 0;

  // -----------------------------------

  Matrix inv = m1.InvertibleMatrix();
  std::cout << "Invertible Matrix: ";
  empty_line();

  inv.show();
  space();

  std::cout << "L matrix: ";
  empty_line();
  for (int i = 0; i < msize; i++) {
    for (int j = 0; j < msize; j++) {
      std::cout << m1.lu_.l_[i][j] << " ";
    }
    empty_line();
  }

  space();

  std::cout << "U matrix: ";
  empty_line();
  for (int i = 0; i < msize; i++) {
    for (int j = 0; j < msize; j++) {
      std::cout << m1.lu_.u_[i][j] << " ";
    }
    empty_line();
  }

  // ---------------------------------------

  space();
  double normA = m1.Norm();
  double normInv = inv.Norm();

  double cond = normA * normInv;

  std::cout << "Norm of our mtrx, norm of invertible, cond number : ";
  empty_line();
  std::cout << normA << " " << normInv << " " << cond;
  empty_line();

  // ------------------------------------------

  vc_dbl coef;
  for (int i = 0; i < msize; i++) coef.push_back(1);

  space();

  vc_dbl tut = m1 * coef;

  std::cout << "Vector we made by Our Matrix * vector of 1 : ";
  empty_line();
  for (auto i : tut) std::cout << i << " ";

  empty_line();
  empty_line();

  vc_dbl result = m1.SolveSystem(tut);
  for (auto i : result) std::cout << i << " ";
  space();

  double norm1 = vnorm(coef);
  double norm2 = vnorm(result);
  double cringe = norm1 - norm2;
  double cringe2 = norm1 / norm2;
  std::cout << "Norm of 1 vc, norm of our solv vc, difference, division";
  empty_line();
  std::cout << norm1 << " " << norm2 << " " << cringe << " " << cringe2
            << std::endl;
  empty_line();
}