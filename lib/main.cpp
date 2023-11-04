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
  int melemcnt = 9;
  int msize = 3;

  // generate matrix
  // system("python3 randmtrx.py");

  std::ifstream fin;
  fin.open("matrix.txt");

  for (int i = 0; i < melemcnt; i++) {
    double a;
    fin >> a;
    vc1.push_back(a);
  }

  Matrix m1(vc1, msize);

  // -----------show_mtrx----------------

  std::cout << "Our Matrix: ";
  empty_line();
  m1.show();
  space();

  // std::cout << "Our Matrix Transp: ";
  // m1.Transponate();
  // m1.CheckIFDegenerateAndAddLU();
  // empty_line();
  // m1.show();
  // space();

  std::cout << "Determinant: " << m1.det_;

  space();

  if (m1.det_ == 0) return 0;

  // -----------------------------------

  // Matrix inv = m1.InvertibleMatrix();
  // std::cout << "Invertible Matrix: ";
  // empty_line();

  // inv.show();
  // space();

  // std::cout << "L matrix: ";
  // empty_line();
  // for (int i = 0; i < msize; i++) {
  //   for (int j = 0; j < msize; j++) {
  //     std::cout << m1.lu_.l_[i][j] << " ";
  //   }
  //   empty_line();
  // }

  // space();

  // std::cout << "U matrix: ";
  // empty_line();
  // for (int i = 0; i < msize; i++) {
  //   for (int j = 0; j < msize; j++) {
  //     std::cout << m1.lu_.u_[i][j] << " ";
  //   }
  //   empty_line();
  // }

  // ---------------------------------------

  // space();
  // double normA = m1.Norm();
  // double normInv = inv.Norm();

  // double cond = normA * normInv;

  // std::cout << "Norm of our mtrx, norm of invertible, cond number : ";
  // empty_line();
  // std::cout << normA << " " << normInv << " " << cond;
  // empty_line();

  // ------------------------------------------

  vc_dbl test;
  test.push_back(531);
  test.push_back(-460);
  test.push_back(193);
  // test.push_back(6);

  // vc_dbl res = m1.SolveSystem(test);

  // space();
  // std::cout << "Result: ";
  // empty_line();
  // for (auto i : res) std::cout << i << " ";
  // empty_line();

  Matrix u = m1.KholetskiyFactorization();

  space();
  std::cout << "Khol_mtrx: ";
  empty_line();
  u.show();
  space();

  vc_dbl res = m1.SolveSystemViaKholetskiy(test, u);
  std::cout << "Result: ";
  empty_line();
  for (auto i : res) std::cout << i << " ";
  empty_line();
}