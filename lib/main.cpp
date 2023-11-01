#include <cmath>
#include <fstream>
#include <iostream>

#include "Matrix.hpp"

void empty_line() { std::cout << std::endl; }

int main() {
  vc_dbl vc1;

  system("python3 randmtrx.py");

  std::ifstream fin;
  fin.open("matrix.txt");

  int melemcnt = 625;
  int msize = 25;

  for (int i = 0; i < melemcnt; i++) {
    double a;
    fin >> a;
    vc1.push_back(a);
  }

  Matrix m1(vc1, msize);
  Matrix inv = m1.InvertibleMatrix();

  double normA = m1.Norm();
  double normInv = inv.Norm();

  double cond = normA * normInv;

  std::cout << normA << " " << normInv << " " << cond;
  empty_line();
}