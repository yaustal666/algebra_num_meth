#include <cmath>
#include <fstream>
#include <iostream>

#include "Matrix.hpp"

int main() {
  //   // generate matrixes
  // system("python3 lab.py");
  int melemcnt = 16;
  int msize = 4;

  std::ifstream fin;
  std::ofstream fout;
  fin.open("matrix.txt");
  fout.open("res.txt");

  std::vector<double> vc1, vc2, vc3, vc4, vc5;

  for (int i = 0; i < melemcnt; i++) {
    double a;
    fin >> a;
    vc1.push_back(a);
  }

  // for (int i = 0; i < melemcnt; i++) {
  //   double a;
  //   fin >> a;
  //   vc2.push_back(a);
  // }

  //   for (int i = 0; i < 625; i++) {
  //     double a;
  //     fin >> a;
  //     vc3.push_back(a);
  //   }

  //   for (int i = 0; i < 625; i++) {
  //     double a;
  //     fin >> a;
  //     vc4.push_back(a);
  //   }

  //   for (int i = 0; i < 625; i++) {
  //     double a;
  //     fin >> a;
  //     vc5.push_back(a);
  //   }

  //   Matrix m1(vc1, 25);
  //   Matrix m2(vc2, 25);
  //   Matrix m3(vc3, 25);
  //   Matrix m4(vc4, 25);
  //   Matrix m5(vc5, 25);

  //   double n1, n2, n3, n4, n5;

  //   n1 = m1.norm();
  //   n2 = m2.norm();
  //   n3 = m3.norm();
  //   n4 = m4.norm();
  //   n5 = m5.norm();

  //   fout << "Norm 1 = " << n1 << "\n";
  //   fout << "Norm 2 = " << n2 << "\n";
  //   fout << "Norm 3 = " << n3 << "\n";
  //   fout << "Norm 4 = " << n4 << "\n";
  //   fout << "Norm 5 = " << n5 << "\n";

  Matrix m1(vc1, msize);
  // Matrix m2(vc2, msize);

  // std::pair<Matrix, Matrix> lu = m1.LU();

  m1.show();
  std::cout << "\n\n";
  vc_dbl vk;
  vk.push_back(4);
  vk.push_back(44);
  vk.push_back(30);
  vk.push_back(8);

  vc_dbl vc = m1.SolveU(vk);

  std::cout << m1.det_ << " ";
  std::cout << std::endl;

  for (auto i : vc) {
    std::cout << i << " ";
  }
}