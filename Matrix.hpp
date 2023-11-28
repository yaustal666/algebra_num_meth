#pragma once

#include "Vec.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

typedef std::vector<std::vector<double>> mtrx;
typedef std::vector<double> vc_dbl;

struct LU {
    mtrx l_;
    mtrx u_;
};

struct Kholetskiy {
    mtrx u_;
    mtrx ut_;
};

struct QR {
    mtrx q_;
    mtrx r_;
};

class Matrix {
   private:
    mtrx data_;

    int width_ = 0;
    int height_ = 0;

    LU lu_;
    Kholetskiy khol_;
    QR qr_;

   public:
    Matrix() = default;
    Matrix(const mtrx&, int, int);

    // Diagonal Matrix
    Matrix(const vc_dbl&, int, int);
    Matrix(double, int, int);

    Matrix(const Matrix&);

    ~Matrix() = default;

    Matrix operator+(const Matrix&);
    Matrix operator-(const Matrix&);
    Matrix operator*(const Matrix&);

    Matrix operator-();

    Matrix operator*(double);
    Matrix operator/(double);

    bool operator==(const Matrix&);
    bool operator!=(const Matrix&);

    double Norm();
    double Det();

    Matrix InvertibleMatrix();
    Matrix Transponate() const;

    void LUFactorization();
    void KholetskiyFactorization();
    void QRFactorization();

    vc_dbl SolveL(const vc_dbl&, const mtrx&);
    vc_dbl SolveU(const vc_dbl&, const mtrx&);

    vc_dbl SolveSystemLU(const vc_dbl&);
    vc_dbl SolveSystemKholetskiy(const vc_dbl&);
    vc_dbl SolveSystemQR(const vc_dbl&);

    void show();
    void showLU();
    void showKholetskiy();
    void showQR();

    bool IsDegenerate();

    double get_det();
    int get_width();
    int get_height();
    vc_dbl getLine(int);
    double get_element(int, int);
};
