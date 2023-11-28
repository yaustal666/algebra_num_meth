#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

typedef std::vector<std::vector<double>> mtrx;
typedef std::vector<double> vc_dbl;


class Vec {
    private:
    vc_dbl data_;
    int size_;

    public:

    Vec() = default;

    Vec(vc_dbl&);
    Vec(const Vec&);

    ~Vec() = default;

    double operator*(const Vec&);
    Vec operator+(const Vec&);
    Vec operator-(const Vec&);

    Vec operator*(double);

    void Normalize();

    vc_dbl get_data();

    void show();

};
