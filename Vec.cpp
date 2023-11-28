#include "Vec.hpp"

Vec::Vec(vc_dbl& vc) {
    data_ = vc;
    size_ = vc.size();
}

Vec::Vec(const Vec& other) {
    data_ = other.data_;
    size_ = other.size_;
}

double Vec::operator*(const Vec& other) {
    double res = 0;

    for (int i = 0; i < size_; i++) {
        res += data_[i] * other.data_[i];
    }

    return res;
}

Vec Vec::operator*(double num) {
    Vec vc(*this);

    for (int i = 0; i < size_; i++) {
        vc.data_[i] *= num;
    }

    return vc;
}

Vec Vec::operator+(const Vec& other) {
    Vec vc(*this);

    for (int i = 0; i < size_; i++) {
        vc.data_[i] += other.data_[i];
    }

    return vc;
}

Vec Vec::operator-(const Vec& other) {
    Vec vc(*this);

    for (int i = 0; i < size_; i++) {
        vc.data_[i] -= other.data_[i];
    }

    return vc;
}

void Vec::show() {
    for (auto i : data_) std::cout << i << " ";
    std::cout << "\n";
}

vc_dbl Vec::get_data() { return data_; }

void Vec::Normalize() {
        double norm = 0;
        for(auto i : data_) norm += i*i;

        norm = sqrt(norm);

        for(int i = 0; i < size_; i++) data_[i] /= norm;
    }