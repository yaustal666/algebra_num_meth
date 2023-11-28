#include "Matrix.hpp"

Matrix::Matrix(const mtrx& m, int width, int height) {
    width_ = width;
    height_ = height;

    data_ = m;
}

Matrix::Matrix(const vc_dbl& diag, int width, int height) {
    // need to check is vc size = width

    width_ = width;
    height_ = height;

    data_.resize(height);

    for (auto i : data_) i.resize(width, 0);

    for (int i = 0; i < width; i++) {
        data_[i][i] = diag[i];
    }
}

Matrix::Matrix(double diag, int width, int height) {
    width_ = width;
    height_ = height;

    for(int i = 0; i < height_; i++) {
        vc_dbl tmp;
        for (int j = 0; j < width_; j++)
        {
            tmp.push_back(0);
        }
        data_.push_back(tmp);        
    }

    int m = height_ < width_ ? height_ : width_;
    for(int i = 0; i < m; i++) {
        data_[i][i] = diag;
    }
}

Matrix::Matrix(const Matrix& other) {
    width_ = other.width_;
    height_ = other.height_;

    data_ = other.data_;
}

//----------------------------------------------------

Matrix Matrix::operator+(const Matrix& other) {
    Matrix res(other);

    for (int i = 0; i < height_; i++) {
        for (int j = 0; j < width_; j++) {
            res.data_[i][j] += data_[i][j];
        }
    }

    return res;
}

Matrix Matrix::operator-(const Matrix& other) {
    Matrix res(other);

    for (int i = 0; i < height_; i++) {
        for (int j = 0; j < width_; j++) {
            res.data_[i][j] -= data_[i][j];
        }
    }

    return res;
}

Matrix Matrix::operator*(const Matrix& other) {
    Matrix res(0, other.width_, height_);

    for (int i = 0; i < res.height_; i++) {
        for (int j = 0; j < res.width_; j++) {
            double tmp = 0;
            for (int k = 0; k < width_; k++) {
                tmp += data_[i][k] * other.data_[k][j];
            }
            if (std::abs(tmp) < 0.00000001) tmp = 0;
            res.data_[i][j] = tmp;
        }
    }

    return res;
}

Matrix Matrix::operator-() {
    Matrix res(*this);

    for (int i = 0; i < height_; i++) {
        for (int j = 0; j < width_; j++) {
            res.data_[i][j] = -res.data_[i][j];
        }
    }

    return res;
}

Matrix Matrix::operator*(double num) {
    Matrix res(*this);

    for (int i = 0; i < height_; i++) {
        for (int j = 0; j < width_; j++) {
            res.data_[i][j] *= num;
        }
    }

    return res;
}

Matrix Matrix::operator/(double num) {
    Matrix res(*this);

    for (int i = 0; i < height_; i++) {
        for (int j = 0; j < width_; j++) {
            res.data_[i][j] /= num;
        }
    }

    return res;
}

bool Matrix::operator==(const Matrix& other) {
    for (int i = 0; i < height_; i++) {
        for (int j = 0; j < width_; j++) {
            if (data_[i][j] != other.data_[i][j]) return false;
        }
    }

    return true;
}

bool Matrix::operator!=(const Matrix& other) { return !(*this == other); }

//-----------------------------------------------

void Matrix::LUFactorization() {
    Matrix L(0, width_, height_);
    Matrix U(0, width_, height_);


    for (int j = 0; j < width_; j++) U.data_[0][j] = data_[0][j];
    for (int i = 0; i < height_; i++)
        L.data_[i][0] = data_[i][0] / U.data_[0][0];
    for (int i = 0; i < width_; i++) L.data_[i][i] = 1;

    for (int i = 1; i < height_; i++) {
        for (int j = 1; j < width_; j++) {
            if (i <= j) {
                double tmp1 = 0;

                for (int k = 0; k <= i - 1; k++) {
                    tmp1 += L.data_[i][k] * U.data_[k][j];
                }

                U.data_[i][j] = data_[i][j] - tmp1;

                // zero-determinant check
                if (U.data_[i][j] != U.data_[i][j]) {
                    std::cout << "Degenerate Matrix" << std::endl;
                    return;
                }

            } else {
                double tmp2 = 0;

                for (int k = 0; k <= j - 1; k++) {
                    tmp2 += L.data_[i][k] * U.data_[k][j];
                }

                L.data_[i][j] = (1 / U.data_[j][j]) * (data_[i][j] - tmp2);

                // zero-determinant check
                if (L.data_[i][j] != L.data_[i][j]) {
                    std::cout << "Degenerate Matrix" << std::endl;
                    return;
                }
            }
        }
    }

    lu_.l_ = L.data_;
    lu_.u_ = U.data_;
}

Matrix Matrix::Transponate() const {
    Matrix copy(0, height_, width_);

    for (int i = 0; i < copy.height_; i++) {
        for (int j = 0; j < copy.width_; j++) {
            copy.data_[i][j] = data_[j][i];
        }
    }

    return copy;
}

void Matrix::KholetskiyFactorization() {
    Matrix U(0, width_, height_);
    U.data_[0][0] = sqrt(data_[0][0]);

    for (int j = 1; j < width_; j++) {
        U.data_[0][j] = data_[0][j] / U.data_[0][0];
    }

    for (int i = 0; i < height_; i++) {
        for (int j = 1; j < width_; j++) {
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

    khol_.u_ = U.data_;

    Matrix ut = U.Transponate();

    khol_.ut_ = ut.data_;
}

void Matrix::QRFactorization() {
    std::vector<Vec> mat;

    Matrix mtrnsp = this->Transponate();

    for(int i = 0; i < width_; i++) {
        vc_dbl tmp = mtrnsp.getLine(i);
        Vec tmpvc(tmp);
        mat.push_back(tmpvc);
    }

    for(int i = 0; i < mtrnsp.get_width(); i++) {
        for(int j = 0; j < i; j++) {
            mat[i] = mat[i] - mat[j] * ((mat[j] * mat[i]) / (mat[j] * mat[j]));
        }
        mat[i].Normalize();
    }

    mtrx m2;
    for(auto i : mat) m2.push_back(i.get_data());

    Matrix mt2(m2, 3, 3);

    qr_.q_ = mt2.Transponate().data_;

    qr_.r_ = (mt2 * (*this)).data_;

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
    while (i < width_) {
        tmp *= lu_.l_[i][j];
        i++;
        j++;
    }

    res = tmp;

    i = 0;
    j = 0;
    tmp = 1;
    while (i < width_) {
        tmp *= lu_.u_[i][j];
        i++;
        j++;
    }

    res *= tmp;

    return res;
}

vc_dbl Matrix::SolveL(const vc_dbl& free_vc, const mtrx& lower) {
    vc_dbl res;
    res.push_back(free_vc[0] / lower[0][0]);

    for (int i = 1; i < width_; i++) {
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
    int last = width_ - 1;
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

vc_dbl Matrix::SolveSystemLU(const vc_dbl& free_vc) {
    vc_dbl tmp = this->SolveL(free_vc, lu_.l_);

    vc_dbl res = this->SolveU(tmp, lu_.u_);

    return res;
}

vc_dbl Matrix::SolveSystemKholetskiy(const vc_dbl& free_vc) {
    vc_dbl tmp = this->SolveL(free_vc, khol_.ut_);

    vc_dbl res = this->SolveU(tmp, khol_.u_);

    return res;
}

vc_dbl Matrix::SolveSystemQR(const vc_dbl& free_vc) { 
    mtrx frvc;
    frvc.push_back(free_vc);

    Matrix fr(frvc, free_vc.size(), 1);

    fr = fr.Transponate();

    Matrix q(qr_.q_, width_, height_);

    Matrix tmp = q.Transponate() * fr;

    vc_dbl tut;
    for(int i = 0; i < tmp.height_; i++) tut.push_back(tmp.data_[i][0]);

    vc_dbl res = this->SolveU(tut, qr_.r_);

    return res;
}

Matrix Matrix::InvertibleMatrix() {
    Matrix idm = Matrix(1, width_, height_);

    mtrx res_data;
    for (int i = 0; i < width_; i++) {
        vc_dbl tmp = this->SolveSystemLU(idm.data_[i]);
        res_data.push_back(tmp);
    }

    Matrix res(res_data, width_, height_);
    res.Transponate();

    return res;
}

void Matrix::show() {
    for (int i = 0; i < height_; i++) {
        for (int j = 0; j < width_; j++) {
            std::cout << data_[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void Matrix::showLU() {
    std::cout << "\nLower matrix: " << "\n";
    for(auto i : lu_.l_) {
        for(auto j : i) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nUpper matrix: " << "\n";

        for(auto i : lu_.u_) {
        for(auto j : i) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }
}

void Matrix::showKholetskiy() {
    std::cout << "\nLower matrix: " << "\n";
    for(auto i : khol_.ut_) {
        for(auto j : i) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nUpper matrix: " << "\n";

        for(auto i : khol_.u_) {
        for(auto j : i) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }
}

void Matrix::showQR() {
    std::cout << "\nOrthogonal Q matrix: " << "\n";
    for(auto i : qr_.q_) {
        for(auto j : i) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nUpper R matrix: " << "\n";

        for(auto i : qr_.r_) {
        for(auto j : i) {
            std::cout << j << " ";
        }
        std::cout << "\n";
    }
}

int Matrix::get_width() { return width_; }

vc_dbl Matrix::getLine(int i) { return data_[i]; }
