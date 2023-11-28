#include "../Matrix.hpp"


Matrix getMatrix(std::ifstream& input) {
    int mwidth;
    input >> mwidth;
    int mheight;
    input >> mheight;

    mtrx m1;
    for (int i = 0; i < mheight; i++) {
        vc_dbl tmp;
        for (int j = 0; j < mwidth; j++) {
            double a;
            input >> a;
            tmp.push_back(a);
        }

        m1.push_back(tmp);
    }

    Matrix mt1(m1, mwidth, mheight);
    return mt1;
}

int main() {
    std::ifstream input;
    input.open("matrix.txt");

    // std::ifstream in;
    // in.open("mt.txt");

    Matrix mt1 = getMatrix(input);
    // Matrix mt2 = getMatrix(in);

    mt1.show();
    std::cout << "\n";
    // mt2.show();

    // Matrix mul = mt1 * mt2;
    
    // std::cout << "\n";

    // mul.show();

    mt1.QRFactorization();

    mt1.showQR();

    vc_dbl free_vc;
    free_vc.push_back(6);
    free_vc.push_back(12);
    free_vc.push_back(24);

    vc_dbl tt =  mt1.SolveSystemQR(free_vc);

    std::cout << "\nResult:\n";
    for(auto i : tt) std::cout << i << " ";
}