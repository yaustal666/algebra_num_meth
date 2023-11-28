#include "../Matrix.hpp"

int main() {

    std::ifstream input;
    input.open("matrix.txt");

    int mwidth;
    input >> mwidth;
    int mheight;
    input >> mheight;


    mtrx m1;
    for(int i = 0; i < mheight; i++) {
        vc_dbl tmp;
        for(int j = 0; j < mwidth; j++) {
            double a;
            input >> a;
            tmp.push_back(a);
        }

        m1.push_back(tmp);
    }

    Matrix mt1(m1, mwidth, mheight);

    mt1.KholetskiyFactorization();

    std::cout << "Given Matrix: " << "\n";
    mt1.show();
    std::cout << "\nKholetskiy factorization: " << std::endl;
    mt1.showKholetskiy();

    vc_dbl freevc1;
    freevc1.push_back(531);
    freevc1.push_back(-460);
    freevc1.push_back(193);

    vc_dbl res = mt1.SolveSystemKholetskiy(freevc1);

    for(auto i : res) std::cout << i << " ";

}