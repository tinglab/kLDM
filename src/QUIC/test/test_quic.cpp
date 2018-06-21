#include "../../utils.h"
#include "../../utils.cpp"
#include "../QUIC.h"
#include "../QUIC.cpp"

int main(int argc, char** argv) {
    MatrixXd temp = MatrixXd::Random(5, 5);
    MatrixXd S = temp;
    S = temp.transpose() + MatrixXd::Identity(5, 5) * 5;
    std::cout << "S: " << S << std::endl;
    MatrixXd Theta = MatrixXd::Zero(5, 5);
    QUICC(S, 0.5, 100, 1e-4, 1, Theta);
    std::cout << "Theta: \n" << Theta << std::endl;

    return 0;
}
