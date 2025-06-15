#include<iostream>
#include "../includes/matrix.h"
#include<cmath>
#include <cstdlib>




using namespace pp;

int main(){

    matrix A = matrix(15,10);
    std::cout<< "Matrix A (15 x 10)" << std::endl;
    A.randomize();
    A.print();

    //Decompose into Q and R
    QR::mtuple qr = QR::decomp(A);
    matrix Q = std::get<0>(qr); // Q matrix
    matrix R = std::get<1>(qr); // R matrix

    std::cout << "\nMatrix R (should be upper triangular)" << std::endl;
    R.print();

    std::cout << (R.is_upper_triangular() ? "True: Upper triangular": "False: Not upper triangular");

    std::cout << "\n \nMatrix product of Q^T Q (Should be identity)";
    matrix QT = Q.transpose();
    matrix QTQ = QT * Q;
    QTQ.print();
    matrix iden = matrix(QTQ.size1(), QTQ.size2());
    iden.setid();
    std::cout << (pp::approx(QTQ, iden) ? "True: Q^T Q is identity" : "False: Q^T Q is not identity");
    
    std::cout << "\n \nMatrix QR (Should be equal to matrix A)";
    matrix M_QR = Q * R;
    M_QR.print();
    std::cout << (pp::approx(M_QR, A) ? "True: QR is equal to A" : "False: QR is not equal to A\n \n");


    //Check solve function:
    vector x_true(A.size2());
    x_true.randomize();
    vector b = A * x_true; // b = A*x_true
    std::cout << "\nVector b (15 x 1) created from A*x_true (To ensure a solution exists)" << std::endl;
    b.print();
    vector x = QR::solve(Q, R, b);
    std::cout << "\nSolution x (R*x should be equal to Q^T b)" << std::endl;
    x.print();
    vector check = Q.transpose() * b;
    vector Rx = R * x; // R*x
    
    std::cout << "\nCheck Q^T b (should be equal to R*x):" << std::endl;
    std::cout << "\nR*x:" << std::endl;
    Rx.print();
    std::cout << "\nQ^T b:" << std::endl;
    check.print();
    std::cout << (pp::approx(R*x, check) ? "True: R*x is equal to Q^T b" : "False: x is not equal to Q^T b\n \n");
    //Check 2:
    vector check2 = A* x; // A*x
    std::cout << "\nCheck A*x (should be equal to b):" << std::endl;
    std::cout << "\nA*x:" << std::endl;
    check2.print();
    std::cout << "\nb:" << std::endl;
    b.print();
    std::cout << (pp::approx(A*x, b) ? "True: A*x is equal to b" : "False: A*x is not equal to b\n \n");
    return 0;
}