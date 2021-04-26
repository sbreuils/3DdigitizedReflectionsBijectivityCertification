

#pragma once

namespace gadg{

    // Compute system A
    Eigen::MatrixXi systemA(const gadg::QuaternionScalar<int>& q){
        Eigen::MatrixXi A(4,6);
        int a = q.value.w;
        int b = q.value.x;
        int c = q.value.y;
        int d = q.value.z;

        A <<  b, c,  d, -b, -c, -d,
                a, -d, c, -a, -d,  c,
                d, a, -b,  d, -a, -b,
                -c, b,  a, -c,  b,  -a;

        return A;
    }

    // compute the Hermite normal form of A with an equivalent algo.
    Eigen::MatrixXi HermiteNormalForm(Eigen::MatrixXi A){
        // hermite normal form of 
        Eigen::MatrixXi H
        Eigen::FullPivLU<Eigen::MatrixXd> lu(A.cast<double>());
        std::cout << "Here is, up to permutations, its LU decomposition matrix:"
             << std::endl << lu.matrixLU() << std::endl;
        return
    }




}

