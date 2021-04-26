
#include <iostream>
#include <algorithm>
#include <tuple>
#include <vector>
#include <cmath> // for approximation
#include <tuple>
#include <Eigen/Dense>

#include "dualQuaternion.h"
#include "util.h"
#include "linearAlgebra.h"
#include "digitalGeometry.h"






int main(int argc, char * argv[]){


//    Eigen::MatrixXi A = gadg::systemA(q);
//    std::cout << "A="<<A<<std::endl;

    // compute the resolution of the system through the hermite normal form of A
    Eigen::MatrixXi A1(3,3);
    A1 << 2,4,4,
         -6,6,12,
          10,4,16; // example from SNF wikipedia page

    std::cout << "A1 = "<<A1 << std::endl;
    gadg::SmithDecomposition snf1 = gadg::SmithDecompositionComputation(A1);
    std::cout << "Smith M1 = \n"<<std::get<0>(snf1) << std::endl;

    // example from https://core.ac.uk/download/pdf/82343294.pdf
    Eigen::MatrixXi A2(5,7);
    A2 <<   1,2,3,4,5,6,7,
            1,0,1,0,1,0,1,
            2,4,5,6,1,1,1,
            1,4,2,5,2,0,0,
            0,0,1,1,2,2,3;

    std::cout << "A2 = "<<A2 << std::endl;
    gadg::SmithDecomposition snf2 = gadg::SmithDecompositionComputation(A2);
    std::cout << "Smith M2 = \n"<<std::get<0>(snf2) << std::endl;
    std::cout << "Smith U2 = \n"<<std::get<1>(snf2) << std::endl;
    std::cout << "Smith V2 = \n"<<std::get<2>(snf2) << std::endl;

    Eigen::VectorXi b(5) ;
    b<<28,4,20,14,9;

    int rank = gadg::getRankSNF(snf2);
    gadg::solve(snf2,b,rank);
    std::cout << "obtained rank ="<<rank << std::endl;

    Eigen::MatrixXf Dinv = Eigen::MatrixXf::Zero(rank,rank);
    for(int i=0;i<rank;++i){
        Dinv(i,i) = 1.0f/(float)(std::get<0>(snf2)(i,i));
    }

    std::cout << "has integer solution -> "<<gadg::hasIntegerSolutions(Dinv,std::get<1>(snf2), b, rank )<<std::endl;

    std::cout << "**************************" << std::endl;

    std::cout << "Compute the set of integer components quaternions inside a region ..."<< std::endl;

    std::cout << "Choose a Lipshitz quaternion ..." << std::endl;
    gadg::QuaternionScalar<int> q(1,7,5,2);
    std::cout <<"input q="<<std::endl;
    q.Print();

    std::cout << "Compute the transformed region ..." << std::endl;
    gadg::IntegerRegion transformedRegion = regionRounding(QuaternionRegionMultiplication(gadg::CZero_Region(), q));

    std::cout << "Compute the bottom right corner" << std::endl;
    gadg::QuaternionScalar<int> bottomCorner = getBottomLeftCornerFromIntegerRegion(transformedRegion);

    std::cout << "Compute the set of Lipschitz quaternions inside integer regions ..." << std::endl;

    std::vector<gadg::QuaternionScalar<int>> setOfLipschitzQ = integerComponentsInsideRegion(transformedRegion,bottomCorner);




    return 0;
}
