
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






#define TESTED 1




int main(int argc, char * argv[]){

    gadg::QuaternionScalar<int> q(1,3,5,2);

    std::cout <<"input q="<<std::endl;
    q.Print();



//    Eigen::MatrixXi A = gadg::systemA(q);
//    std::cout << "A="<<A<<std::endl;

    // compute the resolution of the system through the hermite normal form of A
    Eigen::MatrixXi A(3,3);
    A << 2,4,4,
         -6,6,12,
          10,4,16; // example from SNF wikipedia page

    std::cout << "A = "<<A << std::endl;


    Eigen::MatrixXi M = gadg::SmithDecomposition(A);

    std::cout << "Smith M = \n"<<M << std::endl; 

	return 0;
}