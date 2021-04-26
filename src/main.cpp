
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

    gadg::QuaternionScalar<int> q(2,6,6,9);

    std::cout <<"input q="<<std::endl;
    q.Print();



    Eigen::MatrixXi A = gadg::systemA(q);
    std::cout << "A="<<A<<std::endl;

    // compute the resolution of the system through the hermite normal form of A
    Eigen::MatrixXi Aa(3,4);
    Aa << 2,3,4 ,8,
         -6,6,0 ,0,
         10,0,16,2;
    std::cout <<"Identity of size 4,4 = \n"<<gadg::Identity(4)<<std::endl;

    std::cout << "Aa = "<<Aa << std::endl;

//    std::cout << "implementing SNF functions ..."<<std::endl;
//
    gadg::MatrixIndex Aijmin = gadg::minAij(Aa,2);
    std::cout << "After function, Aa =\n"<<Aa <<std::endl;
    std::cout << "After function, ij min = ("<<Aijmin.first <<","<<Aijmin.second<<")" <<std::endl;


//    std::cout <<"test swap of cols and rows ..."<<std::endl;
//
//    gadg::swapColumns(Aa,0,1); // Ok
//    gadg::swapRows(Aa,1,2);// Ok
//
//    gadg::changeSignRow(Aa, 1);
//    gadg::addToRow(Aa,1,4,0);
//    gadg::addToColumn(Aa,0,-2,1);


    std::cout << "matrix Aa after swap =\n"<<Aa<<"\n"<<std::endl;

    bool isLone = gadg::isLone(Aa,0);

    std::cout <<"is lone = "<<isLone<<std::endl;

//
//    std::cout << "gcd of q = " << gadg::gcd((int)q.value.w,(int)q.value.x,(int)q.value.y,(int)q.value.z)<<std::endl;
//    std::cout <<"reduced quaternion = "<<std::endl;
//    q = q/gadg::gcd((int)q.value.w,(int)q.value.x,(int)q.value.y,(int)q.value.z);
//    q.Print();


	return 0;
}