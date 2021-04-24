
#include <iostream>
#include <algorithm>
#include <tuple>
#include <vector>
#include <cmath> // for approximation
#include <tuple>
#include <Eigen/Dense>

#include "DualQuaternion.h"
#include "util.h"
#include "LinearAlgebra.h"



//Eigen::MatrixXf QuaternionScalarRegionMultiplication(const QuaternionScalar<float>& q,const QuaternionScalar<float>& c )

// compute the Hermite normal form of A with an equivalent algo. than
Eigen::MatrixXi HermiteNormalForm(Eigen::MatrixXi A){

}



// compute the origin cubid voxel C[0]
gadg::Region CZero_Region() {
    return std::make_tuple(
            gadg::QuaternionScalar<float>(-0.5,-0.5,-0.5,0.0),gadg::QuaternionScalar<float>(-0.5,0.5,-0.5,0.0),
            gadg::QuaternionScalar<float>(0.5,-0.5,-0.5,0.0),gadg::QuaternionScalar<float>(0.5,0.5,-0.5,0.0),
            gadg::QuaternionScalar<float>(-0.5,-0.5,0.5,0.0),gadg::QuaternionScalar<float>(-0.5,0.5,0.5,0.0),
            gadg::QuaternionScalar<float>(0.5,-0.5,0.5,0.0),gadg::QuaternionScalar<float>(0.5,0.5,0.5,0.0)
            );
}

// compute the transformation of a region
gadg::Region QuaternionRegionMultiplication(gadg::Region cell, gadg::QuaternionScalar<int>& c) {
    //
    gadg::QuaternionScalar<float> c_cast((float)c.value.x,(float)c.value.y,(float)c.value.z,(float)c.value.w);
    return std::make_tuple(
            std::get<0>(cell)*c_cast,std::get<1>(cell)*c_cast,
            std::get<0>(cell)*c_cast,std::get<1>(cell)*c_cast,
            std::get<0>(cell)*c_cast,std::get<1>(cell)*c_cast,
            std::get<6>(cell)*c_cast,std::get<7>(cell)*c_cast
    );
}




int main(int argc, char * argv[]){

    gadg::QuaternionScalar<int> q(2,6,6,9);

    std::cout <<"q="<<std::endl;
    q.Print();



    Eigen::MatrixXi A = gadg::systemA(q);

    std::cout << "A="<<A<<std::endl;

    std::cout << "gcd of q = " << gadg::gcd((int)q.value.w,(int)q.value.x,(int)q.value.y,(int)q.value.z)<<std::endl;

    std::cout <<"reduced quaternion = "<<std::endl;
    q = q/gadg::gcd((int)q.value.w,(int)q.value.x,(int)q.value.y,(int)q.value.z);
    q.Print();

	return 0;
}