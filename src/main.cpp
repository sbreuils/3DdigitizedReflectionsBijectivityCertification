
#include <iostream>
#include <algorithm>
#include <tuple>
#include <vector>
#include <cmath> // for approximation
#include <tuple>
#include <Eigen/Dense>

#include "DualQuaternion.h"
#include "util.h"

// Compute system A
Eigen::MatrixXi SystemA(const QuaternionScalar<int>& q){
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

//Eigen::MatrixXf QuaternionScalarRegionMultiplication(const QuaternionScalar<float>& q,const QuaternionScalar<float>& c )

// compute the Hermite normal form of A with an equivalent algo. than
Eigen::MatrixXi HermiteNormalForm(Eigen::MatrixXi A){

}

// compute the transformation of a region
QuaternionScalar<int> QuaternionRegionMultiplication(QuaternionScalar<int>& v, QuaternionScalar<int>& c) {
    return (v*c);
}


// compute the origin cubid voxel C[0]
Region CZero_Region() {

}




int main(int argc, char * argv[]){

    QuaternionScalar<int> q(2,6,6,9);

    std::cout <<"q="<<std::endl;
    q.Print();



    Eigen::MatrixXi A = SystemA(q);

    std::cout << "A="<<A<<std::endl;

    std::cout << "gcd of q = " << gcd((int)q.value.w,(int)q.value.x,(int)q.value.y,(int)q.value.z)<<std::endl;

    std::cout <<"reduced quaternion = "<<std::endl;
    q = q/gcd((int)q.value.w,(int)q.value.x,(int)q.value.y,(int)q.value.z);
    q.Print();

	return 0;
}