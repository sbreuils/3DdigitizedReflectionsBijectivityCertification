
#include <iostream>
#include <algorithm>
#include <tuple>
#include <vector>
#include <cmath> // for approximation

#include <Eigen/Dense>
#include "DualQuaternion.h"

// Compute system A
Eigen::MatrixXf SystemA(const QuaternionScalar<float>& q){
    Eigen::MatrixXf A(4,6);
    float a = q.value.w;
    float b = q.value.x;
    float c = q.value.y;
    float d = q.value.z;

    A <<  b, c,  d, -b, -c, -d,
          a, -d, c, -a, -d,  c,
          d, a, -b,  d, -a, -b,
         -c, b,  a, -c,  b,  -a;

    return A;
}

Eigen::MatrixXf QuaternionScalarRegionMultiplication(const QuaternionScalar<float>& q,const QuaternionScalar<float>& c )



int main(int argc, char * argv[]){

    QuaternionScalar<float> q(2,3,4,1);
    Eigen::MatrixXf A = SystemA(q);

    std::cout << "A="<<A<<std::endl;




	return 0;
}