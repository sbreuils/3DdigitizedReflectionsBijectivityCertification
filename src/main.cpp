
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
#include "digitalGeometry.h"











int main(int argc, char * argv[]){

    gadg::QuaternionScalar<int> q(2,6,6,9);

    std::cout <<"input q="<<std::endl;
    q.Print();



    Eigen::MatrixXi A = gadg::systemA(q);
    std::cout << "A="<<A<<std::endl;

    // compute the resolution of the system through the hermite normal form of A
    gadg::HermiteNormalForm(A);

    std::cout << "gcd of q = " << gadg::gcd((int)q.value.w,(int)q.value.x,(int)q.value.y,(int)q.value.z)<<std::endl;
    std::cout <<"reduced quaternion = "<<std::endl;
    q = q/gadg::gcd((int)q.value.w,(int)q.value.x,(int)q.value.y,(int)q.value.z);
    q.Print();


	return 0;
}