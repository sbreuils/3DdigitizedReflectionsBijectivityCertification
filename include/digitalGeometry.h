#pragma once

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
            std::get<2>(cell)*c_cast,std::get<3>(cell)*c_cast,
            std::get<4>(cell)*c_cast,std::get<5>(cell)*c_cast,
            std::get<6>(cell)*c_cast,std::get<7>(cell)*c_cast
    );
}

// compute the transformation of a region
gadg::IntegerRegion regionRounding(gadg::Region cell) {
    //

    return std::make_tuple(
            gadg::toLipschitzQuaternion(std::get<0>(cell)),gadg::toLipschitzQuaternion(std::get<1>(cell)),
            gadg::toLipschitzQuaternion(std::get<2>(cell)),gadg::toLipschitzQuaternion(std::get<3>(cell)),
            gadg::toLipschitzQuaternion(std::get<4>(cell)),gadg::toLipschitzQuaternion(std::get<5>(cell)),
            gadg::toLipschitzQuaternion(std::get<6>(cell)),gadg::toLipschitzQuaternion(std::get<7>(cell))
    );
}


