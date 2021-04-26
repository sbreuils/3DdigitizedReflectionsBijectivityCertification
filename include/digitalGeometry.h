#pragma once
#include <array>
#include <algorithm>

namespace gadg {

// compute the origin cubid voxel C[0]
    Region CZero_Region() {
        return std::make_tuple(
                gadg::QuaternionScalar<float>(-0.5, -0.5, -0.5, 0.0),
                gadg::QuaternionScalar<float>(-0.5, 0.5, -0.5, 0.0),
                gadg::QuaternionScalar<float>(0.5, -0.5, -0.5, 0.0), gadg::QuaternionScalar<float>(0.5, 0.5, -0.5, 0.0),
                gadg::QuaternionScalar<float>(-0.5, -0.5, 0.5, 0.0), gadg::QuaternionScalar<float>(-0.5, 0.5, 0.5, 0.0),
                gadg::QuaternionScalar<float>(0.5, -0.5, 0.5, 0.0), gadg::QuaternionScalar<float>(0.5, 0.5, 0.5, 0.0)
        );
    }


// compute the transformation of a region
    Region QuaternionRegionMultiplication(gadg::Region cell, gadg::QuaternionScalar<int> &c) {
        //
        gadg::QuaternionScalar<float> c_cast((float) c.value.x, (float) c.value.y, (float) c.value.z,
                                             (float) c.value.w);
        return std::make_tuple(
                std::get<0>(cell) * c_cast, std::get<1>(cell) * c_cast,
                std::get<2>(cell) * c_cast, std::get<3>(cell) * c_cast,
                std::get<4>(cell) * c_cast, std::get<5>(cell) * c_cast,
                std::get<6>(cell) * c_cast, std::get<7>(cell) * c_cast
        );
    }

// compute the transformation of a region
    IntegerRegion regionRounding(gadg::Region cell) {
        return std::make_tuple(
                gadg::toLipschitzQuaternion(std::get<0>(cell)), gadg::toLipschitzQuaternion(std::get<1>(cell)),
                gadg::toLipschitzQuaternion(std::get<2>(cell)), gadg::toLipschitzQuaternion(std::get<3>(cell)),
                gadg::toLipschitzQuaternion(std::get<4>(cell)), gadg::toLipschitzQuaternion(std::get<5>(cell)),
                gadg::toLipschitzQuaternion(std::get<6>(cell)), gadg::toLipschitzQuaternion(std::get<7>(cell))
        );
    }

    //
    gadg::QuaternionScalar<int> getBottomLeftCornerFromIntegerRegion(const IntegerRegion& transformedRegion) {
        std::vector<int> wi{ {std::get<0>(transformedRegion).value.w, std::get<1>(transformedRegion).value.w,
                              std::get<2>(transformedRegion).value.w, std::get<3>(transformedRegion).value.w,
                              std::get<4>(transformedRegion).value.w, std::get<5>(transformedRegion).value.w,
                              std::get<6>(transformedRegion).value.w, std::get<7>(transformedRegion).value.w} };

        std::array<int,8> xi{ {std::get<0>(transformedRegion).value.x, std::get<1>(transformedRegion).value.x,
                                      std::get<2>(transformedRegion).value.x, std::get<3>(transformedRegion).value.x,
                                      std::get<4>(transformedRegion).value.x, std::get<5>(transformedRegion).value.x,
                                      std::get<6>(transformedRegion).value.x, std::get<7>(transformedRegion).value.x} };

        std::array<int,8> yi{ {std::get<0>(transformedRegion).value.y, std::get<1>(transformedRegion).value.y,
                                      std::get<2>(transformedRegion).value.y, std::get<3>(transformedRegion).value.y,
                                      std::get<4>(transformedRegion).value.y, std::get<5>(transformedRegion).value.y,
                                      std::get<6>(transformedRegion).value.y, std::get<7>(transformedRegion).value.y} };

        std::array<int,8> zi{ {std::get<0>(transformedRegion).value.z, std::get<1>(transformedRegion).value.z,
                                      std::get<2>(transformedRegion).value.z, std::get<3>(transformedRegion).value.z,
                                      std::get<4>(transformedRegion).value.z, std::get<5>(transformedRegion).value.z,
                                      std::get<6>(transformedRegion).value.z, std::get<7>(transformedRegion).value.z} };

        int wmin = *std::min_element(wi.begin(),wi.end());
        int xmin = *std::min_element(xi.begin(),xi.end());
        int ymin = *std::min_element(yi.begin(),yi.end());
        int zmin = *std::min_element(zi.begin(),zi.end());

        return {xmin,ymin,zmin,wmin};
    }

// return a vector of quaternions inside a region with integer components
    std::vector<QuaternionScalar<int>> integerComponentsInsideRegion(const IntegerRegion& integerArea, const QuaternionScalar<int>& bottomRightQuatBound){
        std::vector<QuaternionScalar<int>> quaternionsInsideRegion;


        // bound of the lipschitz quaternion
        int wmin = bottomRightQuatBound.value.w;
        int wmax = -wmin;
        int xmin = bottomRightQuatBound.value.x;
        int xmax = -xmin;
        int ymin = bottomRightQuatBound.value.y;
        int ymax = -ymin;
        int zmin = bottomRightQuatBound.value.z;
        int zmax = -zmin;

        for(int w=wmin ; w<wmax ; ++w){
            for(int x=xmin ; x<xmax ; ++x){
                for(int y=ymin ; y<ymax ; ++y){
                    for(int z=zmin ; z<zmax ; ++z){
                        quaternionsInsideRegion.emplace_back(QuaternionScalar<int>(x,y,z,w));
                    }
                }
            }
        }


        std::cout << "Print obtained integer set of points" << std::endl;
        for(QuaternionScalar<int> lipschitzQi : quaternionsInsideRegion){
            lipschitzQi.Print();
        }
        return quaternionsInsideRegion;
    }

}