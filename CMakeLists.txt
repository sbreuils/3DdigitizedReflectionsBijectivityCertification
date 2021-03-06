cmake_minimum_required(VERSION 3.5)

project(gaDigitalTransformationCertification)


# include eigen
find_package(Eigen3 REQUIRED)
if(${EIGEN3_FOUND})
    message(STATUS "lib EIGEN3 found")
    message(STATUS "  version " ${EIGEN3_VERSION_STRING})
    message(STATUS "  include " ${EIGEN3_INCLUDE_DIR})
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)




# compilation flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -O2 -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-g")

# files to compile
file(GLOB_RECURSE source_files src/*.cpp include/*.h)
message(STATUS "  files " ${source_files})

# compile
add_executable(gaDigitalTransformationCertification ${source_files})

# link libs
include_directories(${EIGEN3_INCLUDE_DIR}
gaLib/
include/)

