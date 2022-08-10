cmake_minimum_required(VERSION 3.21)
project(compare_analytical_FVM)

set(CMAKE_CXX_STANDARD 14)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

set(EXECUTABLE_OUTPUT_PATH "/media/zyuhao/Media/code/compare_analytical_FVM")



add_executable(compare_analytical_FVM main.cpp MaterialClass.cpp MaterialClass.h LaserClass.cpp LaserClass.h CubeClass.cpp CubeClass.h TimeClass.cpp TimeClass.h TempClass.cpp TempClass.h constants.h PointHistClass.cpp PointHistClass.h PointGroupHistClass.cpp PointGroupHistClass.h MultiTrackLaserClass.cpp MultiTrackLaserClass.h)

# find gsl library (used in numerical integration) for linux
find_package(GSL REQUIRED)
target_link_libraries(compare_analytical_FVM GSL::gsl GSL::gslcblas -lm)

## find gnuplot
#find_package(Boost REQUIRED)
#if(BOOST_FOUND)
#    message(STATUS "find boost")
#endif()
#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/externals/gnuplot-iostream")

#target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

# find gsl library (used in numerical integration) for windows
