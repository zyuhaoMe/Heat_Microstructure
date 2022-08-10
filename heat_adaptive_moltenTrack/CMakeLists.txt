cmake_minimum_required(VERSION 3.20)
project (thermal_multi_mpi)

set(CMAKE_CXX_STANDARD 14)
set(EXECUTABLE_OUTPUT_PATH "/home/zyuhao/Desktop/thermal_multi_mpi")

# Add MPI Package to Project
find_package(MPI REQUIRED)

# Add openmp
find_package(OpenMP REQUIRED)

# Add gsl
find_package(GSL REQUIRED)


add_executable(
        thermal_multi_mpi
        main.cpp constants.h
        LaserClass.cpp LaserClass.h
        MaterialClass.cpp MaterialClass.h
        MultiTrackLaserClass.cpp MultiTrackLaserClass.h
        ComputingCase.cpp ComputingCase.h
        ThermalHistOfCube.cpp ThermalHistOfCube.h
        ThermalHistOfCubeMoltenPoolTrack.cpp ThermalHistOfCubeMoltenPoolTrack.h
        MoltenPoolPoint.cpp MoltenPoolPoint.h
        MicroStructureOfCube.cpp MicroStructureOfCube.h)

# Add libraries for code completion and compiling
target_link_libraries(thermal_multi_mpi PUBLIC MPI::MPI_CXX )
if(OpenMP_CXX_FOUND)
    target_link_libraries(thermal_multi_mpi PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(thermal_multi_mpi PUBLIC GSL::gsl GSL::gslcblas)