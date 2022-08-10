cmake_minimum_required(VERSION 3.21)
project(CA_Heat C)

set(CMAKE_C_STANDARD 17)

add_executable(CA_Heat main.c CA_heat.h CA_functions.c heat_functions.c common_function.c guassian_integration.c)

set(EXECUTABLE_OUTPUT_PATH  "/media/zyuhao/UBUNTU 20_0/code/CA_Heat/bin")

find_package(GSL REQUIRED)
target_link_libraries(CA_Heat GSL::gsl GSL::gslcblas -lm)