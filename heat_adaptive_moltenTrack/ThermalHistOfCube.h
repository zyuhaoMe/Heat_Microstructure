//
// Created by zyuhao on 6/7/22.
// This class calculate the thermal history of a cube
// based upon ComputingCase object

#ifndef THERMAL_MULTI_MPI_THERMALHISTOFCUBE_H
#define THERMAL_MULTI_MPI_THERMALHISTOFCUBE_H
#include <iostream>
#include <fstream>
#include <string>
#include "ComputingCase.h"
#include "MoltenPoolPoint.h"


class ThermalHistOfCube {

  public:
    double xLen, yLen, zLen; // length of the cube
    double meshSize; // size of mesh
    int xsize, ysize, zsize; // number of mesh along x,y,z direction
    int numPlot; // number of plot
    double timeLen; // length of the total time

    ComputingCase* computingCaseObject; // object of computing case

    vector<vector<vector<double>>> tempField; // 3d temperature field

    vector<double> timeList; // list of timeCurrent to calculate thermal history
    vector<double> xc; // list of x-axis coordinates
    vector<double> yc; // list of y-axis coordinates
    vector<double> zc; // list of z-axis coordinates


    // default ctor
    ThermalHistOfCube();


    // ctor
    ThermalHistOfCube(
         double inXLen,
         double inYLen,
         double inZLen,
         double inMeshSize, // size of mesh
         double inTimeLen, // time length
         int inNumPlot,
         ComputingCase* inComputingCase
         );

    // dtor
    ~ThermalHistOfCube() = default;

    // function that calculate thermal field
    // based on 1d array: timeList
    // store the thermal field with a give name
    void calThermalHist(string fileName);

    // function that calculate the thermal field at a time
    // this function will be called at each time step by calThermalHist function
    void calThermalFieldAtOneTime(double timeCurrent);


    // function that store thermal field at each time step
    void storeThermalField(
         string fileName,
         double solutionTime,
         vector<vector<vector<double>>> inTempField
         );


};
#endif //THERMAL_MULTI_MPI_THERMALHISTOFCUBE_H
