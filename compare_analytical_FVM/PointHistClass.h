/* 
 * Author: zyuhao 
 * Date: 26/04/22
 * Content:  head file of the PointHistClass
 * contains the variables and functions of creating, writing
 * the thermal history of one points in the cube
 */

//

#ifndef COMPARE_ANALYTICAL_FVM_POINTHISTCLASS_H
#define COMPARE_ANALYTICAL_FVM_POINTHISTCLASS_H

#include "CubeClass.h"
#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

class PointHistClass {
  public:
    int numPlot; // number of time steps to record the temperature of one point
    double* timeArray; // 1d array to store the time of each time step
    double* tempArray; // 1d array to store the temp of each time step


    int xSize, ySize, zSize; // size of index along 3 directions
    double pointPosX; // position of the point
    double pointPosY;
    double pointPosZ;

    int pointIndexX; // calculate position of the point
    int pointIndexY;
    int pointIndexZ;

    double tempLimit; // the upper limit of temperature to be stored

    // function that initialize the object
    // given the point 3d position (w.r.t the local corner of the cube,
    // not the starting position of the cube)
    // and cube information, calculate the 3d index of the position
    // given the numPlot, dynamically allocate two arrarys
    void initializePointHist(
         double inPointPosX,
         double inPointPosY,
         double inPointPosZ,
         const CubeClass* inMyCube,
         int numPlot,
         double inTempLimit
         );

    // function that store time and temperature of one time step into 2 arrays
    // given the timeStepCurrent, which has already been calculated in
    // writeTemp2D function
    // exit the program if index out of bound
    // this function will be called within the "writeTemp2D" function
    void storePointOneTime(
         int timeStepCurrent,
         double timeCurrent,
         double* temp);

    // function that check and check the point

    // function that write the total array of time and temperature
    // which is thermal history
    // into a dat file with a name fileName
    void writePointHist(string fileName);
};


#endif //COMPARE_ANALYTICAL_FVM_POINTHISTCLASS_H
