/* 
 * Author: zyuhao 
 * Date: 27/04/22
 * Content:  head file of PointGroupClass
 * create a group of points to store thermal history
 * iteratively call the function of each PointHistClass object
 * to store and write
 * 
 */

//

#ifndef COMPARE_ANALYTICAL_FVM_POINTGROUPHISTCLASS_H
#define COMPARE_ANALYTICAL_FVM_POINTGROUPHISTCLASS_H

#include "CubeClass.h"
#include "PointHistClass.h"
#include <string>

using namespace std;

class PointGroupHistClass {
  public:
    int numPoint; // number of points
    double *pointPosArrayX;// array of positions of points
    double *pointPosArrayY;
    double *pointPosArrayZ;

    PointHistClass *pointArray; // array of the PointHistClass object
    string *fileNameArray; // array of file name to write into

    // ctor, give start point and end point of a straight line
    // w.r.t the corner of the cube and the number of points
    // and the fileName of TempClass who call it.
    // First, this will generate uniformly-distributed points along that line
    // including starting and end points
    // then it will initialize all the fileName array into "to be given"
    // after that it will define all the variables for each point object.
    // add more functions if necessary
    PointGroupHistClass(
         double inStartPosX, // start position of the line
         double inStartPosY,
         double inStartPosZ,
         double inEndPosX, // end position of the line
         double inEndPosY,
         double inEndPosZ,
         int inNumPoint, // number of uniform point along the line
         int inNumPlot, // number of time steps to store the temperature
         double inTempLimit, // upper limit of temperature
         const CubeClass *inMyCube // myCube object
         );

    // function that change the name of file names based on the fileName of
    // TempClass, which is related to other conditions s.t. the numerical solver
    void setFileName(const string tempObjectFileName);


    // function that store time and temp at one time step for all point objects
    void storePointsOneTime(
         int timeStepCurrent,
         double timeCurrent,
         double* temp
         );


    // function that write thermal hist of all point objects into dat files
    void writePointsHist();

    // friend function of ostream
    friend ostream& operator<<(
         ostream& os,
         const PointGroupHistClass& myPointGroup
         );
};

// overloading << for PointGroupHistClass
ostream& operator<<(
     ostream& os,
     const PointGroupHistClass& myPointGroup
     );

#endif //COMPARE_ANALYTICAL_FVM_POINTGROUPHISTCLASS_H
