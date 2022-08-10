/* 
 * Author: zyuhao 
 * Date: 2022-04-09
 * Content:  head file of TimeClass,
 * contain total time, time step and functions of determining time step in
 * different ways
 */

//

#ifndef COMPARE_ANALYTICAL_FVM_TIMECLASS_H
#define COMPARE_ANALYTICAL_FVM_TIMECLASS_H

#include <iostream>
#include <cmath>
#include "CubeClass.h"
#include "MaterialClass.h"

using namespace std;

class TimeClass {
  public:
    double totalTime; // total time of simulation
    double timeStep; // size of time step
    int numTimeStep; // number of time steps, should be an integer

    CubeClass* cubePtr; // will use the meshSize defined in the cube geometry
    MaterialClass* materialPtr;

    // default ctor
    TimeClass(CubeClass *myCube, MaterialClass* myMaterial);


    // set the time step using a self-defined timeStep
    // return true if the timeStep is proper
    bool setTimeStep(double myTimeStep);

    // set the time step based on the Von Neumann stability
    // this is specifically for FVM solver
    // the result depends on the mesh size and thermal diffusivity
    // a safety coefficient between 0 and 1 can be used to tune timeStep
    // the smaller, more accuracy but slower
    void setTimeStepWithVonNeu(double safeCoeff);

    // set the time step based on the numPlots and total time
    // this is specifically for analytical solver
    // the analytical method requires significantly more time to calculate
    // thermal history, but can give prediction of temp field at a specific
    // time independently, therefore, to save time, we had better just calculate
    // the temp field os specific time step exactly when it is plotted
    void setTimeStepWithExactPlot(int inNumPlot);


    // function that set the totalTime
    void setTotalTime(double myTimeTotal);

    // function that recalculate the numTimeStep after resetting the totalTime
    // or timeStep
    void calNumTimeStep();


    // friend function for <<
    friend ostream& operator<<(ostream& os, const TimeClass& myTimePlan);
};

// overloading << for TimeClass
ostream& operator<<(ostream& os, const TimeClass& myTimePlan);

#endif //COMPARE_ANALYTICAL_FVM_TIMECLASS_H
