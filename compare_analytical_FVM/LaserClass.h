/* 
 * Author: zyuhao 
 * Date: 2022-04-08
 * Content: head file of the LaserClass
 * contain the physic parameters of the laser
 * and the function describe the laser position wrt current time
 */

//

#ifndef COMPARE_ANALYTICAL_FVM_LASERCLASS_H
#define COMPARE_ANALYTICAL_FVM_LASERCLASS_H

#include "constants.h"
#include <cmath>
#include <iostream>
using namespace std;

class LaserClass {
  public:
    double laserPower; // laser power
    double beamRadiusX; // the 3d-dimension of heat source
    double beamRadiusY;
    double beamRadiusZ;
    double scanTime; // time when the laser is active

    // index of scanning path (stationary, linear, spiral)
    int scanPatternIndex;

    // if laserClass is multi line track, default ctor
    LaserClass()
    {
      laserPower = 280.;
      beamRadiusX = 50e-6;
      beamRadiusY = 50e-6;
      beamRadiusZ = 50e-7;

      scanPatternIndex = INDEX_SCANNING_PATH_MULTI_LINE;
      scanTime = 0.;
    }

    // default ctor, input the index of scanning path
    LaserClass(int inSPIndex);

    // function that set the number of model

    // function calculate the position of the laser
    // based on the given time timeVal
    // return true if laser is still on and update laserPositionX/Y by reference
    // return false if the laser is off
    bool calLaserPosition(
         double &laserPositionX,
         double &laserPositionY,
         double timeVal
         );

    // function calculate the instentaneous speed of the laser beam
    // based on the given time timeVal
    // return true if the laser is still on and update laserSpeed by reference
    // return false if the laser if off
    bool calLaserSpeed(
         double &laserSpeed,
         double timeVal
         );

    // virtual function updating xy positions of laser for multi line scanning
    virtual bool calLaserPositionMultiLine(
         double &laserPositionX,
         double &laserPositionY,
         double timeVal
         );

    // virtual function updating abs speed of laser for multi line scanning
    virtual bool calLaserSpeedMultiLine(
         double &laserSpeed,
         double timeVal
         );


    // friend function of ostream
    friend ostream& operator<<(ostream& os, const LaserClass& myLaser);
};

// overloading << for LaserClass
ostream& operator<<(ostream& os, const LaserClass& myLaser);


#endif //COMPARE_ANALYTICAL_FVM_LASERCLASS_H
