//
// Created by zyuhao on 6/4/22.
// This class is inherited from the LaserClass.
// This class contains information describing a multi straight lines scan pattern
// and functions that builds and plots this scan pattern
//
// !!Be careful, a multi line scanning tracks should be constructed using
// function add1LaserOffIntervalAnd1LineTrack chronologically! Because I
// did not use any sorted data structure in this function.
//

#ifndef COMPARE_ANALYTICAL_FVM_MULTITRACKLASERCLASS_H
#define COMPARE_ANALYTICAL_FVM_MULTITRACKLASERCLASS_H

#include "LaserClass.h"
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>

using namespace std;

class MultiTrackLaserClass: public LaserClass{
  public:
    int numLineTrack = 0; // number of line track
    vector<double> timeStartList; // list containing starting time of each track
    vector<double> timeEndList; // list containing ending time of each track
    vector<double> posStartXList; // list containing starting positions along X
    vector<double> posStartYList; // list containing starting positions along Y
    vector<double> posEndXList; // list containing ending positions along X
    vector<double> posEndYList; // list containing ending positions along Y
    vector<double> speedList; // list containing average speed of each segment


    MultiTrackLaserClass() = default; // default ctor

    ~MultiTrackLaserClass() = default; // default dtor

    // update the position of laser position based on the given time
    // return false if the laser is off
    // return true if the laser is still on
    bool calLaserPositionMultiLine(
         double &laserPositionX,
         double &laserPositionY,
         double timeVal
         );

    // update the speed of the laser based on the given time
    // return false if the laser is off
    // return true if the laser is still on
    bool calLaserSpeedMultiLine(
         double &speed,
         double timeVal
         );

    // add a laser off interval and a line segment
    // to total multi line scan pattern
    void add1LaserOffIntervalAnd1LineTrack(
         double inLaserOffInterval,
         double inPosStartX,
         double inPosStartY,
         double inPosEndX,
         double inPosEndY,
         double speed
         );

    // write the scan track data into the file
    void writeScanTrackInto(
         string fileName
         );

};


#endif //COMPARE_ANALYTICAL_FVM_MULTITRACKLASERCLASS_H
