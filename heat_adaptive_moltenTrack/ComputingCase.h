//
// Created by zyuhao on 6/6/22.
// Class that store the information of
// and function that calculates temperature of one point
//

#ifndef THERMAL_MULTI_MPI_COMPUTINGCASE_H
#define THERMAL_MULTI_MPI_COMPUTINGCASE_H
#include <iostream>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <omp.h>
#include <vector>
#include <set>
#include "MaterialClass.h"
#include "MultiTrackLaserClass.h"
#include "MoltenPoolPoint.h"

using namespace std;

class ComputingCase {

  public:

    // current time
    double timeCurrent;


    // nd form of current time
    double timeCurrent_nd;


    // pointer to material object and laser object
    MaterialClass* materialObject;
    MultiTrackLaserClass* trackLaserObject;


    // list of abssicas of non-dimensional conduction time and weight
    // for gaussian integration
    vector<double> abssicaList;
    vector<double> weightList;


    // list of corresponding laser position w.r.t abssicaList
    vector<double> laserPosXList;
    vector<double>  laserPosYList;


    // list of three inner constants of each node
    // shared at the same time step by all POI
    vector<double> innerConst1List;
    vector<double> innerConst2List;
    vector<double> innerConst3List;


    // physical constant regarding material and laser
    double phyConstMaterialLaser;


    // physical constant used for nondimensionalization of time
    double phyConstTime;


    // physical constant used for nondimensionalization of speed
    double phyConstSpeed;


    // ambient temperature used for initialize temperature at one point
    double tempAmbient;


    // ctor
    ComputingCase(
         MultiTrackLaserClass* inTrackObject,
         MaterialClass* inMaterialObject
         );


    // ~ctor
    ~ComputingCase() = default;


    // member function that set the value timeCurrent
    // used to set a new time during loop w.r.t time
    void setTimeCurrent(double inTimeCurrent);


    // member function that determine the list of abssicas, weights,
    // three inner constants and laser positions
    // based on current time and scanning track history before current time
    // before calculating the list of current time, all three lists will be
    // "cleared"
    void calAbssicaAndPositionList();


    // member function which is used by calAbssicaAndPositionList()
    // given the start and end non-dimensional conduction time
    // this function will push a sub vector
    // containing abssicas of nd conduction time into the total list.
    // So does the laser position.
    void calAbssicaAndPositionListOneTrack(
         double timeDiff_nd_start,
         double timeDiff_nd_end
         );


    // calculate the temperature of one point based on given position of POI
    // and calculated AbssicaList shared by all POI at current time
    double calTempOnePoint(
         double inXc,
         double inYc,
         double inZc
         );


    // give the timeCurrent, meshSize
    // based on beam size and laser position
    // add indexs of potential melt point near laser into the molten pool list
    void addPotentialPointNearLaser(
         double inTimeCurrent,
         double inMeshSize,
         vector<MoltenPoolPoint> &potentialMoltenPoolList,
         int xIndex_min,
         int xIndex_max,
         int yIndex_min,
         int yIndex_max
         );

};


#endif //THERMAL_MULTI_MPI_COMPUTINGCASE_H
