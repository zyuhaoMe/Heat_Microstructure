//
// Created by zyuhao on 6/8/22.
// class that track molten pool or any temperature contour.
// First, the whole temperature field at starting time will be calculated.
// Then, thermal history above specific temperature will be
// calculated using molten pool track algorithm.
// This algorithm is used to minimize computations by calculating only the
// data of interest.
//

#ifndef THERMAL_MULTI_MPI_THERMALHISTOFCUBEMOLTENPOOLTRACK_H
#define THERMAL_MULTI_MPI_THERMALHISTOFCUBEMOLTENPOOLTRACK_H

#include "ThermalHistOfCube.h"
#include "MoltenPoolPoint.h"
#include "constants.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <string>

using namespace std;

class ThermalHistOfCubeMoltenPoolTrack: public ThermalHistOfCube{

  public:
    double timeStepSize; // size of time step
    double numTimeStep; // number of time step
    double startingTime; // starting time of calculation


    // temperature limit of molten pool, can be set to other value to
    // track a local temp field downwards to a specific contour.
    double tempBound;


    // vector of molt pool point objects at old time step and new time step
    vector<MoltenPoolPoint> moltenPoolList_old;
    vector<MoltenPoolPoint> moltenPoolList_new;


    // vector of points to be computed
    vector<MoltenPoolPoint> listToCompute;


    // 3d bool arrays indicating whether points have been updated
    vector<vector<vector<bool>>> tempHasComputed;


    // reconstructed 3d tempField data containing temperature above tempLimit
    // and with ambient temp in other points
    vector<vector<vector<double>>> tempToPlot;


    // ctor
    ThermalHistOfCubeMoltenPoolTrack(
         double inXLen,
         double inYLen,
         double inZLen,
         double inMeshSize,
         double inStartingTime,
         double inTimeLen,
         double inTimeStepSize,
         double inTempBond,
         int inNumPlot,
         ComputingCase* inComputingCase
         );


    // dtor
    ~ThermalHistOfCubeMoltenPoolTrack() = default ;


    // member function that detect molten pool points from a full 3d temperature
    // field and add those points objects into the moltenPoolList_new
    void detectAndAddMoltenPoolFromTempField();


    // member function track molten pool at a new time step based on the
    // molten pool information from the last time step
    void trackMoltenPoolAtNewTime(double timeCurrent);


    // member function that clear moltenPoolList_new and fill moltenPoolList_new
    // with moltenPoolList_old
    void clearAndSwapMoltenPoolList();


    // member functions that track molten pool.
    // It starts from total temperature field starting time,
    // increase time and track molten pool neighbour by neighbour until
    // all the possible molten pool has been tracked and added into the list
    // Based on the number of plots, write data into temperature field.
    void trackMoltenPoolUntilEnd(string fileName);


    // member functions that transform the molten pool vector and
    // the corresponding point in tempField into a new 3d tempField
    // and fill the undefined data in the 3d vector with ambient temperature.
    // Only do this transform when writing data into .dat file.
    // will initialize all the temperature field before assigning value
    void transformMoltenPoolDataIntoTempField();

};


#endif //THERMAL_MULTI_MPI_THERMALHISTOFCUBEMOLTENPOOLTRACK_H
