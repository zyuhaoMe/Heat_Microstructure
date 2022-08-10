/*
 * Author: zyuhao
 * Date: 2022-04-15
 *
 * Main function to define the parameters: Material, Laser, Geometry,
 * then the thermal history will be calculated using both analytical solver
 * and FVM solver
 *
 * Then several parameter of cube position and laser scan pattern
 * will be changed to investigate the effect of the boundaries and laser scan
 * pattern on the error between analytical model and FVM model
 *
 * Thermal history of specific time step and specific points
 * will be written as a file for post-processing
 *
 */

#include "LaserClass.h"
#include "MultiTrackLaserClass.h"
#include "MaterialClass.h"
#include "CubeClass.h"
#include "TimeClass.h"
#include "TempClass.h"
#include "PointGroupHistClass.h"
#include "constants.h"
#include <iostream>
#include <string>

using namespace std;

int main() {

  // build a material object
  MaterialClass* steel = new MaterialClass("steel");
  // print material parameter
  cout << *steel << endl;

  // build a material object named IN625
  MaterialClass* IN625 = new MaterialClass("IN625");
  // set the properties
  IN625->setDens(8440);
  IN625->setSpecHeat(500);
  IN625->setThermalCond(16);
  IN625->setTempLiquid(1623);
  IN625->setAbsorp(0.43);
  IN625->setConvecCoeff(10);
  IN625->setEmissivity(0.8);
  // print info
  cout << IN625 << endl;

  // build a material object named Steel155
  MaterialClass* steel155 = new MaterialClass("Steel155");
  // set the properties
  steel155->setDens(7770);
  steel155->setSpecHeat(420);
  steel155->setThermalCond(17.8);
  steel155->setTempLiquid(1700);
  steel155->setAbsorp(0.22);
  steel155->setConvecCoeff(10);
  steel155->setEmissivity(0.8);
  // print info
  cout << steel155 << endl;

  // build a cube (contain meshes) object
  CubeClass* myCube = new CubeClass();
  // print cube parameter
  cout << *myCube << endl;

  // first, we will compute the result of FVM and analytical solution
  // of a stationary heat source
  // for the analytical method, it is no need to set the time step based on
  // Von Noemann stability,
  // therefore, it is better that, FVM Temp and Analy Temp use two different
  // Time objects with same totalTime

  // total time of interested
  double timeTotal_stationary = 2. * TIME_LASER_ACTIVE_STATIONARY;

  // build a laser object of stationary scanning path
  LaserClass* myLaser = new LaserClass(INDEX_SCANNING_PATH_STATIONARY);
  cout << *myLaser << endl;

  int numPlot = 10; // number of plot we are interested in

  // build a point group object
  double startPosX = 1e-3; // center line position
  double startPosY = 0.;
  double startPosZ = 0.;
  double endPosX = 1e-3;
  double endPosY = 2e-3;
  double endPosZ = 0.;
  int numPoint = 10;

  PointGroupHistClass* myPointGroup0 = new PointGroupHistClass(
       startPosX,
       startPosY,
       startPosZ,
       endPosX,
       endPosY,
       endPosZ,
       numPoint,
       numPlot,
       steel->tempLiquid,
       myCube
       );
  cout << *myPointGroup0 << endl;
//
//
//  // stationary heat source FVM at position 0
//  TimeClass* timeStationaryFVM = new TimeClass(myCube, steel155);
//  timeStationaryFVM->setTotalTime(timeTotal_stationary);
//  timeStationaryFVM->setTimeStepWithVonNeu(0.8);
//  cout << *timeStationaryFVM << endl;
//
//  TempClass* caseStationaryFVM0 = new TempClass(
//       steel,
//       myCube,
//       myLaser,
//       timeStationaryFVM,
//       myPointGroup0
//       );
//  caseStationaryFVM0->setNumPlot(numPlot);
//  caseStationaryFVM0->calTempHistFVM("stationaryFVM_0_");
//  caseStationaryFVM0->calTempHistFVMwithConvAndRadiaBC("stationaryFVM3BC_0_");

//
//  // stationary heat source analytical at position 0
//  TimeClass* timeStationaryAna = new TimeClass(myCube, steel);
//  timeStationaryAna->setTotalTime(timeTotal_stationary);
//  timeStationaryAna->setTimeStepWithExactPlot(numPlot);
//  cout << *timeStationaryAna << endl;
//
//  TempClass* caseStationaryAna0 = new TempClass(
//       steel,
//       myCube,
//       myLaser,
//       timeStationaryAna,
//       myPointGroup0
//       );
//  caseStationaryAna0->setNumPlot(numPlot);
//  //caseStationaryAna0->calTempHistAnalytical("stationaryAna_0_");]


//
////
//  // stationary heat source FVM at position 1
//  // reuse the time of FVM
//  // move the cube to make stationary point close to boundary
//  myCube->setStartPositionXYZ(0.5e-3, 0.5e-3, 0);
//  cout << *myCube << endl;
//  // new group point
//  // build a point group object
//  startPosX = 0.5e-3; // center line position
//  startPosY = 0.;
//  startPosZ = 0.;
//  endPosX = 0.5e-3;
//  endPosY = 2e-3;
//  endPosZ = 0.;
//  numPoint = 10;
//  PointGroupHistClass* myPointGroup1 = new PointGroupHistClass(
//       startPosX,
//       startPosY,
//       startPosZ,
//       endPosX,
//       endPosY,
//       endPosZ,
//       numPoint,
//       numPlot,
//       steel->tempLiquid,
//       myCube
//       );
//  cout << *myPointGroup0 << endl;
//  TempClass* caseStationaryFVM1 = new TempClass(
//       steel155,
//       myCube,
//       myLaser,
//       timeStationaryFVM,
//       myPointGroup1
//       );
//  caseStationaryFVM1->setNumPlot(numPlot);
//  caseStationaryFVM1->calTempHistFVM("stationaryFVM_1_");
//  caseStationaryFVM1->calTempHistFVMwithConvAndRadiaBC("stationaryFVM3BC_1_");

//
//
//  // stationary heat source analytical at position 1
//  // reuse the cube of FVM1
//  // reuse the time of analytical
//  TempClass* caseStationaryAna1 = new TempClass(
//       steel,
//       myCube,
//       myLaser,
//       timeStationaryAna,
//       myPointGroup1
//       );
//  caseStationaryAna1->setNumPlot(numPlot);
//  caseStationaryAna1->calTempHistAnalyMirror("stationaryAnaMirror_1_");


  // linear heat source FVM at the center
  // set the linear scan plan
  LaserClass* myLaserLinear = new LaserClass(INDEX_SCANNING_PATH_LINEAR);
  cout << *myLaserLinear << endl;
  // reset the cube position
  myCube->setStartPositionXYZ(0, 0, 0);
  // set the scan time
  TimeClass* timeLinearFVM = new TimeClass(myCube, steel);
  timeLinearFVM->setTotalTime(2e-3);
  timeLinearFVM->setTimeStepWithVonNeu(0.4);
  cout << *timeLinearFVM << endl;
  // reset the group point
  startPosX = 0.; // center line position
  startPosY = 1e-3;
  startPosZ = 0.;
  endPosX = 2e-3;
  endPosY = 1e-3;
  endPosZ = 0.;
  numPoint = 10;
  PointGroupHistClass* myPointGroup2 = new PointGroupHistClass(
       startPosX,
       startPosY,
       startPosZ,
       endPosX,
       endPosY,
       endPosZ,
       numPoint,
       numPlot,
       steel->tempLiquid,
       myCube
       );
  TempClass* caseLinearFVM0 = new TempClass(
       steel,
       myCube,
       myLaserLinear,
       timeLinearFVM,
       myPointGroup2
       );
  caseLinearFVM0->setNumPlot(numPlot);
  //caseLinearFVM0->calTempHistFVM("linearFVM_0_");


  // linear heat source analytical at the center
  // time plan for linear analytical
  TimeClass* timeLinearAna = new TimeClass(myCube, steel);
  timeLinearAna->setTotalTime(2e-3);
  timeLinearAna->setTimeStepWithExactPlot(numPlot);
  cout << *timeLinearAna << endl;

  TempClass* caseLinearAna0 = new TempClass(
       steel,
       myCube,
       myLaserLinear,
       timeLinearAna,
       myPointGroup2
       );
  caseLinearAna0->setNumPlot(numPlot);
  //caseLinearAna0->calTempHistAnalyticalAdaptive("linearAnaAdap_0_");

  // multi line scan track
  MultiTrackLaserClass *multiTrackLaser1  = new MultiTrackLaserClass();
  multiTrackLaser1->add1LaserOffIntervalAnd1LineTrack(
       0,
       0.5e-3, 0.5e-3, 1.5e-3, 0.5e-3,
       1
  );
  multiTrackLaser1->add1LaserOffIntervalAnd1LineTrack(
       1e-3,
       1.5e-3, 1e-3, 0.5e-3 , 1e-3,
       1
  );
  multiTrackLaser1->add1LaserOffIntervalAnd1LineTrack(
       1e-3,
       0.5e-3, 1.5e-3, 1.5e-3, 1.5e-3,
       1
  );

  multiTrackLaser1->writeScanTrackInto("test_data");

  TimeClass* timeMultiLineAna = new TimeClass(myCube, steel);
  timeMultiLineAna->setTotalTime(6e-3);
  int numPlotMultiLine = 30;
  timeMultiLineAna->setTimeStepWithExactPlot(numPlotMultiLine);

  TempClass* caseMultiLineAna0 = new TempClass(
       steel,
       myCube,
       multiTrackLaser1,
       timeMultiLineAna,
       myPointGroup2
       );

  caseMultiLineAna0->setNumPlot(numPlotMultiLine);
  caseMultiLineAna0->calTempHistAnalyticalAdaptive("multiLinearAnaAdap_0_");


//
//
//  // linear heat source near the boundary
//  // set the linear scan plan
//  cout << *myLaserLinear << endl;
//  // reset the cube position
//  myCube->setStartPositionXYZ(0, 0.8e-3, 0);
//
//  cout << *timeLinearFVM << endl;
//  // reset the group point
//  startPosX = 0.; // center line position
//  startPosY = 0.2e-3;
//  startPosZ = 0.;
//  endPosX = 2e-3;
//  endPosY = 0.2e-3;
//  endPosZ = 0.;
//  numPoint = 10;
//  PointGroupHistClass* myPointGroup3 = new PointGroupHistClass(
//       startPosX,
//       startPosY,
//       startPosZ,
//       endPosX,
//       endPosY,
//       endPosZ,
//       numPoint,
//       numPlot,
//       steel->tempLiquid,
//       myCube
//       );
//  TempClass* caseLinearFVM1 = new TempClass(
//       steel,
//       myCube,
//       myLaserLinear,
//       timeLinearFVM,
//       myPointGroup3
//       );
//  caseLinearFVM1->setNumPlot(numPlot);
//  caseLinearFVM1->calTempHistFVM("linearFVM_1_");
//
//  // linear heat source analytical at the center
//  // time plan for linear analytical
//  cout << *timeLinearAna << endl;
//
//  TempClass* caseLinearAna1 = new TempClass(
//       steel,
//       myCube,
//       myLaserLinear,
//       timeLinearAna,
//       myPointGroup3
//       );
//  caseLinearAna1->setNumPlot(numPlot);
//  caseLinearAna1->calTempHistAnalytical("linearAna_1_");




//
////  string name1 = "statCenter";
////  myTempCase1->calTempHistFVM(name1);
//  myTempCase1->calTempHistAnalytical("testAnaStationary_");
//


//  TempClass* myTempCase2 = new TempClass(steel, myCube, myLaser, myTimePlan);
//
//  string name2 = "stationaryAna_1_";
//  myTempCase2->calTempHistFVM(name2);


//  // more close to boundary
//  myCube->setStartPositionXYZ(0.5e-3, 0.5e-3, 0);
//  cout << *myCube << endl;
//
//  TempClass* myTempCase3 = new TempClass(steel, myCube, myLaser, myTimePlan);
//  string name3 = "statNotCenter2_";
//
//  myTempCase3->calTempHistFVM(name3);
//
//  // more close to boudary
//  myCube->setStartPositionXYZ(0.75e-3, 0.75e-3, 0);
//  cout << *myCube << endl;
//
//  TempClass* myTempCase4 = new TempClass(steel, myCube, myLaser, myTimePlan);
//  string name4 = "statNotCenter3_";
//
//  myTempCase4->calTempHistFVM(name4);

//  // test linear scan path
//  LaserClass* myLaserLinear = new LaserClass(scanPatternIndexLinear);
//  cout << *myLaserLinear << endl;
//
//  // set the scan time
//  myTimePlan->setTotalTime(2e-3);
//  cout << *myTimePlan << endl;
//
//  TempClass* myTempLinear1 = new TempClass(steel, myCube, myLaserLinear, myTimePlan);
//  myTempLinear1->calTempHistFVM("FVM_linear1_");
//
//  // make the scan more close to the boundary
//  myCube->setStartPositionXYZ(0., 0.1e-3, 0.);
//  TempClass* myTempLinear2 = new TempClass(steel, myCube, myLaserLinear, myTimePlan);
//  myTempLinear2->calTempHistFVM("FVM_linear2_");
//
//  // make the scan more close to the boundary
//  myCube->setStartPositionXYZ(0., 0.2e-3, 0.);
//  TempClass* myTempLinear3 = new TempClass(steel, myCube, myLaserLinear, myTimePlan);
//  myTempLinear3->calTempHistFVM("FVM_linear3_");


  // build and test multi line scan track



  return 0;
}
