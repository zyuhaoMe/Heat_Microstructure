#include <iostream>
#include "constants.h"
#include "LaserClass.h"
#include "MaterialClass.h"
#include "MultiTrackLaserClass.h"
#include "ComputingCase.h"
#include "ThermalHistOfCube.h"
#include "ThermalHistOfCubeMoltenPoolTrack.h"

using namespace std;

int MoltenPoolPoint::xsize;
int MoltenPoolPoint::ysize;
int MoltenPoolPoint::zsize;

int main() {

  /*==========================================================================*/

  /* [1]
   * Build several material object
   */

  // build a material object
  auto* steel = new MaterialClass("steel");
  // print material parameter
  cout << *steel << endl;

  // build a material object named IN625
  auto* IN625 = new MaterialClass("IN625");
  // set the properties
  IN625->setDens(8440);
  IN625->setSpecHeat(500);
  IN625->setThermalCond(16);
  IN625->setTempLiquid(1623);
  IN625->setAbsorp(0.43);
  IN625->setConvecCoeff(10);
  IN625->setEmissivity(0.8);
  // print info
  cout << *IN625 << endl;

  // build a material object named Steel155
  auto* steel155 = new MaterialClass("Steel155");
  // set the properties
  steel155->setDens(7770);
  steel155->setSpecHeat(420);
  steel155->setThermalCond(17.8);
  steel155->setTempLiquid(1700);
  steel155->setAbsorp(0.22);
  steel155->setConvecCoeff(10);
  steel155->setEmissivity(0.8);
  // print info
  cout << *steel155 << endl;

  /*==========================================================================*/

  /* [2]
   * Build a multi line track MultiTrackLaserClass object
   * Write the path data into multiTrackLaser1 to pre-check the path
   */

  auto* multiTrackLaser1  = new MultiTrackLaserClass();
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

  multiTrackLaser1->writeScanTrackInto("Multi_track_data");

  /*==========================================================================*/

  /* [3]
   * Build a computing case which compose material object and laser object
   */

  auto* computingCaseSteel = new ComputingCase(
       multiTrackLaser1,
       steel
       );

  /* [3.1]
   * use a cube to test the thermal history
   */

  double startTime = omp_get_wtime();

  int numPlot = 100;
  auto* cubeTestCase = new ThermalHistOfCube(
       XLEN_CUBE_TEST,
       YLEN_CUBE_TEST,
       ZLEN_CUBE_TEST,
       MESHSIZE_CUBE_TEST,
       TOTAL_TIME_CUBE_TEST,
       numPlot,
       computingCaseSteel
       );
  cubeTestCase->calThermalHist("cubeTestCase");

  double endTime = omp_get_wtime();

  cout << "Running time: " << endTime - startTime << endl;

  /* [4]
   * use molten-tracking strategy to track the thermal history above
   * specific temperature with very fine time step
   */

  startTime = omp_get_wtime();

  numPlot = 100;
  double timeComputeStart = 0.;
  double timeStepSize = TOTAL_TIME_CUBE_TEST / numPlot;
  double tempBond = steel->tempLiquid;

  auto* trackMoltenPoolTestCube = new ThermalHistOfCubeMoltenPoolTrack(
       XLEN_CUBE_TEST,
       YLEN_CUBE_TEST,
       ZLEN_CUBE_TEST,
       MESHSIZE_CUBE_TEST,
       timeComputeStart,
       TOTAL_TIME_CUBE_TEST,
       timeStepSize,
       tempBond,
       numPlot,
       computingCaseSteel
       );

  trackMoltenPoolTestCube->trackMoltenPoolUntilEnd("moltenPoolCube");

  endTime = omp_get_wtime();
  cout << "Running time: " << endTime - startTime << endl;

  return 0;
}
