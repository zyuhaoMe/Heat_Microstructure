//
// Created by zyuhao on 6/6/22.
// implementation file of ComputingCase class

#include "ComputingCase.h"

ComputingCase::ComputingCase(
     MultiTrackLaserClass *inTrackObject,
     MaterialClass *inMaterialObject
     )
{
  trackLaserObject = inTrackObject;
  materialObject = inMaterialObject;
  timeCurrent = 0.;

  phyConstTime = materialObject->thermalDiff /
       (trackLaserObject->beamRadiusX * trackLaserObject->beamRadiusX);

  phyConstSpeed = trackLaserObject->beamRadiusX / materialObject->thermalDiff;

  phyConstMaterialLaser = 2. * materialObject->absorp * trackLaserObject->laserPower
       / (pow(PI/3., 1.5) * materialObject->thermalCond * trackLaserObject->beamRadiusX);

  timeCurrent_nd = timeCurrent * phyConstTime;
  tempAmbient = TEMP_AMBIENT_DEFAULT;
}


void ComputingCase::setTimeCurrent(double inTimeCurrent)
{
  timeCurrent = inTimeCurrent;

  // recalculate timeCurrent_nd
  timeCurrent_nd = timeCurrent * phyConstTime;
}


void ComputingCase::calAbssicaAndPositionList()
{

  // clear all lists
  abssicaList.clear();
  weightList.clear();
  laserPosXList.clear();
  laserPosYList.clear();
  innerConst1List.clear();
  innerConst2List.clear();
  innerConst3List.clear();

  // non-dimensional conduction time
  double timeDiff_nd;

  // pass the current time to the laser object
  // ask for current scan information
  double laserPositionX, laserPositionY, speed; // position of laser and scan speed
  int indexLineSegmentCurrent; // index of current line segment
  bool laserIsOnCurrently = true; // whether the laser is on or off currently

  cout << endl;

  trackLaserObject->calLaserInfoMultiLine(
       laserPositionX,
       laserPositionY,
       speed,
       indexLineSegmentCurrent,
       laserIsOnCurrently,
       timeCurrent
       );

  cout << endl;

  // calculate the thermal abssica of current line segment until the previous
  // full segment

  // if currently laser is off
  // all the previous line segment should have been scanned fully, nothing needs to do
  // if currently laser is on, current line segment has not been finished
  // add the abssica, weight, laser position and inner constant of current line segment

  if (laserIsOnCurrently) {
    timeDiff_nd = 0.;
    double timeDiff_nd_End = timeCurrent_nd -
                             trackLaserObject->timeStartList[indexLineSegmentCurrent] *
                             phyConstTime;
    calAbssicaAndPositionListOneTrack(
         timeDiff_nd,
         timeDiff_nd_End
    );
    // thermal contribution of current unfinished line segment has been
    // considered, minus line segment index by one
    indexLineSegmentCurrent--;
  }

  // for all the fully-scanned line segment, calculate their node information
  // Since the conduction time is in the inverse sequence of scan time,
  // do the for loop of each scanning track inversely
  for (int i = indexLineSegmentCurrent; i >= 0; i--) {
    calAbssicaAndPositionListOneTrack(
         timeCurrent_nd - trackLaserObject->timeEndList[i] * phyConstTime,
         timeCurrent_nd - trackLaserObject->timeStartList[i] * phyConstTime
         );
  }

//  // used for test: print out all the node information for current time step
//
//  cout << "TimeCurrent: " << timeCurrent << endl;
//  cout << "TimeCurrent_nd: " << timeCurrent_nd << endl;
//  unsigned vecLen = abssicaList.size();
//  for (unsigned i = 0; i < vecLen; i++)
//  {
//    cout << "Node Num: " << i << " / " << vecLen << ":" << endl;
//    cout << "Abssica: " << abssicaList[i] << endl;
//    cout << "Weight: " << weightList[i] << endl;
//    cout << "---------------------------------" << endl;
//  }
//  cout << "==================================" << endl;
}

void ComputingCase::calAbssicaAndPositionListOneTrack(
     double timeDiff_nd_start,
     double timeDiff_nd_end
     )
{
  double timeDiff_nd = timeDiff_nd_start; // current nd conduction time
  double speed; // current speed
  double speed_nd; // current nondimensional speed
  double prelimTimeSegSize_nd; // nd preliminary time segment size
  int numNodes; // number of nodes for each nd time segment
  double timeSegSize_nd;  // nd time segment size
  int i; // loop variable
  double abssicaCurrent, weightCurrent; // abssica and weight of current node
  double laserPosXCurrent, laserPosYCurrent; // laser pos of current node


  while (timeDiff_nd < timeDiff_nd_end)
  {
    // update the current speed
    trackLaserObject->calLaserSpeed(
         speed,
         (timeCurrent_nd - timeDiff_nd) / phyConstTime
         );

    // nd speed
    speed_nd = speed * phyConstSpeed;

    // determine preliminary segment size
    prelimTimeSegSize_nd = pow(2.0, floor(log2(1. + timeDiff_nd)));

    // number of nodes for current sub interval
    numNodes = int(max(2.0, 16. / prelimTimeSegSize_nd));

    // determine if within diffusion controlled regime
    // and set appropriate segment size
    if (speed_nd < 0.59)
    {
      timeSegSize_nd = prelimTimeSegSize_nd;
    }
    else
    {
      timeSegSize_nd = sqrt((12. * timeDiff_nd + 1.) * log(sqrt(2.)))
           / speed_nd;
    }

    // check timeDiff_nd_end
    if (timeDiff_nd + timeSegSize_nd > timeDiff_nd_end)
    {
      timeSegSize_nd = timeDiff_nd_end - timeDiff_nd;
    }

    // for current nd time segment, generate a gaussian integration case to get
    // node abssicas and weights
    gsl_integration_glfixed_table *table
         = gsl_integration_glfixed_table_alloc(numNodes);

    for (i = 0; i < numNodes; i++)
    {
      // decide the abssica and weight
      gsl_integration_glfixed_point(
           timeDiff_nd,
           timeDiff_nd + timeSegSize_nd,
           i,
           &abssicaCurrent,
           &weightCurrent,
           table
           );

      // laser pos of current node
      trackLaserObject->calLaserPosition(
           laserPosXCurrent,
           laserPosYCurrent,
           (timeCurrent_nd - abssicaCurrent) / phyConstTime
           );

      // push info of current integral node into vector
      abssicaList.push_back(abssicaCurrent);
      weightList.push_back(weightCurrent);
      laserPosXList.push_back(laserPosXCurrent);
      laserPosYList.push_back(laserPosYCurrent);

      // calculate and push inner constants of current integral node
      innerConst2List.push_back(12. * abssicaCurrent + 1.);
      innerConst3List.push_back(12. * abssicaCurrent);
      innerConst1List.push_back(1. /
      (innerConst2List.back() * sqrt(innerConst3List.back())));

    }
    // free table
    gsl_integration_glfixed_table_free(table);

    // increment conduction time
    timeDiff_nd += timeSegSize_nd;
  }
}


double ComputingCase::calTempOnePoint(
     double inXc,
     double inYc,
     double inZc
     )
{
  omp_set_num_threads(8);

  unsigned vecLen = abssicaList.size(); // total number of nodes 1-based

  double temp = 0.; // temperature at this point

  // too many calls of beamRadius
  // copy it locally
  double beamRadius = trackLaserObject->beamRadiusX;


#pragma omp parallel for reduction(+:temp)
  for (unsigned i = 0; i < vecLen; i++)
  {
    temp += weightList[i] * innerConst1List[i] *
            exp(-3. *(((inXc - laserPosXList[i]) / beamRadius
                  * (inXc - laserPosXList[i]) / beamRadius +
                  (inYc - laserPosYList[i]) / beamRadius
                  * (inYc - laserPosYList[i]) / beamRadius) / innerConst2List[i] +
                 (inZc * inZc) / beamRadius / beamRadius / innerConst3List[i]));
  }

  temp = temp * phyConstMaterialLaser + tempAmbient;

  return temp;
}


void ComputingCase::addPotentialPointNearLaser(
     double inTimeCurrent,
     double inMeshSize,
     vector<MoltenPoolPoint> &potentialMoltenPoolList,
     int xIndex_min,
     int xIndex_max,
     int yIndex_min,
     int yIndex_max
     )
{
  double laserPosX, laserPosY;
  int xIndexClose, yIndexClose;

  // if laser is on
  if(trackLaserObject->calLaserPositionMultiLine(
       laserPosX,
       laserPosY,
       inTimeCurrent
       ))
  {
    xIndexClose = int(round(laserPosX / inMeshSize));
    yIndexClose = int(round(laserPosY / inMeshSize));

    // check whether the point index is within the cube
    if (
         xIndexClose >= xIndex_min &&
         xIndexClose <= xIndex_max &&
         yIndexClose >= yIndex_min &&
         yIndexClose <= yIndex_max
         )
    {
      potentialMoltenPoolList.emplace_back(MoltenPoolPoint(
           xIndexClose,
           yIndexClose,
           0
           ));

      // check whether current close-to-beam point has already existed in
      // moltenPoolList_old using "unique"

      auto it = unique(potentialMoltenPoolList.begin(), potentialMoltenPoolList.end());

      potentialMoltenPoolList.resize(distance(potentialMoltenPoolList.begin(), it));
    }
  }
}

