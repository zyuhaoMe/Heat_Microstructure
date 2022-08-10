//
// Created by zyuhao on 6/4/22.
//

#include "MultiTrackLaserClass.h"

bool MultiTrackLaserClass::calLaserPositionMultiLine(
     double &laserPositionX,
     double &laserPositionY,
     double timeVal
     )
{
  int i = 0; // index

  bool laserIsOn;

  if (timeVal >= timeEndList.back())
  {
    // if the all the scanning has already been finished
    laserIsOn = false;
    return laserIsOn;
  }
  else
  {
    // decide which interval the timeVal is in
    while (true)
    {
      if (timeVal <= timeEndList[i])
      {
        // laser is in interval i and laser is on
        laserIsOn = true;
        laserPositionX = posStartXList[i] + (posEndXList[i] - posStartXList[i])
             * (timeVal - timeStartList[i]) / (timeEndList[i] - timeStartList[i]);
        laserPositionY = posStartYList[i] + (posEndYList[i] - posStartYList[i])
             * (timeVal - timeStartList[i]) / (timeEndList[i] - timeStartList[i]);
        return laserIsOn;
      }
      // if timeVal > timeEndList[i] and timeVal < timeEndList.back(),
      // this means i + 1 is still a valid index
      else if (timeVal >= timeStartList[i + 1])
      {
        // laser is not in the current time interval, increase index
        i++;
      }
      else
      {
        // laser spot is moving to next line segment but currently, laser is off
        laserIsOn = false;
        return laserIsOn;
      }
    }
  }

}

bool MultiTrackLaserClass::calLaserSpeedMultiLine(
     double &speed,
     double timeVal
     )
{
  int i = 0; // index

  bool laserIsOn;

  if (timeVal > timeEndList.back())
  {
    // if the all the scanning has already been finished
    laserIsOn = false;
    return laserIsOn;
  }
  else
  {
    // decide which interval the timeVal is in
    while (true)
    {
      if (timeVal <= timeEndList[i])
      {
        // laser is in interval i and laser is on
        laserIsOn = true;
        speed = speedList[i];
        return laserIsOn;
      }
        // if timeVal > timeEndList[i] and timeVal < timeEndList.back(),
        // this means i + 1 is still a valid index
      else if (timeVal >= timeStartList[i + 1])
      {
        // laser is not in the current time interval, increase index
        i++;
      }
      else
      {
        // laser spot is moving to next line segment but currently, laser is off
        laserIsOn = false;
        return laserIsOn;
      }
    }
  }
}


void MultiTrackLaserClass::add1LaserOffIntervalAnd1LineTrack(
     double inLaserOffInterval,
     double inPosStartX,
     double inPosStartY,
     double inPosEndX,
     double inPosEndY,
     double inSpeed
     )
{
  // add position node
  posStartXList.push_back(inPosStartX);
  posStartYList.push_back(inPosStartY);
  posEndXList.push_back(inPosEndX);
  posEndYList.push_back(inPosEndY);

  // add speed
  speedList.push_back(inSpeed);

  // add number of line segment
  numLineTrack++;

  // time of line scan
  double timeScan = sqrt( (inPosEndX - inPosStartX) * (inPosEndX - inPosStartX)
       + (inPosEndY - inPosStartY) * (inPosEndY - inPosStartY) ) / inSpeed;

  // push the new off time into the vector
  // for the first line segment,
  if (timeStartList.empty())
  {
    timeStartList.push_back(inLaserOffInterval);
    timeEndList.push_back(inLaserOffInterval + timeScan);
  }
  // for non-first line segment, add time based on the last element of timeList
  else
  {
    timeStartList.push_back(timeEndList.back() + inLaserOffInterval);
    timeEndList.push_back(timeStartList.back() + timeScan);
  }
}

void MultiTrackLaserClass::writeScanTrackInto(
     string fileName
     )
{
  // because the overall speed do not vary too much for all line segment,
  // it is suitable to define the resolution to discrete points
  // according to time interval length of each segment

  int numPointMinLine = 200;
  double timeInterval;
  double timeIntervalMin = timeEndList[0] - timeStartList[0];

  // traverse all the line segment, find the shortest time interval
  for (int i = 0; i < numLineTrack; i++)
  {
    timeInterval = timeEndList[i] - timeStartList[i];
    timeIntervalMin = min(timeIntervalMin, timeInterval);
  }

  double timeStep = timeIntervalMin / numPointMinLine;

  // write all the time, position along x and y into a file
  string fileNameWithPath = "./plot_scan_track/" + fileName + ".dat";
  ofstream outFile;
  outFile.open(fileNameWithPath.c_str());

  if (outFile.fail())
  {
    cout << "failed to open file" << endl;
    exit(3);
  }

  double timeCurrent = 0.;
  double laserPositionX, laserPositionY;

  while (timeCurrent < timeEndList.back())
  {
    if (calLaserPositionMultiLine(
         laserPositionX,
         laserPositionY,
         timeCurrent
         ))
    {
      outFile << timeCurrent << " "
           << laserPositionX << " "
           << laserPositionY << endl;
    }
    timeCurrent += timeStep;
  }

  outFile.close();
}