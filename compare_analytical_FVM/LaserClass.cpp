/* 
 * Author: zyuhao 
 * Date: 2022-04-08
 * Content: source code of the LaserClass
 * contain the physic parameters of the laser
 * and the function describing the laser position wrt current time
 * 
 */

//

#include "LaserClass.h"

// default ctor
LaserClass::LaserClass(int inSPIndex)
{
  laserPower = 280.;
  beamRadiusX = 50e-6;
  beamRadiusY = 50e-6;
  beamRadiusZ = 50e-7;

  scanPatternIndex = inSPIndex;
  scanTime = 0.;

  if(scanPatternIndex == INDEX_SCANNING_PATH_STATIONARY)
  {
    scanTime = TIME_LASER_ACTIVE_STATIONARY;
  }
  else if (scanPatternIndex == INDEX_SCANNING_PATH_LINEAR)
  {
    scanTime = TIME_LASER_ACTIVE_LINEAR;
  }
}

// function calculate the position of the laser
// based on the given time timeVal
// return true if laser is still active and update laserPositionX/Y by reference
// return false if the laser is inactive
bool LaserClass::calLaserPosition(
     double &laserPositionX,
     double &laserPositionY,
     double timeVal
     )
{
  // if the scan path is stationary point
  if (scanPatternIndex == INDEX_SCANNING_PATH_STATIONARY)
  {
    if (timeVal <= scanTime)
    {
      laserPositionX = 1e-3;
      laserPositionY = 1e-3;
      return true;
    }
    else
    {
      return false;
    }
  }
  // if the scan path is a straight line
  else if (scanPatternIndex == INDEX_SCANNING_PATH_LINEAR)
  {
    if (timeVal <= scanTime)
    {
      laserPositionX = 0.5e-3 + 1 * timeVal;
      laserPositionY = 1e-3;
      return true;
    }
    else
    {
      return false;
    }
  }
//  // if the scan path is a spiral shape (to be continued if necessary)
//  else if (scanPatternIndex == scanPatternIndexSpiral)
//  {
//    if (timeVal <= 1e-3)
//  }
//  else
//  {
//    return false;
//  }
  else if (scanPatternIndex == INDEX_SCANNING_PATH_MULTI_LINE)
  {
    return calLaserPositionMultiLine(
         laserPositionX,
         laserPositionY,
         timeVal
         );
  }
}

// function calculate the instentaneous speed of the laser beam
// based on the given time timeVal
// return true if the laser is still on and update laserSpeed by reference
// return false if the laser if off
bool LaserClass::calLaserSpeed(
     double &laserSpeed,
     double timeVal
     )
{
  // if the scan pattern is stationary point
  if (scanPatternIndex == INDEX_SCANNING_PATH_STATIONARY)
  {
    if (timeVal <= scanTime)
    {
      laserSpeed = 0.;
      return true;
    }
    else
    {
      return false;
    }
  }
  // if the scan path is a straight line
  else if (scanPatternIndex == INDEX_SCANNING_PATH_LINEAR)
  {
    if (timeVal <= scanTime)
    {
      laserSpeed = 1.;
      return true;
    }
    else
    {
      return false;
    }
  }
  // if the scan path is a spiral shape (to be continued if necessary)
//  else if (scanPatternIndex == INDEX_SCANNING_PATH_SPIRAL)
//  {
//    if (timeVal <= scanTime)
//    {
//      laserSpeed = 0.;
//      return true;
//    }
//    else
//    {
//      return false;
//    }
//  }
  else if (scanPatternIndex == INDEX_SCANNING_PATH_MULTI_LINE)
  {
    return calLaserSpeedMultiLine(
         laserSpeed,
         timeVal
         );
  }
}

// virtual function updating xy positions of laser for multi line scanning
bool LaserClass::calLaserPositionMultiLine(
     double &laserPositionX,
     double &laserPositionY,
     double timeVal
     ){}

// virtual function updating abs speed of laser for multi line scanning
bool LaserClass::calLaserSpeedMultiLine(
     double &laserSpeed,
     double timeVal
     ){}


// overloading << for LaserClass
ostream& operator<<(ostream& os, const LaserClass& myLaser)
{
  os << endl << "Laser parameters list: " << endl
       << "  Laser power:      " << myLaser.laserPower << endl
       << "  Beam radius of X: " << myLaser.beamRadiusX << endl
       << "  Beam radius of Y: " << myLaser.beamRadiusY << endl
       << "  Beam radius of Z: " << myLaser.beamRadiusZ << endl
       << "  Scan time:        " << myLaser.scanTime << endl
       << "  Scan pattern:     " << myLaser.scanPatternIndex << endl;

  return os;
}