/* 
 * Author: zyuhao 
 * Date: 26/04/22
 * Content:  Source file of the PointHistClass
 * contains the variables and functions of creating, writing
 * the thermal history of one points in the cube
 * 
 */

//

#include "PointHistClass.h"

// initialize
void PointHistClass::initializePointHist(
     double inPointPosX,
     double inPointPosY,
     double inPointPosZ,
     const CubeClass *inMyCube,
     int inNumPlot,
     double inTempLimit
     )
{

  numPlot = inNumPlot;

  // set the upper limit of the temperature to write
  tempLimit = inTempLimit;

  // set xyzSize
  xSize = inMyCube->xsize;
  ySize = inMyCube->ysize;
  zSize = inMyCube->zsize;

  // allocate arrays
  tempArray = new double[numPlot+1];
  timeArray = new double[numPlot+1];
  // initialize to zero
  for (int i = 0; i < numPlot + 1; i++)
  {
    tempArray[i] = 0.;
    timeArray[i] = 0.;
  }

  // calculate the index position
  // check range of the point position
  if (inPointPosX < 0. || inPointPosX > inMyCube->xLen)
  {
    cout << "Point position in X direction is out of range!" << endl;
    cout << "Your value is: " << inPointPosX << endl;
    exit(3);
  }
  else if (inPointPosY < 0. || inPointPosY > inMyCube->yLen)
  {
    cout << "Point position in Y direction is out of range!" << endl;
    cout << "Your value is: " << inPointPosY << endl;
    exit(3);
  }
  else if (inPointPosX < 0. || inPointPosZ > inMyCube->zLen)
  {
    cout << "Point position in Z direction is out of range!" << endl;
    cout << "Your value is: " << inPointPosZ << endl;
    exit(3);
  }
  else
  {
    // calculate the index of the point
    pointIndexX = int(inPointPosX / inMyCube->meshSize);
    pointIndexY = int(inPointPosY / inMyCube->meshSize);
    pointIndexZ = int(inPointPosZ / inMyCube->meshSize);

    pointPosX = inPointPosX;
    pointPosY = inPointPosY;
    pointPosZ = inPointPosZ;
  }
}

// store temp and time of one point into one point
void PointHistClass::storePointOneTime(
     int timeStepCurrent,
     double timeCurrent,
     double* temp)
{
  if (timeStepCurrent >= 0 && timeStepCurrent < numPlot+1)
  {
    timeArray[timeStepCurrent] = timeCurrent;
    // record the temp blow the upper limit
    tempArray[timeStepCurrent] =
         min(
              temp[pointIndexX *  ySize * zSize + pointIndexY * zSize
              + pointIndexZ],
              tempLimit
              );
//    cout << pointIndexX *  ySize * zSize + pointIndexY * zSize
//            + pointIndexZ << endl;
  }
  else
  {
    cout << "timeStepCurrent is out of range!" << endl;
    cout << "timeStepCurrent: " << timeStepCurrent << endl;
    cout << "max valid index: " << numPlot << endl;
    //exit(3);
  }
}

// write thermal hist of one point into dat file
// the thermal hist of one point will be several line of values with a space
// between each value
// 1st line: 3 elements are the position of point along x,y,z coordinates
// 2ed line: number of time steps
// 3rd line: each value of time array
// 4th line: each value of temperature
void PointHistClass::writePointHist(string fileName)
{
  ofstream outFile;
  outFile.open(fileName);

  if (outFile.fail())
  {
    cout << "Fail to write data into " << fileName << endl;
    exit(3);
  }
  else
  {
    // commented line are used for plotting one point thermal history in tecplot
    // .dat for tecplot
//    outFile << "VARIABLES=\"Time\" \"Temp\" "<< endl;
//    outFile << "zone t=\" point \", i=" << numPlot+1
//           << ", j=" << numPlot+1 <<", f=point" << endl;
//    for (int i = 0; i < numPlot + 1; i++)
//    {
//      outFile << timeArray[i] << " " << tempArray[i] << endl;
//    }


    //.dat for python

    // first line xyz position
    outFile << pointPosX << " " << pointPosY << " " << pointPosZ << endl;
    // second line num of time steps
    outFile << numPlot+1 << endl;
    // third line element of time array
    for (int i = 0; i < numPlot + 1; i++)
    {
      outFile << timeArray[i] << " ";
    }
    outFile << endl;
    // forth line element of temperature array
    for (int i = 0; i < numPlot + 1; i++)
    {
      outFile << tempArray[i] << " ";
    }
    outFile << endl;
  }

  outFile.close();
}