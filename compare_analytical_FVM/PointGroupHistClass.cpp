/* 
 * Author: zyuhao 
 * Date: 27/04/22
 * Content:  source file of PointGroupClass
 * create a group of points to store thermal history
 * iteratively call the function of each PointHistClass object
 * to store and write
 * 
 */

//

#include "PointGroupHistClass.h"


// ctor
PointGroupHistClass::PointGroupHistClass(
     double inStartPosX,
     double inStartPosY,
     double inStartPosZ,
     double inEndPosX,
     double inEndPosY,
     double inEndPosZ,
     int inNumPoint,
     int inNumPlot,
     double inTempLimit,
     const CubeClass *inMyCube
     )
{
  numPoint = inNumPoint;

  // dynamically allocate object arrays and other arrays
  pointArray = new PointHistClass[numPoint];
  pointPosArrayX = new double[numPoint];
  pointPosArrayY = new double[numPoint];
  pointPosArrayZ = new double[numPoint];
  fileNameArray = new string[numPoint];


  // distance between neighbouring points along 3d
  double nbDistX, nbDistY, nbDistZ;
  nbDistX = (inEndPosX - inStartPosX) / (numPoint - 1);
  nbDistY = (inEndPosY - inStartPosY) / (numPoint - 1);
  nbDistZ = (inEndPosZ - inStartPosZ) / (numPoint - 1);

  // initialize arrays
  for (int i = 0; i < numPoint; i++)
  {
    // initialize the file name string
    fileNameArray[i] = "To be given";

    // position along 3d inclusively;
    pointPosArrayX[i] = inStartPosX + i * nbDistX;
    pointPosArrayY[i] = inStartPosY + i * nbDistY;
    pointPosArrayZ[i] = inStartPosZ + i * nbDistZ;

    // initialize each PointHistClass
    pointArray[i].initializePointHist(
         pointPosArrayX[i],
         pointPosArrayY[i],
         pointPosArrayZ[i],
         inMyCube,
         inNumPlot,
         inTempLimit
         );
  }
}

// set the file name
void PointGroupHistClass::setFileName(const string tempObjectFileName)
{
  delete [] fileNameArray;
  // re allocate
  fileNameArray = new string[numPoint];
  for (int i = 0; i < numPoint; i++)
  {
    // file name of a point
    fileNameArray[i] = "./tecplot_points/" +
         tempObjectFileName + "Point_" + to_string(i) + ".dat";
  }
}

// store time and temp of all points at one time step
void PointGroupHistClass::storePointsOneTime(
     int timeStepCurrent,
     double timeCurrent,
     double *temp
     )
{
  for (int i = 0; i < numPoint; i++)
  {
    pointArray[i].storePointOneTime(timeStepCurrent, timeCurrent, temp);
  }
}

// write thermal history of all points into different dat file
void PointGroupHistClass::writePointsHist()
{
  for (int i = 0; i < numPoint; i++)
  {
    pointArray[i].writePointHist(fileNameArray[i]);
  }
}

// overloading << for PointGroupHistClass
ostream& operator<<(
     ostream& os,
     const PointGroupHistClass& myPointGroup
     )
{
  os << endl << "Point group parameters:" << endl;
  os << "  Total number of the points: " << myPointGroup.numPoint << endl;
  for (int i = 0; i < myPointGroup.numPoint; i++)
  {
    os << "  Point number: " << i << endl;
    os << "  Point position: X:" << myPointGroup.pointPosArrayX[i]
         << " Y: " <<  myPointGroup.pointPosArrayY[i]
         << " Z: " << myPointGroup.pointPosArrayZ[i] << endl;
    os << "  Output file name: " << myPointGroup.fileNameArray[i] << endl;
  }

  return os;
}