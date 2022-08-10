/* 
 * Author: zyuhao 
 * Date: 2022-04-09
 * Content:  source code file of TimeClass,
 * contain total time, time step and functions of determining time step in
 * different ways
 */

//

#include "TimeClass.h"


// default ctor
TimeClass::TimeClass(CubeClass *myCube, MaterialClass *myMaterial)
{
  cubePtr = myCube;
  materialPtr = myMaterial;
  totalTime = 1e-3;
  timeStep = 1.5/200;
  calNumTimeStep();

}


// set the time step using a self-defined timeStep
bool TimeClass::setTimeStep(double myTimeStep)
{
  if (myTimeStep >= totalTime || myTimeStep <= 0.)
  {
    return false;
  }
  else
  {
    timeStep = myTimeStep;
    calNumTimeStep();
    return true;
  }
}

// set the time step based on the Von Neumann stability
// the result depends on the mesh size and thermal diffusivity
void TimeClass::setTimeStepWithVonNeu(double safeCoeff)
{
  timeStep = safeCoeff * (cubePtr->meshSize) * (cubePtr->meshSize)
       / (6 * materialPtr->thermalDiff);
  calNumTimeStep();
}

// set the time step based on the number of plot
void TimeClass::setTimeStepWithExactPlot(int inNumPlot)
{
  timeStep = totalTime / (inNumPlot);
  calNumTimeStep();
}


// function that set the totalTime
void TimeClass::setTotalTime(double myTimeTotal)
{
  totalTime = myTimeTotal;
  calNumTimeStep();
}

// function that recalculate the numTimeStep after resetting the totalTime
// or timeStep
void TimeClass::calNumTimeStep()
{
  numTimeStep = int(ceil(totalTime / timeStep));
}

// overloading << for TimeClass
ostream& operator<<(ostream& os, const TimeClass& myTimePlan)
{
  os << endl << "Time parameters list:" << endl
      << "  Total time:    " << myTimePlan.totalTime << endl
      << "  Time step:     " << myTimePlan.timeStep << endl
      << "  Num time step: " << myTimePlan.numTimeStep << endl;

  return os;
}