//
// Created by zyuhao on 6/8/22.
//

#include "ThermalHistOfCubeMoltenPoolTrack.h"

ThermalHistOfCubeMoltenPoolTrack::ThermalHistOfCubeMoltenPoolTrack(
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
     ):
     ThermalHistOfCube(
          inXLen,
          inYLen,
          inZLen,
          inMeshSize, // size of mesh
          inTimeLen, // time length
          inNumPlot,
          inComputingCase
          )
{
  startingTime = inStartingTime;
  timeStepSize = inTimeStepSize;
  tempBound = inTempBond;

  // calculate total number of timeStep,
  // will exclude the final time step (timeStarting + timeLen)
  numTimeStep = int(round(timeLen / timeStepSize)) - 1;

  // allocate tempHasCompted vector
  tempHasComputed = vector<vector<vector<bool>>>
       (xsize, vector<vector<bool>>(ysize, vector<bool>(zsize, false)));

  // allocate tempToPlot
  tempToPlot = vector<vector<vector<double>>>
       (xsize, vector<vector<double>>(ysize, vector<double>(zsize, TEMP_AMBIENT_DEFAULT)));
}


void ThermalHistOfCubeMoltenPoolTrack::detectAndAddMoltenPoolFromTempField()
{
  int i, j, k; // loop variable
  for (i = 0; i < xsize; i++)
  {
    for (j = 0; j < ysize; j++)
    {
      for (k = 0; k < zsize; k++)
      {
        if(tempField[i][j][k] >= tempBound)
        {
          // add into the molten pool list
          moltenPoolList_new.emplace_back(MoltenPoolPoint(
               i,
               j,
               k
               ));
        }
      }
    }
  }
}


void ThermalHistOfCubeMoltenPoolTrack::trackMoltenPoolAtNewTime(
     double inTimeCurrent
     )
{
  // variable used for debugging, commented after finishing debugging
  int numWhileLoop = 0.;


  // tempHasComputed 3d vector of last time step need to be initialized to false
  for (int i = 0; i < xsize; i++)
  {
    for (int j = 0; j < ysize; j++)
    {
      for (int k = 0; k < zsize; k++)
      {
       tempHasComputed[i][j][k] = false;
      }
    }
  }

  // set the time of computingCase as timeCurrent
  computingCaseObject->setTimeCurrent(inTimeCurrent);

  // update the integration nodes and inner constants
  computingCaseObject->calAbssicaAndPositionList();

  // check the laser beam position and add the point close to the laser beam
  // into old molten pool list (after this step, moltenPoolList_old become
  // potential molten pool list)
  computingCaseObject->addPotentialPointNearLaser(
       inTimeCurrent,
       meshSize,
       moltenPoolList_old,
       0,
       xsize - 1,
       0,
       ysize -1
       );

  // check whether at least one point has been added into the molten pool
  cout << "Num of potential molten points: " << moltenPoolList_old.size() << endl;


  // calculate temperature of points in old molten pool list
  // 1. label it to be computed
  // 2. if above tempBond, push into new vector
  for (
       auto it = moltenPoolList_old.begin();
       it != moltenPoolList_old.end();
       ++it
       )
  {
    double tempLocal; // local point
    tempLocal =  computingCaseObject->calTempOnePoint(
         xc[it->indexX],
         yc[it->indexY],
         zc[it->indexZ]
         );

    // label this point to be computed
    tempHasComputed[it->indexX][it->indexY][it->indexZ] = true;

    // whether to add to new molten pool list
    if (tempLocal > tempBound)
    {
      moltenPoolList_new.emplace_back(MoltenPoolPoint(
           it->indexX,
           it->indexY,
           it->indexZ
           ));

      // set the corresponding value in tempField
      tempField[it->indexX][it->indexY][it->indexZ] = tempLocal;
    }
  }


  do
  {
    // iterate the listToCompute, calculate the temperature and add into
    // the moltenPoolList_New
    for (
         auto it = listToCompute.begin();
         it != listToCompute.end();
         ++it
         )
    {
      double tempLocal;

      // calculate the temp of POI
      tempLocal = computingCaseObject->calTempOnePoint(
           xc[it->indexX],
           yc[it->indexY],
           zc[it->indexZ]
      );

      // label this point has been calculated
      tempHasComputed[it->indexX][it->indexY][it->indexZ] = true;


      // if in the listToCompute, this point is above the tempLimit
      // put the index into the molten pool set
      // put the index into the molten pool set
      if (tempLocal > tempBound)
      {
        moltenPoolList_new.emplace_back(MoltenPoolPoint(
             it->indexX,
             it->indexY,
             it->indexZ
        ));

        // set the corresponding value in tempField
        tempField[it->indexX][it->indexY][it->indexZ] = tempLocal;
      }
    }

    // used for debugging, print the number of while loop and number of elements
    // in the molten pool
    cout << "Num of while loop: " << numWhileLoop << endl;
    numWhileLoop++;
    cout << "Len of listToCompute: " << listToCompute.size() << endl;

    // clear listToCompute
    listToCompute.clear();


    // detect the neighbour of new molten pool and add to the list to compute
    for (
         auto it = moltenPoolList_new.begin();
         it != moltenPoolList_new.end();
         ++it
         )
    {
      // check neighbour of i-1
      if (
           it->indexX != 0 &&
           !tempHasComputed[it->indexX - 1][it->indexY][it->indexZ]
           )
      {
        listToCompute.emplace_back(MoltenPoolPoint(
             it->indexX - 1,
             it->indexY,
             it->indexZ
        ));
      }

      // check neighbour of i+1
      if (
           it->indexX != xsize &&
           !tempHasComputed[it->indexX + 1][it->indexY][it->indexZ]
           )
      {
        listToCompute.emplace_back(MoltenPoolPoint(
             it->indexX + 1,
             it->indexY,
             it->indexZ
        ));
      }

      // check neighbour of j-1
      if (
           it->indexY != 0 &&
           !tempHasComputed[it->indexX][it->indexY - 1][it->indexZ]
           )
      {
        listToCompute.emplace_back(MoltenPoolPoint(
             it->indexX,
             it->indexY - 1,
             it->indexZ
        ));
      }

      // check neighbour of j+1
      if (
           it->indexY != ysize &&
           !tempHasComputed[it->indexX][it->indexY + 1][it->indexZ]
           )
      {
        listToCompute.emplace_back(MoltenPoolPoint(
             it->indexX,
             it->indexY + 1,
             it->indexZ
        ));
      }

      // check neighbour of k-1
      if (
           it->indexZ != 0 &&
           !tempHasComputed[it->indexX][it->indexY][it->indexZ - 1]
           )
      {
        listToCompute.emplace_back(MoltenPoolPoint(
             it->indexX,
             it->indexY,
             it->indexZ - 1
        ));
      }

      // check neighbour of k+1
      if (
           it->indexZ != zsize &&
           !tempHasComputed[it->indexX][it->indexY][it->indexZ + 1]
           )
      {
        listToCompute.emplace_back(MoltenPoolPoint(
             it->indexX,
             it->indexY,
             it->indexZ + 1
        ));
      }
    }
  } while (!listToCompute.empty()); //stop when there is no neighbour to calculate

  // Now in tempField, only the index in the moltenPoolList is valid
  // other points is not valid for current time step


}


void ThermalHistOfCubeMoltenPoolTrack::clearAndSwapMoltenPoolList()
{
  moltenPoolList_old.clear(); // clear old molten pool list
  moltenPoolList_new.swap(moltenPoolList_old); // exchange old list with new list
}


void ThermalHistOfCubeMoltenPoolTrack::trackMoltenPoolUntilEnd(string fileName)
{
  // assume that numTimeStep = n * numTimeStep_plot, where n is an integer
  int n = int(round((numTimeStep + 1) / numPlot));

  double timeCurrent = startingTime;

  // current plot number
  int plotNumCurrent = 0;

  // at the starting time, calculate the whole 3d tempField
  calThermalFieldAtOneTime(timeCurrent);

  // process the 3d tempField at the starting time, extract molten pool
  detectAndAddMoltenPoolFromTempField();

  storeThermalField(
       (fileName + "_" + to_string(plotNumCurrent) + ".dat"),
       timeCurrent,
       tempField
       );

  cout << "TimeStep " << 1 << " / " << numTimeStep + 1
       << " calculating complete!" << endl;

  for (int i = 1; i < numTimeStep; i++)
  {
    // clear and swap moltenPoolList
    clearAndSwapMoltenPoolList();

    timeCurrent = timeStepSize * i; // update current time

    // track molten pool at this time,
    // will call computingCaseObject to compute points temperatures
    trackMoltenPoolAtNewTime(timeCurrent);

    // print in terminal, current time has finished
    cout << "TimeStep " << i + 1 << " / " << numTimeStep + 1
         << " calculating complete!" << timeCurrent << endl;


    // transform and plot tempField at specific temp step to write data to plot
    if (i % n == 0)
    {
      plotNumCurrent++;
      transformMoltenPoolDataIntoTempField();
      storeThermalField(
           ("./tecplot2d/" + fileName + "_" + to_string(plotNumCurrent) + ".dat"),
           timeCurrent,
           tempToPlot
           );

    }
  }
}


void ThermalHistOfCubeMoltenPoolTrack::transformMoltenPoolDataIntoTempField()
{
  // initialize the tempToPlot to ambient temperature
  for (int i = 0; i < xsize; i++)
  {
    for (int j = 0; j < ysize; j++)
    {
      for (int k = 0; k < zsize; k++)
      {
        tempToPlot[i][j][k] = TEMP_AMBIENT_DEFAULT;
      }
    }
  }

  // according to the moltenPoolList_new, set the point
  for (
       auto it = moltenPoolList_new.begin();
       it != moltenPoolList_new.end();
       ++it
       )
  {
    tempToPlot[it->indexX][it->indexY][it->indexZ] =
         tempField[it->indexX][it->indexY][it->indexZ];
  }
}







































