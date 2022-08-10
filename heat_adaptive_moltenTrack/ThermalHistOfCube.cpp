//
// Created by zyuhao on 6/7/22.
// implementation file ThermalHistOfCube class
//

#include "ThermalHistOfCube.h"


ThermalHistOfCube::ThermalHistOfCube()
{
  xLen = 0;
  yLen = 0;
  zLen = 0;
  meshSize = 0;
  xsize = 0;
  ysize = 0;
  zsize = 0;
  numPlot = 0;
  computingCaseObject = nullptr;
}


ThermalHistOfCube::ThermalHistOfCube(
     double inXLen,
     double inYLen,
     double inZLen,
     double inMeshSize, // size of mesh
     double inTimeLen, // time length
     int inNumPlot,
     ComputingCase* inComputingCase
     )
{
  xLen = inXLen;
  yLen = inYLen;
  zLen = inZLen;
  meshSize = inMeshSize;
  timeLen = inTimeLen;
  numPlot = inNumPlot;
  computingCaseObject = inComputingCase;

  // calculate xyz size
  xsize = int(round(xLen / meshSize));
  ysize = int(round(yLen / meshSize));
  zsize = int(round(zLen / meshSize));

  // initialize the static member in each 3d index MoltenPoolPoint class
  // used for operator == and < overwriting (sort and unique)
  MoltenPoolPoint::xsize = xsize;
  MoltenPoolPoint::ysize = ysize;
  MoltenPoolPoint::zsize = zsize;


  double timeStep = inTimeLen / (inNumPlot - 1); // step size of time

  // allocate vectors
  tempField = vector<vector<vector<double>>>
       (xsize, vector<vector<double>>(ysize, vector<double>(zsize, 0.)));

  for (int i = 0; i < xsize; i++)
  {
    xc.push_back(( i + 0.5) * meshSize);
  }

  for (int i = 0; i < ysize; i++)
  {
    yc.push_back(( i + 0.5) * meshSize);
  }

  for (int i = 0; i < zsize; i++)
  {
    zc.push_back(( i + 0.5) * meshSize);
  }

  for (int i = 0; i < numPlot; i++)
  {
    timeList.push_back(i * timeStep);
  }
}


void ThermalHistOfCube::calThermalHist(string fileName)
{
  int l; // loop variables

  // loop over timeList
  for (l  = 0; l < numPlot; l++)
  {
    calThermalFieldAtOneTime(timeList[l]);

    cout << "TimeStep " << l+1 << " / " << numPlot
         << " calculating complete!" << endl;

    // store the plot of current plot
    string fileNameCurrent = "./tecplot2d/" +
         fileName + "_" + to_string(l) + ".dat";
    storeThermalField(fileNameCurrent, timeList[l], tempField);
  }
}


void::ThermalHistOfCube::calThermalFieldAtOneTime(double timeCurrent)
{
  int i, j, k; // loop variable

  // update the current time for computing case
  computingCaseObject->setTimeCurrent(timeCurrent);

  // pre-calcualte all the node info used for integration shared by all POI
  computingCaseObject->calAbssicaAndPositionList();

  // calculate temperature field
  for (i = 0; i < xsize; i++)
  {
    for (j = 0; j < ysize; j++)
    {
      for (k = 0; k < zsize; k++)
      {
        tempField[i][j][k] =
             computingCaseObject->calTempOnePoint(
                  xc[i],
                  yc[j],
                  zc[k]
             );
      }
    }
  }
}





void ThermalHistOfCube::storeThermalField(
     string fileName,
     double solutionTime,
     vector<vector<vector<double>>> inTempField
     )
{
  cout << "Starting writing data into " << fileName << ". . ." << endl;

  ofstream outFile;
  outFile.open(fileName.c_str());

  if (outFile.fail())
  {
    cout << "Failed to open the file!!" << endl;
  }
  else
  {
    outFile << "VARIABLES=\"X\" \"Y\" \"Temp\" "<< endl;
    outFile << "zone t=\"" << solutionTime << "\", i=" << xsize
            << ", j=" << ysize <<", f=point" << endl;
    outFile << "solutiontime=" << solutionTime << endl;

    for (int i = 0; i < xsize; i++)
    {
      for (int j = 0; j < ysize; j++)
      {
        outFile << xc[i] << " " << yc[j] << " "
                << inTempField[i][j][0] << endl;
      }
    }
  }
  outFile.close();
  cout << "Write finished!" << endl;
}
