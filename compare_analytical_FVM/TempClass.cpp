/* 
 * Author: zyuhao 
 * Date: 2022-04-11
 * Content:  source code file of TempClass, compute temperature field
 * contains old temperature 3d array and new temperature 3d array
 * defines relative functions of calculating temperature filed
 * using FVM or analytical method
 * defines functions that write temperature field into .dat file
 * 
 */

#include "TempClass.h"

// default ctor, load data from dependent object to and store it locally
// to prevent extra work
TempClass::TempClass(
     MaterialClass* myMaterial,
     CubeClass* myCube,
     LaserClass* myLaser,
     TimeClass* myTimePlan,
     PointGroupHistClass* myPointGroup
     )
{
  tempAmbient = 300.; // ambient temperature
  SBConstant = 5.67e-8; // stefan-bolzmann constant
  numPlot = 50; // number of plot

  // number of threads used for computing
  numThreads = 10;

  // make a copy of frequently used variables
  xSize = myCube->xsize;
  ySize = myCube->ysize;
  zSize = myCube->zsize;
  meshSize = myCube->meshSize;
  thermalDiff = myMaterial->thermalDiff;
  timeStep = myTimePlan->timeStep;
  numTimeStep = myTimePlan->numTimeStep;
  xc = myCube->xc;
  yc = myCube->yc;
  zc = myCube->zc;
  xLen = myCube->xLen;
  yLen = myCube->yLen;
  zLen = myCube->zLen;
  xStart = myCube->xStart;
  yStart = myCube->yStart;
  zStart = myCube->zStart;
  thisLaser = myLaser;
  thisMaterial = myMaterial;
  beamRadius = myLaser->beamRadiusX;
  scanTime = myLaser->scanTime;
  totalTime = myTimePlan->totalTime;

  thisPointGroup = myPointGroup;

  // calculate the physical constant
  ConstantFlux1 = 3 * myLaser->laserPower / (PI * beamRadius * beamRadius);
  ConstantFlux2 = thisMaterial->absorp * timeStep
       / (thisMaterial->dens * thisMaterial->specHeat * meshSize);
  ConstantAnaly = 2 * thisLaser->laserPower * myMaterial->absorp /
       (myMaterial->dens * myMaterial->specHeat * pow(PI / 3., 1.5));
  phyConst_nd = 2. * thisMaterial->absorp * thisLaser->laserPower
       / (pow(PI/3., 1.5) * thisMaterial->thermalCond * beamRadius);

  // calculate coefficient especially for myFunc
  coeff_a = myMaterial->emissivity * SBConstant;
  coeff_b = myMaterial->convecCoeff + 2.0 * myMaterial->thermalCond
       / (meshSize / 2.0);
  coeff_c0 = -1.0 * myMaterial->emissivity * SBConstant * pow(tempAmbient, 4.0)
             - 1.0 * myMaterial->convecCoeff * tempAmbient;
  coeff_c1 = -2.0 * myMaterial->thermalCond / (meshSize / 2.0);

  // dynamically allocate 3d array of tempOld
  tempOld = new double[xSize * ySize * zSize];
  tempOldX = new double[xSize * ySize * zSize];
  tempOldY = new double[xSize * ySize * zSize];

  // dynamically allocate 3d array of tempNew
  tempNew = new double[xSize * ySize * zSize];

  // initialize tempOld and tempNew
  initializeTemp(tempAmbient);
}

// function that initialize the tempOld and tempNew to the ambient temp
void TempClass::initializeTemp(double inTemp)
{
  int i; // loop var

  for (i = 0; i < xSize * ySize * zSize; i++)
  {
    tempOld[i] = inTemp;
    tempNew[i] = inTemp;
  }
}

// dtor, delete two 3d temp array
TempClass::~TempClass()
{
  delete [] tempOld;
  delete [] tempNew;
}

// function that update the new temperature using FVM
// based on old temperature, loaded time, laser, material
// use FVM and openMP
// after that update the tempOld with tempNew
void TempClass::calTempHistFVM(const string fileName)
{
  // initialize the temp
  initializeTemp(tempAmbient);

  // set the fileName for the pointGroup object
  thisPointGroup->setFileName(fileName);
  cout << *thisPointGroup << endl;

  int i,j,k,t; // loop variable
  double xLaserPos, yLaserPos;
  double timeCurrent = 0.;
  double beginTime = omp_get_wtime();
  bool laserIsActive; // whether laser is still active

  // set the number of threads used in parallel computing
  omp_set_num_threads(16);

  for (t = 0; t < numTimeStep; t++)
  {
    timeCurrent = t * timeStep;

    // update the laser position and state (if active)
    laserIsActive =
         thisLaser->calLaserPosition(
         xLaserPos,
         yLaserPos,
         timeCurrent
         );


    // update the BC based on the time
#pragma omp parallel for collapse(3)
    for (i = 0; i < xSize; i++) {
      for (j = 0; j < ySize; j++) {
        for (k = 0; k < zSize; k++) {
          int yzSize = ySize * zSize;
          int iyzSize = i * yzSize;
          int jzSize = j * zSize;
          double tempLocal = 0.;

          // BC when x = 0
          if (i != 0)
          {
            tempLocal +=
                 thermalDiff * (tempOld[iyzSize - yzSize + jzSize + k]
                 - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }

          // BC when x = xLen
          if (i != xSize - 1) {
            tempLocal +=
                 thermalDiff * (tempOld[iyzSize + yzSize + jzSize + k]
                 - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }

          // BC when y = 0
          if (j != 0)
          {
            tempLocal +=
                 thermalDiff * (tempOld[iyzSize + jzSize - zSize + k]
                 - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }

          // BC when y = yLen
          if (j != ySize - 1)
          {
            tempLocal +=
                 thermalDiff * (tempOld[iyzSize + jzSize + zSize + k]
                 - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }

          // BC when z = zLen
          if (k != zSize - 1)
          {
            tempLocal +=
                 thermalDiff * (tempOld[iyzSize + jzSize + k + 1]
                 - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }

          // special BC relative to heat flux when z = 0
          if (k == 0)
          {
            double xDistToLaser = xc[i] - xLaserPos;
            double yDistToLaser = yc[j] - yLaserPos;
            // if within laser spot and laser is still active
            if (
                 xDistToLaser * xDistToLaser + yDistToLaser * yDistToLaser
                 <= beamRadius * beamRadius &&
                 laserIsActive
                 )
            {
              tempLocal += (
                   ConstantFlux1 * ConstantFlux2 *
                   exp(-3 * (xDistToLaser * xDistToLaser +
                   yDistToLaser * yDistToLaser) /
                   (beamRadius * beamRadius))
              );
            }
            // if on the top boundary but outside the laser spot
            // there will be zero heat flux, therefore tempLocal += 0
          }
          else // if not on the top surfaces
          {
            tempLocal += (
                 thermalDiff * (tempOld[iyzSize + jzSize + k - 1]
                 - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep
                 );
          }
          // calculate the temperature at current position
          tempNew[iyzSize + jzSize + k] =
               tempLocal + tempOld[iyzSize + jzSize + k];
        }
      }
    }

#pragma omp barrier
    // write into dat
    writeTemp2D(timeCurrent, t, fileName);
    thisPointGroup->writePointsHist();

    // update the temperature
#pragma omp parallel for
    for (i = 0; i < xSize * ySize * zSize; i++)
    {
      tempOld[i] = tempNew[i];
    }
#pragma omp barrier


    //cout << "time step: " << t << " calculated finished." << endl;
  }

  double endTime = omp_get_wtime();

  runTime = endTime - beginTime;

  cout << "Calculation finished!" << endl;

  printRunTime();

  // initialize the temp
  initializeTemp(tempAmbient);
}

// function that calculate thermal history using FVM method
// all the surface except the bottom face is set to be
// convective and radiation boundary condition
// the bottom surface will still be isolated boundary condition
// after calculating the thermal history, will still automatically
// write temperature into .dat file with fileName
void TempClass::calTempHistFVMwithConvAndRadiaBC(const string fileName)
{
  // initialize the temp
  initializeTemp(tempAmbient);

  // set the fileName for the pointGroup object
  thisPointGroup->setFileName(fileName);
  cout << *thisPointGroup << endl;

  int i,j,k,t; // loop variable
  double xLaserPos, yLaserPos;
  double timeCurrent = 0.;
  double beginTime = omp_get_wtime();
  bool laserIsActive; // whether laser is still active

  // set the number of threads used in parallel computing
  omp_set_num_threads(14);

  for (t = 0; t < numTimeStep; t++) {
    timeCurrent = t * timeStep;

    // update the laser position and state (if active)
    laserIsActive =
         thisLaser->calLaserPosition(
              xLaserPos,
              yLaserPos,
              timeCurrent
         );


    // update the BC based on the time
#pragma omp parallel for collapse(3)
    for (i = 0; i < xSize; i++) {
      for (j = 0; j < ySize; j++) {
        for (k = 0; k < zSize; k++) {
          int yzSize = ySize * zSize;
          int iyzSize = i * yzSize;
          int jzSize = j * zSize;
          double tempLocal = 0.;

          // check whether local position is west face: x = meshSize/2
          // if so, apply the convective and radiation BD
          // using equivalent temperature of the boundary
          if (i != 0) {
            tempLocal +=
                 thermalDiff * (tempOld[iyzSize - yzSize + jzSize + k]
                                - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          } else {
            tempLocal +=
                 thermalDiff * (
                      calcEquTempBound(tempOld[iyzSize + jzSize + k])
                      - tempOld[iyzSize + jzSize + k])
                      / (meshSize * meshSize) * timeStep;
          }

          // check whether local position is east surface: x = xlen - meshSize/2
          // if so, apply the equivalent boundary temperature
          if (i != xSize - 1) {
            tempLocal +=
                 thermalDiff * (tempOld[iyzSize + yzSize + jzSize + k]
                                - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          } else {
            tempLocal +=
                 thermalDiff * (
                      calcEquTempBound(tempOld[iyzSize + jzSize + k])
                      - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }


          // check whether local position is east surface: y = meshSize/2
          // if so, apply the equivalent boundary temperature
          if (j != 0) {
            tempLocal +=
                 thermalDiff * (tempOld[iyzSize + jzSize - zSize + k]
                                - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }else {
            tempLocal +=
                 thermalDiff * (
                      calcEquTempBound(tempOld[iyzSize + jzSize + k])
                      - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }

          // check whether local position is east surface: y = ylen - meshSize/2
          // if so, apply the equivalent boundary temperature
          if (j != ySize - 1) {
            tempLocal +=
                 thermalDiff * (tempOld[iyzSize + jzSize + zSize + k]
                                - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }else{
            tempLocal +=
                 thermalDiff * (
                      calcEquTempBound(tempOld[iyzSize + jzSize + k])
                      - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }

          // check whether local position is at the bottom: x = 0
          // if so, apply the isolated BD: tempLocal += 0
          if (k != zSize - 1) {
            tempLocal +=
                 thermalDiff * (tempOld[iyzSize + jzSize + k + 1]
                                - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }else{
            tempLocal +=
                 thermalDiff * (
                      calcEquTempBound(tempOld[iyzSize + jzSize + k])
                      - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep;
          }

          // special BC relative to heat flux when z = 0
          // for the top surface
          // if the local point position is within the heat source,
          // it will have heat flux boundary.
          // otherwise, it will be convection and radiation boundary
          if (k == 0) {
            double xDistToLaser = xc[i] - xLaserPos;
            double yDistToLaser = yc[j] - yLaserPos;
            // if within laser spot and laser is still active
            if (
                 xDistToLaser * xDistToLaser + yDistToLaser * yDistToLaser
                 <= beamRadius * beamRadius &&
                 laserIsActive
                 ) {
              tempLocal += (
                   ConstantFlux1 * ConstantFlux2 *
                   exp(-3 * (xDistToLaser * xDistToLaser +
                             yDistToLaser * yDistToLaser) /
                       (beamRadius * beamRadius))
              );
            }
            // if on the top boundary but outside the laser spot
            // there will still be convective and radiation effect
            else {
              tempLocal +=
                   thermalDiff * (
                        calcEquTempBound(tempOld[iyzSize + jzSize + k])
                        - tempOld[iyzSize + jzSize + k])
                   / (meshSize * meshSize) * timeStep;
            }
          } else // if not on the top surface, normal FVM formula
          {
            tempLocal += (
                 thermalDiff * (tempOld[iyzSize + jzSize + k - 1]
                                - tempOld[iyzSize + jzSize + k])
                 / (meshSize * meshSize) * timeStep
            );
          }
          // calculate the temperature at current position
          tempNew[iyzSize + jzSize + k] =
               tempLocal + tempOld[iyzSize + jzSize + k];
        }
      }
    }

#pragma omp barrier
    // write into dat
    writeTemp2D(timeCurrent, t, fileName);
    thisPointGroup->writePointsHist();

    // update the temperature
#pragma omp parallel for
    for (i = 0; i < xSize * ySize * zSize; i++) {
      tempOld[i] = tempNew[i];
    }
#pragma omp barrier

  }
    //cout << "time step: " << t << " calculated finished." << endl;
}

// function calculating the boundary temperature transformed from
// convective and radiation boundary condition
// solve myFunc using brent method (a hybrid root-finding method combining
// bisection, newton and scant method) provided by gsl library
// coefficient a, b of myFunc could be shared by all the myFunc
// for coefficient c = c0 + c1 * temp_p, c0 could be pre-calculated
double TempClass::calcEquTempBound(
     double temp_p // temperature at the point next to the boudnary
     )
{
  gsl_function F;

  struct myFuncParams params ={
       coeff_a,
       coeff_b,
       coeff_c0 + coeff_c1 * temp_p
       };

  F.function =&myFunc;
  F.params = &params;

  int status;
  int iter = 0, max_iter = 100;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0, r_expected = 3998;
  double x_lo = tempAmbient-1, x_hi = temp_p+1;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

//    printf ("using %s method\n",
//            gsl_root_fsolver_name (s));
//
//    printf ("%5s [%9s, %9s] %9s %10s %9s\n",
//            "iter", "lower", "upper", "root",
//            "err", "err(est)");

  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi,
                                     0, 0.01);

//    if (status == GSL_SUCCESS)
//    {
//      printf ("Converged result:\n");
//      printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n",
//              iter, x_lo, x_hi,
//              r, r - r_expected,
//              x_hi - x_lo);
//    }
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  // if brent method has not meet condition, print the situation
  if (iter == max_iter)
  {
    cout << "Bad condition:" << endl
         << "  temp_p: " << temp_p
         << "  temp_b: " << (x_lo + x_hi) / 2.0;
  }

  gsl_root_fsolver_free (s);

  //cout << temp_p - r << endl;

  return r;
}

// value of my function
double TempClass::myFunc (double x, void *params)
{
  struct myFuncParams *p
       = (struct myFuncParams *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  return (a * x * x * x + b) * x + c;
}

// first derivative of my function
double TempClass::myFunc_deriv (double x, void *params)
{
  struct myFuncParams *p
       = (struct myFuncParams *) params;

  double a = p->a;
  double b = p->b;

  return 4.0 * a * x * x * x + b;
}

// both the value and the first derivative of my function
void TempClass::myFunc_fdf (double x, void *params,
                 double *y, double *dy)
{
  struct myFuncParams *p
       = (struct myFuncParams *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  *y = (a * x * x * x + b) * x + c;
  *dy = 4.0 * a * x * x * x + b;
}


// function that calculate thermal history using analytical method
// based on loaded time, laser, material
// for each time step, decide the abssica of nodes,
// calculate function value at each nodes, times weights and sum together
// will automatically write temperature into .dat file with fileName
void TempClass::calTempHistAnalytical(const string fileName)
{
  // initialize the temp
  initializeTemp(tempAmbient);

  // set the file name for the pointGroup object
  thisPointGroup->setFileName(fileName);
  cout << *thisPointGroup << endl;

  int i, j, k, t, l; // loop variable
  //double xLaserPos, yLaserPos; // position of laser
  double timeCurrent = 0.; // current time

  double beginTime = omp_get_wtime();
  //bool laserIsActive; // whether laser is still active
  int numNodes = 16; // number of nodes used for calculating

  // set the number of threads used in parallel computing
  omp_set_num_threads(16);


  for (t = 0; t <= numTimeStep; t++)
  {
    timeCurrent = t * timeStep;

    gsl_integration_glfixed_table *table
         = gsl_integration_glfixed_table_alloc(numNodes);

    // for the same time step, all the cube domain share the same
    // 1. abssicaList
    // 2. weightList
    // 2. innerConstant for each node
    // 3. phiX for each node
    // 4. 1/sqrt(phiX * phiY * phiZ)
    // therefore, the following loop precalculate these value
    // before loop w.r.t 3D space

    // allocate these shared value
    double abssicaList[numNodes];
    double weightList[numNodes];
    double innerConstant[numNodes];
    double phiX[numNodes];
    double invSqrtPhiX3[numNodes];
    double xLaserPosList[numNodes];
    double yLaserPosList[numNodes];

    // precalculate
    for (l = 0; l < numNodes; l++)
    {
      // based on the active time of laser, decide the abssica and weight
      // routine to retrieve one abssica and weight
      gsl_integration_glfixed_point(
           0.,
           min(scanTime, timeCurrent),
           l,
           abssicaList + l, // pointer arithematic
           weightList + l, // pointer arithematic
           table
           );

      innerConstant[l] = 12. * thermalDiff * (timeCurrent - abssicaList[l]);
      phiX[l] = innerConstant[l] + (beamRadius * beamRadius);
      invSqrtPhiX3[l] = 1. / sqrt(phiX[l] * phiX[l] * phiX[l]);
      // update the laser position and state (if active)
      thisLaser->calLaserPosition(
           xLaserPosList[l],
           yLaserPosList[l],
           abssicaList[l]
           );
    }

//#pragma omp parallel for collapse(3) private(weightList, invSqrtPhiX3, abssicaList)
    // loop wrt 3d position
    for (i = 0; i < xSize; i++)
    {
      for (j = 0; j < ySize; j++)
      {
        for (k = 0; k < zSize; k++)
        {
          double tempLocal = 0.; // initial local temp

          // calculate the approx numerical integration value
          for (l = 0; l < numNodes; l++)
          {
            tempLocal += weightList[l] * invSqrtPhiX3[l] *
                 exp(-3 * ((xc[i] - xLaserPosList[l]) * (xc[i] - xLaserPosList[l]) / phiX[l] +
                 (yc[j] - yLaserPosList[l]) * (yc[j] - yLaserPosList[l]) / phiX[l] +
                 (zc[k] * zc[k]) / phiX[l]));

          }

          // update the current temperature point
          tempOld[i * ySize * zSize + j * zSize + k] =
               tempLocal * ConstantAnaly + tempAmbient;
        }
      }
    }

//#pragma omp barrier
    // write into .dat
    writeTemp2D(timeCurrent, t, fileName);
    thisPointGroup->writePointsHist();

    // free table at each time step
    gsl_integration_glfixed_table_free(table);

    // cout << "time step: " << t << " calculated finished." << endl;
  }

  double endTime = omp_get_wtime();

  runTime = endTime - beginTime;

  printRunTime();

  // initialize the temp
  initializeTemp(tempAmbient);
}

// function that calculate thermal history using analytical method
// with adaptive numerical integration
void TempClass::calTempHistAnalyticalAdaptive(const string fileName)
{
  // initialize the temp
  initializeTemp(0.);

  // set the file name for the pointGroup object
  thisPointGroup->setFileName(fileName);
  cout << *thisPointGroup << endl;

  int i, j, k, t, l; // loop variables

  double timeCurrent = 0.; // current time

  double beginTime = omp_get_wtime();

  int numNodes = 0; // number of integration nodes for each sub-interval

  omp_set_num_threads(14);

  // non-dimensional parameters

  // nd total time of next path segment
  double timeTotal_nd = totalTime * thermalDiff / (beamRadius * beamRadius);
  // nd current time (s)
  double timeCurrent_nd = timeStep * numTimeStep
       * thermalDiff / (beamRadius * beamRadius);
  // nd conduction time (s'')
  double timeDiff_nd;
  // nd scan speed
  double speed_nd;
  // nd preliminary segment size of time
  double prelimTimeSegSize_nd;
  // nd segment size of time
  double timeSegSize_nd;


  // loop w.r.t conduction time, which is the inverse of scan time
  for (t = 0; t <= numTimeStep; t++)
  {
    // update current time and nd current time
    timeCurrent = t * timeStep;
    timeCurrent_nd = timeCurrent * thermalDiff / (beamRadius * beamRadius);


    // non-dimensional diffusion time
    timeDiff_nd = 0.;

    // loop over conduction time to figure out abssicas of nodes and weights
    while (timeDiff_nd < timeCurrent_nd)
    {

      // for the current sub interval of non-dimensional time
      // calculate the velocity of current sub interval
      if (thisLaser->calLaserSpeed(
           speed_nd,
           (timeCurrent_nd - timeDiff_nd) * beamRadius * beamRadius / thermalDiff
           ))
      {
        speed_nd = speed_nd * beamRadius / thermalDiff;
      }
      else // if currently laser is not active, be more convectional
      {
        speed_nd = 0.58;
      }

      //cout << "speed_nd: " << speed_nd << endl;

      // preliminary segment size
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

      // do not skip path segments
      if (timeDiff_nd + timeSegSize_nd > timeCurrent_nd)
      {
        timeSegSize_nd = timeCurrent_nd - timeDiff_nd;
      }

      // determine quadrature nodes and weights for current integration segment
      // and add to list of all nodes and weights

      gsl_integration_glfixed_table *table
           = gsl_integration_glfixed_table_alloc(numNodes);

      // for the same time step, all the cube domain share the same
      // 1. abssicaList
      // 2. weightList
      // 3. innerConstantList
      // 4. xy position of laser beam
      // therefore, the following loop precalculate these value
      // before loop w.r.t 3D space

      // allocate these shared value
      double* abssicaList = new double[numNodes];
      double* weightList = new double[numNodes];

      double* innerConst1 = new double[numNodes];
      double* innerConst2 = new double[numNodes];
      double* innerConst3 = new double[numNodes];

      double* xLaserPosList = new double[numNodes];
      double* yLaserPosList = new double[numNodes];

      // array that store the status of laser at each integration nodes
      // if laser is active, element will be 1., if off, will be 0.
      double* laserIsOn = new double[numNodes];


      // precalculate
      for (l = 0; l < numNodes; l++)
      {
        // based on the active time of laser, decide the abssica and weight
        // routine to retrieve one abssica and weight
        gsl_integration_glfixed_point(
             timeDiff_nd,
             timeDiff_nd + timeSegSize_nd,
             l,
             abssicaList + l, // pointer arithematic
             weightList + l, // pointer arithematic
             table
             );

        // 3 constants lists
        innerConst2[l] = 12. * abssicaList[l] + 1.;
        innerConst3[l] = 12. * abssicaList[l];
        innerConst1[l] = 1. / (innerConst2[l] * sqrt(innerConst3[l])) ;

        // update the laser position and state
        if (!thisLaser->calLaserPosition(
             xLaserPosList[l],
             yLaserPosList[l],
             (timeCurrent_nd - abssicaList[l]) * beamRadius * beamRadius / thermalDiff
             ))
        {
          // if laser is off, zeros thermal contribution
          laserIsOn[l] = 0.;
//          cout << "time is out of range for:" << endl
//              << "abssicaList[l]: " << abssicaList[l] << endl
//              << "equ_time:" << (timeCurrent_nd - abssicaList[l]) * beamRadius * beamRadius / thermalDiff << endl;
        }
        else
        {
          // if laser is on, full thermal contribution
          laserIsOn[l] = 1.;
        }
      }
//      if (t == 2)
//      {
//        cout << "Time: " << timeCurrent_nd  << "/ "
//            << "s'' sub interval: [" << timeDiff_nd << " ,"
//            << timeDiff_nd + timeSegSize_nd << "]" << endl;
//
//        for (l = 0; l < numNodes; l++)
//        {
//          cout << "equ_time:" << (timeCurrent_nd - abssicaList[l]) * beamRadius * beamRadius / thermalDiff << endl;
//          cout << "  [xPos: " << xLaserPosList[l] << " | "
//              << "yPos: " << yLaserPosList[l] << "] " << endl;
////            cout << "timeDiff_nd_currentNode: " << abssicaList[l] << "  ";
////            cout << "  prelimTimeSegSize_nd: " << prelimTimeSegSize_nd;
////            cout << "  numNodes: " << numNodes << " || ";
//        }
//      }


      // loop wrt 3d position to
      for (i = 0; i < xSize; i++)
      {
        for (j = 0; j < ySize; j++)
        {
          for (k = 0; k < zSize; k++)
          {
            double tempLocal = 0.; // summation of thermal contribution of this sub interval

            // calculate the approx numerical integration value
            for (l = 0; l < numNodes; l++)
            {
              tempLocal += laserIsOn[l] * weightList[l] * innerConst1[l] *
                           exp(-3. *
                           (((xc[i] - xLaserPosList[l]) / beamRadius
                           * (xc[i] - xLaserPosList[l]) / beamRadius +
                           (yc[j] - yLaserPosList[l]) / beamRadius
                           * (yc[j] - yLaserPosList[l]) / beamRadius) / innerConst2[l]  +
                           (zc[k] * zc[k]) / beamRadius / beamRadius / innerConst3[l]));
            }
            tempOld[i * ySize * zSize + j * zSize + k] +=
                 phyConst_nd * tempLocal;

          }
        }
      }
      delete [] abssicaList;
      delete [] weightList;
      delete [] innerConst1;
      delete [] innerConst2;
      delete [] innerConst3;
      delete [] xLaserPosList;
      delete [] yLaserPosList;
      delete [] laserIsOn;

      // free table at each time step
      gsl_integration_glfixed_table_free(table);

      // increment conduction time
      timeDiff_nd += timeSegSize_nd;
    }

    // add the ambient temperature
    for (i = 0; i < xSize * ySize * zSize; i++)
    {
      tempOld[i] += tempAmbient;
    }


//#pragma omp barrier
    // write into .dat
    writeTemp2D(timeCurrent, t, fileName);
    // thisPointGroup->writePointsHist();
    initializeTemp(0.);

    // cout << "time step: " << t << " calculated finished." << endl;

  }

  double endTime = omp_get_wtime();

  runTime = endTime - beginTime;

  printRunTime();

  // initialize the temp
  initializeTemp(0.);

  // calculation finished
  cout << "Calculation finished!" << endl;
}

// function that calculate thermal history using analytical method
// and consider the edge through folding the temperature filed outside the
// boundary back into the field
// will automatically write temperature into .dat file with fileName
void TempClass::calTempHistAnalyMirror(const string fileName)
{
  // initialize the temp
  initializeTemp(tempAmbient);

  // set the file name for the pointGroup object
  thisPointGroup->setFileName(fileName);
  cout << *thisPointGroup << endl;

  int i, j, k, t, l; // loop variable
  //double xLaserPos, yLaserPos; // position of laser
  double timeCurrent = 0.; // current time

  double beginTime = omp_get_wtime();
  //bool laserIsActive; // whether laser is still active
  int numNodes = 16; // number of nodes used for calculating

  // set the number of threads used in parallel computing
  omp_set_num_threads(16);


  for (t = 0; t <= numTimeStep; t++)
  {
    timeCurrent = t * timeStep;

    gsl_integration_glfixed_table *table
         = gsl_integration_glfixed_table_alloc(numNodes);

    // for the same time step, all the cube domain share the same
    // 1. abssicaList
    // 2. weightList
    // 2. innerConstant for each node
    // 3. phiX for each node
    // 4. 1/sqrt(phiX * phiY * phiZ)
    // therefore, the following loop precalculate these value
    // before loop w.r.t 3D space

    // allocate these shared value
    double abssicaList[numNodes];
    double weightList[numNodes];
    double innerConstant[numNodes];
    double phiX[numNodes];
    double invSqrtPhiX3[numNodes];
    double xLaserPosList[numNodes];
    double yLaserPosList[numNodes];

    // precalculate
    for (l = 0; l < numNodes; l++)
    {
      // based on the active time of laser, decide the abssica and weight
      // routine to retrieve one abssica and weight
      gsl_integration_glfixed_point(
           0.,
           min(scanTime, timeCurrent),
           l,
           abssicaList + l, // pointer arithematic
           weightList + l, // pointer arithematic
           table
           );

      innerConstant[l] = 12. * thermalDiff * (timeCurrent - abssicaList[l]);
      phiX[l] = innerConstant[l] + (beamRadius * beamRadius);
      invSqrtPhiX3[l] = 1. / sqrt(phiX[l] * phiX[l] * phiX[l]);
      // update the laser position and state (if active)
      thisLaser->calLaserPosition(
           xLaserPosList[l],
           yLaserPosList[l],
           abssicaList[l]
           );
    }

//#pragma omp parallel for collapse(3) private(weightList, invSqrtPhiX3, abssicaList)
    // loop wrt 3d position
    for (i = 0; i < xSize; i++)
    {
      for (j = 0; j < ySize; j++)
      {
        for (k = 0; k < zSize; k++)
        {
          double tempLocal = 0.; // initial local temp
          double tempLocalX = 0.; // extra place
          double tempLocalY = 0.; // extra place

         // calculate the approx numerical integration value
          for (l = 0; l < numNodes; l++)
          {
            tempLocal += weightList[l] * invSqrtPhiX3[l] *
                         exp(-3 * ((xc[i] - xLaserPosList[l]) * (xc[i] - xLaserPosList[l]) / phiX[l] +
                                   (yc[j] - yLaserPosList[l]) * (yc[j] - yLaserPosList[l]) / phiX[l] +
                                   (zc[k] * zc[k]) / phiX[l]));
            tempLocalX += weightList[l] * invSqrtPhiX3[l] *
                         exp(-3 * ((xc[i] - (xLen) - xLaserPosList[l]) * (xc[i] - (xLen)- xLaserPosList[l]) / phiX[l] +
                                   (yc[j] - yLaserPosList[l]) * (yc[j] - yLaserPosList[l]) / phiX[l] +
                                   (zc[k] * zc[k]) / phiX[l]));
            tempLocalY += weightList[l] * invSqrtPhiX3[l] *
                         exp(-3 * ((xc[i] - xLaserPosList[l]) * (xc[i] - xLaserPosList[l]) / phiX[l] +
                                   (yc[j] - yLaserPosList[l] - (yLen)) * (yc[j] - yLaserPosList[l] - (yLen)) / phiX[l] +
                                   (zc[k] * zc[k]) / phiX[l]));
          }
          tempOld[i * ySize * zSize + j * zSize + k] =
               tempLocal * ConstantAnaly + tempAmbient;

          tempOldX[i * ySize * zSize + j * zSize + k] =
               tempLocalX * ConstantAnaly;

          tempOldY[i * ySize * zSize + j * zSize + k] =
               tempLocalY * ConstantAnaly;
        }
      }
    }

    // update every point temp of current time step
    for (i = 0; i < xSize; i++) {
      for (j = 0; j < ySize; j++) {
        for (k = 0; k < zSize; k++) {
          tempOld[i * ySize * zSize + j * zSize + k] +=
               tempOldX[(xSize - 1 - i) * ySize * zSize + j * zSize + k] +
               tempOldY[i * ySize * zSize + (ySize - 1 - j) * zSize + k];
        }
      }
    }

    // after calculating the whole field, calculate the extra field to be folded
    // based on the analytical solution.

    // define the extra field parameters

//#pragma omp barrier
    // write into .dat
    writeTemp2D(timeCurrent, t, fileName);
    thisPointGroup->writePointsHist();

    // free table at each time step
    gsl_integration_glfixed_table_free(table);

    // cout << "time step: " << t << " calculated finished." << endl;
  }

  double endTime = omp_get_wtime();

  runTime = endTime - beginTime;

  printRunTime();

  // initialize the temp
  initializeTemp(tempAmbient);
}


// function that write the current top layer of tempNew into dat file
// will be called after updating temp field
void TempClass::writeTemp2D(
     double timeCurrent,  // current time
     int timeStepCurrent, // current time step
     string fileName
     )
{
  int timeStepPlot = int (numTimeStep / (numPlot));

  // judge whether it is the time to writing
  if ((timeStepCurrent % timeStepPlot) == 0.)
  {
    // write 2D plane of temperature
    ofstream outFile;
    int indexPlotCurrent = (int(timeStepCurrent / timeStepPlot));
    fileName = "./tecplot_2D/" + fileName + to_string(indexPlotCurrent) + ".dat";
    outFile.open(fileName.c_str());

    if (outFile.fail())
    {
      cout << "Fail to write data into " << fileName << endl;
      exit(3);
    }
    else
    {
      outFile << "VARIABLES=\"X\" \"Y\" \"Temp\" "<< endl;
      outFile << "zone t=\"" << indexPlotCurrent << "\", i=" << xSize
           << ", j=" << ySize <<", f=point" << endl;
      outFile << "solutiontime=" << indexPlotCurrent << endl;
      cout << "write output at time step: " << timeStepCurrent
           << " , solutiontime: " << timeCurrent
           << " into " << fileName << endl;

      for (int i = 0; i < xSize; i++)
      {
        for (int j = 0; j < ySize; j++)
        {
          outFile << xc[i] << " " << yc[j] << " "
               << tempOld[i * ySize * zSize + j * zSize] << endl;
        }
      }
    }
    outFile.close();
    cout << "Write finished!" << endl;

    // store the point's temperature
    thisPointGroup->storePointsOneTime(
         indexPlotCurrent,
         timeCurrent,
         tempOld
         );
    cout << "Store point Hist finished!" << endl;
  }
}

// function that print the run time
void TempClass::printRunTime()
{
  cout << "The run time for this case is: " << runTime << "s" << endl;
}

// function of setting numPlot
void TempClass::setNumPlot(int inNumPlot)
{
  numPlot = inNumPlot;
}
