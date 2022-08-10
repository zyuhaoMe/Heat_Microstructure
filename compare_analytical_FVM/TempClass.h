/* 
 * Author: zyuhao 
 * Date: 2022-04-11
 * Content:  Head file of TempClass, compute temperature field
 * contains old temperature 3d array and new temperature 3d array
 * defines relative functions of calculating temperature filed
 * using FVM or analytical method
 * defines functions that write temperature field into .dat file
 * 
 */

//

#ifndef COMPARE_ANALYTICAL_FVM_TEMPCLASS_H
#define COMPARE_ANALYTICAL_FVM_TEMPCLASS_H

#include "LaserClass.h"
#include "CubeClass.h"
#include "MaterialClass.h"
#include "TimeClass.h"
#include "PointGroupHistClass.h"
#include "constants.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <vector>
#include <omp.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

class TempClass {

  private:
    double coeff_a; // coefficient a of myFunc
    double coeff_b; // coefficient b of myFunc
    double coeff_c0; // c = c0 + c1 * temp_p
    double coeff_c1;

    double SBConstant; // stefan-bolzmann constant (bad abbreviation in Chinese)

    // function and its derivative whose root are temperature at boundary
    // transformed from convective and radiation boundary condition

    // struct to pack the parameters function
    struct myFuncParams
    {
      double a, b, c;
    };

    // value of my function
    static double myFunc (double x, void *params);

    // first derivative of my function
    static double myFunc_deriv (double x, void *params);

    // both the value and the first derivative of my function
    static void myFunc_fdf (double x, void *params,
                     double *y, double *dy);

    // function calculating the boundary temperature transformed from
    // convective and radiation boundary condition
    // using brent method (a hybrid root-finding method combining
    // bisection, newton and scant method) provided by gsl library
    // coefficient a, b of myFunc could be shared by all the myFunc
    // for coefficient c = c0 + c1 * temp_p, c0 could be pre-calculated
    double calcEquTempBound(
         double temp_p // temperature at the point next to the boudnary
         );

  public:
    double tempAmbient; // ambient temperature

    double* tempOld; // dynamic 3d array of old temperature
    double* tempOldX;
    double* tempOldY;
    double* tempNew; // dynamic 3d array of new temperature

    // To increase the memory locality, therefore boost performance
    // I create copy of frequently-used variables within the temp object
    // in FVM and ana equ from other objects locally and set them in ctor.

    // Honestly, it is an ugly way, more future optimization could be made
    // if it is necessary, or the whole program architecture could be
    // re-designed (from composition to inherit etc).
    // It is totally up to you.
    // Please email me: zhouyuhao2020@gmail.com if you are struggling
    // understanding the code. I will be glad to walk you through
    // the script line by line. Good luck!

    double meshSize; // meshSize
    double thermalDiff; // diffusivity
    double timeStep; // time step size
    int numTimeStep; // total number of time step
    double beamRadius; // beam radius of laser
    double scanTime; // time when laser is active
    double totalTime; // total time of simulation

    LaserClass* thisLaser; // will use the function of laser power

    MaterialClass* thisMaterial;
    // will use absorp, density, specHeat, convecCoeff and emissivity

    PointGroupHistClass* thisPointGroup; // used for writing points' thermalHist

    double* xc;
    double* yc;
    double* zc;
    int xSize, ySize, zSize; // number of mesh along 3d dimensions
    double xLen, yLen, zLen;
    double xStart, yStart, zStart;

    double ConstantFlux1; // physical constant in heat flux BC
    double ConstantFlux2; // physical constant in heat flux BC of FVM

    double ConstantAnaly; // physical constant used in the analytical model

    double phyConst_nd; // physical constant used in non-dimensional form

    int numPlot; // number of plots to be written
    
    int numPlotPoint;

    double runTime; // time of running the case

    int numThreads; // number of threads used for parallel computing

    // default ctor, load data from dependent object to and store it locally
    // to prevent extra work
    TempClass(
         MaterialClass* myMaterial,
         CubeClass* myCube,
         LaserClass* myLaser,
         TimeClass* myTimePlan,
         PointGroupHistClass* myPointGroup
         );

    // function that initialize the tempOld and tempNew to the ambient temp
    void initializeTemp(double inTemp);

    // dtor, delete two 3d temp array
    ~TempClass();

    // function that calculate thermal history using FVM method
    // all the boundary condition is set to be isolated condition
    // which mean zeros heat flux on all the surfaces
    // except the tiny area within the active heat source
    // properties will be based on loaded time, laser, material
    // will automatically write temperature into .dat file with fileName
    void calTempHistFVM(const string fileName);

    // function that calculate thermal history using FVM method
    // all the surface except the bottom face is set to be
    // convective and radiation boundary condition
    // the bottom surface will still be isolated boundary condition
    // due to the forth-power term in radiation equation
    // an iterative root finding algorithm is used
    // to find the equivalent boundary temperature due to
    // heat convection and radiation
    // after calculating the thermal history, will still automatically
    // write temperature into .dat file with fileName
    void calTempHistFVMwithConvAndRadiaBC(const string fileName);

    // function that calculate thermal history using analytical method
    // based on loaded time, laser, material
    // will automatically write temperature into .dat file with fileName
    // using gaussian legendre integration of fixed number nodes
    void calTempHistAnalytical(const string fileName);

    // function that calculate thermal history using analytical method
    // with adaptive numerical integration
    void calTempHistAnalyticalAdaptive(const string fileName);

    // function that calculate thermal history using analytical method
    // and consider the edge through folding the temperature filed outside the
    // boundary back into the field
    // will automatically write temperature into .dat file with fileName
    void calTempHistAnalyMirror(const string fileName);

    // function that write the current tempNew into the .dat file
    // will be called after updating Temp field
    void writeTemp2D(
         double timeCurrent,
         int timeStepCurrent,
         string fileName
         );


    // function that print the run time
    void printRunTime();


    // function that set the num of plot, which is the number of .dat file
    // to be store during calculation
    void setNumPlot(int inNumPlot);

};


#endif //COMPARE_ANALYTICAL_FVM_TEMPCLASS_H
