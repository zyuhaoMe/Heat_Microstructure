// This file stores functions about fixed-nodes Gaussian-Legendre integration
// using GSL library directly, for detail documentation, please refer to the following website:
// gnu.org/software/gsl/doc/html/integration.html#c.gsl_integration_glfixed
// Numerical integration is the "mechanism" of heat transfer.

#include "CA_heat.h"
#include <gsl/gsl_integration.h>

double computeThermalFieldWithGaussianLegendre(
        int numNodes,
        int xsize,
        int ysize,
        int zsize,
        double* xc,
        double* yc,
        double* zc,
        double dtca,
        double timeca,
        double ***temp,
        double phyConst,
        double laserTime) {
    int i, j, k; // loop variable

    // determine the Gauss-Legendre abscissae and weights necessary for an n-point fixed order integration scheme
    gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(numNodes);


    // at each time step, calculate the G-L point and relative weight
    double abssicaList[numNodes];
    double weightList[numNodes];

    if (timeca <= laserTime) { // if the laser source is still on
        for (i = 0; i < numNodes; i++) {
            double *xi = abssicaList + i;
            double *wi = weightList + i;
            gsl_integration_glfixed_point(0.0, timeca, i, xi, wi, t);
        }
    }else{ // if the laser is off, it is during the cooling time
        for (i = 0; i < numNodes; i++) {
            double *xi = abssicaList + i;
            double *wi = weightList + i;
            gsl_integration_glfixed_point(0.0, laserTime, i, xi, wi, t);
        }
    }
    // for each point of interest, calculate temperature by summing heat contribution of each node
    for (i = 0; i < xsize; i++) {
        for (j = 0; j < ysize; j++) {
            for (k = 0; k < zsize; k++) {
                temp[i][j][k] = phyConst * summationOfGaussianLegendre(
                        xc[i],
                        yc[j],
                        zc[k],
                        timeca,
                        abssicaList,
                        weightList,
                        numNodes);
            }
        }
    }
    gsl_integration_glfixed_table_free(t); // free the memory associated with the table
}



// evaluate the numerical integration by summing wi*f(xi) over i for each point of interest
double summationOfGaussianLegendre(
        double xVal,
        double yVal,
        double zVal,
        double timeca,
        double* abssicaList,
        double* weightList,
        int numNodes)
{
    int i;
    double tempOnePoint = 0.0;
    for (i = 0; i < numNodes; i++)
    {
        tempOnePoint += weightList[i] * calculateIntegrandVal(abssicaList[i], timeca, xVal, yVal, zVal) + 300.0;
    }
    return tempOnePoint;
}


