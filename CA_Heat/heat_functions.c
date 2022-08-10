//
// Created by zyuhao on 2022-01-18.
// This file contains functions about heat transfer of moving heat source,
// such as the function of defining laser scanning path and scanning velocities.
// Therefore, this file is about the "policy" of heat transfer and is seperated
// from the numerical integration implemention file.

// functions used in analytical heat transfer computation

#include "CA_heat.h"

//void import_data(double* X, double* Y, double* L, double* D, double* T, int nrows) {
//    // read_csv
//    char filename[] = "DAQ_T50us_MS_Cube_fused_P04_layer0101.csv";
//    FILE* the_file;
//    int err;
//    err = fopen(&the_file, filename, "r");
//    if (err != 0) {
//        perror("Unable to open the file.");
//    }
//    char line[nrows];
//    int row = 0;
//    while ((fgets(line, sizeof(line), the_file)) && (row < nrows)) {
//        char* token;
//        token = strtok(line, ",");
//        X[row] = atof(token);
//        token = strtok(NULL, ",");
//        Y[row] = atof(token);
//        token = strtok(NULL, ",");
//        L[row] = atof(token);
//        token = strtok(NULL, ",");
//        D[row] = atof(token);
//        token = strtok(NULL, ",");
//        T[row] = atof(token);
//        token = strtok(NULL, ",");
//
//        printf("%lf  ", X[row]);
//        printf("%lf  ", Y[row]);
//        printf("%lf  ", L[row]);
//        printf("%lf  ", D[row]);
//        printf("%lf  ", T[row]);
//        printf("\n");
//        row++;
//    }
//    fclose(the_file);
//}

void storeTempOld(int xsize, int ysize, int zsize, double*** temp, double*** tempold)
{
    int i, j, k;

    for (i = 0; i < xsize; i++){
        for (j = 0; j < ysize; j++){
            for (k = 0; k < zsize; k++){
                tempold[i][j][k] = temp[i][j][k];
            }
        }
    }
}

void Update_Temperature_Field(double ***temp, double ***temp_old, int xsize, int ysize, int zsize, double *yc, double timeca, double Gy, double Vy){
    int i, j, k;
    double SF;

    SF = 0.0 + Vy*timeca;
    for (i = 0; i < xsize; i++) {
        for (j = 0; j < ysize; j++) {
            for (k = 0; k < zsize; k++) {

                temp_old[i][j][k] = temp[i][j][k];
                temp[i][j][k] = metalts + (yc[j]-SF)*Gy;

            }
        }
    }

}

double calculateIntegrandVal(double timeVal, double timeca, double xVal, double yVal, double zVal)
{
    double integrandVal;
    double innerConstant; // to avoid repeated calculation, calculate a constant first
    double phiX, phiY, phiZ; // inner constant

    innerConstant = 12 * thermalDiffusivity * (timeca - timeVal);

    phiX = innerConstant + pow(beamRadiusX, 2);
    phiY = innerConstant + pow(beamRadiusY, 2);
    phiZ = innerConstant + pow(beamRadiusZ, 2);

    double laserPosition[2];

//    // calculate heat source position of a spiral track
//    calculateHeatSourcePositionSpiral(timeVal, laserPosition);

    // calculate heat source position of a single track
    calculateHeatSourcePositionSingleTrack(timeVal, laserPosition);

    // integrandVal has a long expression, please check the analytical formula in the report
    integrandVal = 1 / sqrt(phiX * phiY * phiZ) *
                   exp(-3 * pow(xVal - laserPosition[0], 2) / phiX -
                       3 * pow(yVal - laserPosition[1], 2) / phiY -
                       3 * pow(zVal, 2) /phiZ) ;

    return  integrandVal;
}

void calculateHeatSourcePositionSpiral(double timeVal, double* laserPosition)
{
    laserPosition[0] = 0.001  - (spiralRadius - radiusDifference / cycleTime * timeVal)
                                   * cos(2 * PI / cycleTime * timeVal);
    laserPosition[1] = 0.001  - (spiralRadius - radiusDifference / cycleTime * timeVal)
                                   * sin(2 * PI / cycleTime * timeVal); // these two expression should be checked with the output of python
}

void calculateHeatSourcePositionSingleTrack(double timeVal, double* laserPosition)
{
    laserPosition[0] = linearCenterX + singleTrackVelocity * timeVal;
    laserPosition[1] = linearCenterY;
}

double calculatePhysicalConstant()
{
    double physicalConstant;
    physicalConstant = 2 * powerAbsorptivity * powerLaser / (materialDensity * specificHeat *
                                                             pow(PI / 3, 1.5));
    return  physicalConstant;
}

double computeVelocitySpiral(double timeca)
{
    double linearVelocity;
    if(timeca <= cycleTime)
    {
        // add both the components of velocity tangent to radius and perpendicular to radius
        linearVelocity = sqrt(pow((2 * PI / cycleTime, 2) * (spiralRadius - timeca / cycleTime) * radiusDifference
                , 2) + pow(radiusDifference / cycleTime, 2));
    }else{
        linearVelocity = 0.0;
    }
    return linearVelocity;
}


