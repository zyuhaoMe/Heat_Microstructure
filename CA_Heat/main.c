// the main function based on CA code
// the source code of heat transfer will be stored in heat.h
#include "CA_heat.h"

int main(int argc, char *argv[]) {

    int i, j, k, ii, jj, kk;
    int xsize, ysize, zsize;
    double startingx, startingy, startingz;
    int timeend;
    double dtca, dx;
    double *pS, *pS2;
    int isin_old, isin_new;

    double Gy, Vy, Tdot;
    double nmax_s, nmax_b, TN_s, TN_b, Tsig_s, Tsig_b;
    int count_max_vel, count_growing_cell, stop_grow;
    double dist, SF, grain_size;
    double timeTotal = 0.1 * 1.8 * singleTrackTime; // total range of time of interest
    char filename[100]= ""; // initialize filename
    int phase_index;
    FILE* fw;

    double length, width, height; // size of area of interest
    length = 1.5e-3;
    width = 5e-4;
    height = 1e-4;
    dx = 3e-6; // mesh size

    xsize = (int)(length / dx);
    ysize = (int)(width / dx);
    zsize = (int)(height / dx);

    startingx = 0.;
    startingy = 0.;
    startingz = 0.;

    Gy = 1.5e5;

    // the velocity of solidification front equals to the laser speed
//    Vy = 2 * PI * spiralRadius / cycleTime; // for spiral path
    Vy = singleTrackVelocity; // for single track path

    dtca = 0.1 * dx / Vy; // time step size
    // dtca = 0.1 * dx / 1.0;  // approximate interface velocity is 1.0 m/s
    timeend = (int)(timeTotal / dtca);

    // surface nucleation condition at the bottom surface
    nmax_s = 1.5e9;
    TN_s = 0.1;
    Tsig_s = 0.01;

    // bulk nucleation condition
    nmax_b = 1e10;
    TN_b = 0.1;
    Tsig_b = 0.01;
    grain_size = 30e-6;  // initial grain size if want to have solid regions (not used here)

    Tdot = Gy * Vy;
    //srand (time(NULL));

    double *xc = Allocate_1D_Double(xsize, 0);
    double *yc = Allocate_1D_Double(ysize, 0);
    double *zc = Allocate_1D_Double(zsize, 0);

    for (i = 0; i < xsize; i++)
        xc[i] = (1. * i) * dx + startingx;
    for (j = 0; j < ysize; j++)
        yc[j] = (1. * j) * dx + startingy;
    for (k = 0; k < zsize; k++)
        zc[k] = (1. * k) * dx + startingz;
    double ***temp = Allocate_3D_Double(xsize, ysize, zsize, 300);
    double ***tempold = Allocate_3D_Double(xsize, ysize, zsize, 300);
    double ***tempnuc = Allocate_3D_Double(xsize, ysize, zsize, 300);
    double ***phase = Allocate_3D_Double(xsize, ysize, zsize, 0);
    double ***nucindex = Allocate_3D_Double(xsize, ysize, zsize, 0);
    double ***nucx = Allocate_3D_Double(xsize, ysize, zsize, 0);
    double ***nucy = Allocate_3D_Double(xsize, ysize, zsize, 0);
    double ***nucz = Allocate_3D_Double(xsize, ysize, zsize, 0);
    double ***Phi = Allocate_3D_Double(xsize, ysize, zsize, 0);
    double ***phi2 = Allocate_3D_Double(xsize, ysize, zsize, 0);
    double ***phi1 = Allocate_3D_Double(xsize, ysize, zsize, 0);
    double ***nucr = Allocate_3D_Double(xsize, ysize, zsize, 0);
    double ***nucrold = Allocate_3D_Double(xsize, ysize, zsize, 0);
    double ***nucv = Allocate_3D_Double(xsize, ysize, zsize, 0);

    //**********************************
    /*
     * phase = 0 >> solid
     * phase = 1 >> liquid
     * phase = 2 >> growing cells
     */
    //***********************************

//    //the initial conditions
//    //temperature gradient in Y direction
//    SF = yc[0] - dx;
//    for (i = 0; i < xsize; i++) {
//        for (j = 0; j < ysize; j++) {
//            for (k = 0; k < zsize; k++) {
//                temp[i][j][k] = metalts + (yc[j]-SF)*Gy;
//                tempold[i][j][k] = temp[i][j][k];
//            }
//        }
//    }

    //****************************
    // seeding
    Seeding(tempnuc, nmax_s, nmax_b, TN_s, TN_b, Tsig_s, Tsig_b, xsize, ysize, zsize, dx);
    tempnuc[xsize/2][0][zsize/2] = metalts - 1e-6;

    Initialization(phi1, Phi, phi2, nucindex, xc, yc, zc, nucx, nucy, nucz, grain_size, dx, xsize, ysize, zsize);
    //     tempnuc[12][0][12] = metalts - 1e-6;
    //     tempnuc[37][0][37] = metalts - 1e-6;
    //time marching
    int timect = 0;
    double timeca = 0;

    double phyConst = calculatePhysicalConstant();// calculate the physical constant for temperature calculation
    double laserTime = singleTrackTime; // time range when the heat source is on

    //  time advancement
    for (timect = 0; timect < timeend; timect++) {

        timeca = timeca + dtca;
//        SF = SF + Vy*dtca;

        //substep1: update the temperature and phase field, update the nucindex field

        // store the temperature field from temp into the tempold
        storeTempOld(xsize, ysize, zsize, temp, tempold);

        // Analytical model to update temperature field
        // this is the ONLY function to calculate thermal history
        computeThermalFieldWithGaussianLegendre(
                16,
                xsize,
                ysize,
                zsize,
                xc,
                yc,
                zc,
                dtca,
                timeca,
                temp,
                phyConst,
                laserTime
                );
        printf("Finish calculating temperature distribution of time step %d \n \n", timect);


        for (i=0;i<xsize;i++){
            for (j=0;j<ysize;j++){
                for (k=0;k<zsize;k++){
                    if (phase[i][j][k]<=0 && temp[i][j][k]>metalts)
                    {
                        phase[i][j][k] = 1;
                        nucindex[i][j][k] =0;
                        nucx[i][j][k]=0;
                        nucy[i][j][k]=0;
                        nucz[i][j][k]=0;
                        phi1[i][j][k]=0;
                        Phi[i][j][k]=0;
                        phi2[i][j][k]=0;
                        nucr[i][j][k]=0;
                        nucrold[i][j][k]=0;
                        nucv[i][j][k]=0;
                    }
                }
            }
        }
        //////////////////
        for (i=0;i<xsize;i++){
            for (j=0;j<ysize;j++){
                for (k=0;k<zsize;k++){

                    if (temp[i][j][k]<=metalts && tempold[i][j][k]>metalts && phase[i][j][k] == 1){
                        int flag=0;
                        for (ii=i-1;ii<=i+1;ii++){
                            for (jj=j-1;jj<=j+1;jj++){
                                for (kk=k-1;kk<=k+1;kk++){

                                    double dist=fabs(1.*ii-1.*i)+fabs(1.*jj-1.*j)+fabs(1.*kk-1.*k);
                                    if (dist>0.1 && dist<1.1 && ii>=0 && ii<xsize && jj>=0 && jj<ysize && kk>=0 && kk<zsize){
                                        if (phase[ii][jj][kk]==0 && fabs(nucindex[ii][jj][kk])>1e-9) {
                                            flag=1;
                                            phi1[i][j][k]=phi1[ii][jj][kk];
                                            Phi[i][j][k]=Phi[ii][jj][kk];
                                            phi2[i][j][k]=phi2[ii][jj][kk];
                                            nucindex[i][j][k]=nucindex[ii][jj][kk];
                                        }
                                    }
                                }}}

                        if (flag>0){

                            phase[i][j][k]=2; //phase=2 is nuclei
                            nucrold[i][j][k]=0;
                            nucr[i][j][k]=1e-3*dx;
                            nucx[i][j][k]=xc[i];
                            nucy[i][j][k]=yc[j];
                            nucz[i][j][k]=zc[k];

                        }
                    }

                }}}

        //substep2: nucleation

        for (i = 0; i < xsize; i++) {
            for (j = 0; j < ysize; j++) {
                for (k = 0; k < zsize; k++) {
                    if (tempold[i][j][k] > tempnuc[i][j][k] && temp[i][j][k] < tempnuc[i][j][k] && phase[i][j][k] == 1) {

                        phase[i][j][k] = 2;
                        nucrold[i][j][k] = 0;
                        nucr[i][j][k] = 1e-3 * dx;
                        nucx[i][j][k] = xc[i];
                        nucy[i][j][k] = yc[j];
                        nucz[i][j][k] = zc[k];
                        nucindex[i][j][k] = 100. * rand() / (1. * RAND_MAX);

                        // random orientation
                        // printf("nucleation \n");
                        Phi[i][j][k]=acos(2.*rand()/(1.*RAND_MAX)-1.)*180./PI;
                        phi2[i][j][k] = (360.*rand()/(1.*RAND_MAX));
                        phi1[i][j][k] = (360.*rand()/(1.*RAND_MAX));

                        // 	      Phi[i][j][k]= 20.;
                        // 	      phi2[i][j][k] = 20.;
                        // 	      phi1[i][j][k] = 20.;

                        // Y-fiber
                        // 	      phi1[i][j][k] = 90;
                        // 	      phi2[i][j][k] = 90;
                        // 	      Phi[i][j][k]=acos(2.*rand()/(1.*RAND_MAX)-1.)*180./PI;

                        /*if(i==12 && j==0 && k==12)
                          {
                          nucindex[i][j][k] = 25.;
                          Phi[i][j][k]= 110.;
                          phi2[i][j][k] = 92.;
                          phi1[i][j][k] = 173.;
                          }
                          else
                          {
                          nucindex[i][j][k] = 75.;
                          Phi[i][j][k]= 316.;
                          phi2[i][j][k] = 10.;
                          phi1[i][j][k] = 208.;
                          }	*/

                        if(i==xsize/2 && j==0 && k==zsize/2)
                        {
                            printf("Y-fiber seeding \n");
                            Phi[i][j][k] = 45.;
                            phi2[i][j][k] = 90.;
                            phi1[i][j][k] = 90.;
                            nucindex[i][j][k] = 50;
                        }

                    }
                }
            }
        }

        count_growing_cell = 0;
        count_max_vel = 0;
        //substep3: calculate growth of existing interface cell
        for (i = 0; i < xsize; i++) {
            for (j = 0; j < ysize; j++) {
                for (k = 0; k < zsize; k++) {

                    // octahedron growing stop condition: nucr < 3*dx
                    if (phase[i][j][k] > 1.1 && nucr[i][j][k] > 0)
                    {
                        // octahedron growing stop condition: all neighbor cells are mushy cells
                        stop_grow = 1;
                        for (ii = i - 1; ii <= i + 1; ii++) {
                            for (jj = j - 1; jj <= j + 1; jj++) {
                                for (kk = k - 1; kk <= k + 1; kk++) {
                                    dist = fabs(1. * ii - 1. * i) + fabs(1. * jj - 1. * j) + fabs(1. * kk - 1. * k);
                                    if (dist > 0.1 && ii >= 0 && ii < xsize && jj >= 0 && jj < ysize && kk >= 0 && kk < zsize)
                                    {
                                        if(nucr[ii][jj][kk] <= 0 && phase[ii][jj][kk] > 0.1)
                                        {
                                            stop_grow = 0;
                                            goto endloop1;
                                        }
                                    }
                                }
                            }
                        }

                        endloop1:

                        if(stop_grow != 1)
                        {
                            if(nucr[i][j][k] / dx > 10)
                            {
                                printf("warning: large octahedron still growing \n");
                            }

                            double undercooling = metalts - temp[i][j][k];

                            double vel = 1.0909E-05 * pow(undercooling, 3)
                                         - 2.0336E-04 * pow(undercooling, 2)
                                         + 2.7397E-03 * undercooling
                                         + 1.1504E-04;

                            //double vel = 2.976247708e-5 * pow(undercooling, 2) + 5.171907896e-5 * pow(undercooling, 3);

                            if (vel < 0) vel = 0;
                            nucv[i][j][k] = vel;
                            double dl = nucv[i][j][k] * dtca;

                            if (dl > 0.6 * dx) {
                                //printf("\n Growing too fast %f %f %f %f %f", tempold[i][j][k], temp[i][j][k], nucrold[i][j][k] / dx, nucr[i][j][k] / dx, vel * dtca / dx);
                                count_max_vel = count_max_vel + 1;
                                dl = 0.6 * dx;
                            }

                            if(vel > 0)
                            {
                                count_growing_cell = count_growing_cell + 1;
                            }

                            nucrold[i][j][k] = nucr[i][j][k];
                            nucr[i][j][k] = nucrold[i][j][k] + dl;
                        }
                        else
                        {
                            phase[i][j][k] = 0;
                        }
                    }

                }
            }
        }


        if (count_max_vel > 20) {
            printf("************************************ \n");
            printf("warning: # of cells growing too fast = %d # of growing cells = %d \n", count_max_vel, count_growing_cell);
            printf("************************************ \n");
        }

        printf("time step =  %d # of growing cells = %d \n", timect, count_growing_cell);
        //substep4: find new interface cell
        for (i = 0; i < xsize; i++) {
            for (j = 0; j < ysize; j++) {
                for (k = 0; k < zsize; k++) {
                    // huntee cell: (i,j,k)
                    // hunter cell: (ii,jj,kk)
                    double min = 1e6;

                    if (phase[i][j][k] > 0.1 && nucr[i][j][k] <= 1e-12) {
                        for (ii = i - 1; ii <= i + 1; ii++) {
                            for (jj = j - 1; jj <= j + 1; jj++) {
                                for (kk = k - 1; kk <= k + 1; kk++) {

                                    double dist = fabs(1. * ii - 1. * i) + fabs(1. * jj - 1. * j) + fabs(1. * kk - 1. * k);
                                    if (dist > 0.1 && ii >= 0 && ii < xsize && jj >= 0 && jj < ysize && kk >= 0 \
                                 && kk < zsize && nucrold[ii][jj][kk] > 0){

                                        isin_old = Is_in_Octahedron(phi1[ii][jj][kk], Phi[ii][jj][kk], phi2[ii][jj][kk], \
                                    nucx[ii][jj][kk], nucy[ii][jj][kk], nucz[ii][jj][kk], nucrold[ii][jj][kk], \
                                    xc[i], yc[j], zc[k]);

                                        isin_new = Is_in_Octahedron(phi1[ii][jj][kk], Phi[ii][jj][kk], phi2[ii][jj][kk], \
                                    nucx[ii][jj][kk], nucy[ii][jj][kk], nucz[ii][jj][kk], nucr[ii][jj][kk], \
                                    xc[i], yc[j], zc[k]);

                                        if(isin_old <= 0 && isin_new > 0)
                                        {
                                            double win;

                                            // competition criterion: distance to (nucx, nucy, nucz)
                                            win = sqrt(pow(xc[i] - nucx[ii][jj][kk], 2) + pow(yc[j] - nucy[ii][jj][kk], 2) \
                                       + pow(zc[k] - nucz[ii][jj][kk], 2));

                                            if(min > win)
                                            {
                                                min = win;

                                                Calculate_New_Octahedron(phi1[ii][jj][kk], Phi[ii][jj][kk], phi2[ii][jj][kk], \
                                          nucx[ii][jj][kk], nucy[ii][jj][kk], nucz[ii][jj][kk], nucr[ii][jj][kk], \
                                          xc[i], yc[j], zc[k], dx, &nucx[i][j][k], &nucy[i][j][k], &nucz[i][j][k], \
                                          &nucr[i][j][k]);

                                                nucrold[i][j][k] = 0;
                                                phase[i][j][k] = 2;
                                                nucindex[i][j][k] = nucindex[ii][jj][kk];
                                                phi1[i][j][k] = phi1[ii][jj][kk];
                                                Phi[i][j][k] = Phi[ii][jj][kk];
                                                phi2[i][j][k] = phi2[ii][jj][kk];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                }
            }
        }

//        if(timect%100 == 0 || timect == timeend-1)
//        {
//            // output to tecplot
//            sprintf(filename, "./tecplot_output_2/tecplot_%d.dat", timect);
//            printf("writing output \n");
//            fw = fopen(filename, "w");
//
//            fprintf(fw,"VARIABLES=\"X\", \"Y\", \"Z\", \"nucindex\", \"temperature\" \n");
//            fprintf(fw,"ZONE T=\"%d\",I=%i, J=%i, K=%i, DATAPACKING=POINT \n", timect, xsize, ysize, zsize);
//            fprintf(fw,"SOLUTIONTIME=%le \n", timect * dtca);
//            for(k = 0; k < zsize; k++){
//                for(j = 0; j < ysize; j++){
//                    for(i = 0; i < xsize; i++){
//                        fprintf(fw, "%le %le %le %le %le \n", xc[i], yc[j], zc[k], nucindex[i][j][k], temp[i][j][k]);
//                    }
//                }
//            }
//            fclose(fw);
//        }
    }

//    // write out data for mtex analysis
//    sprintf(filename, "./mtex/Euler.dat");
//    printf("writing Euler angles \n");
//    fw = fopen(filename, "w");
//    fprintf(fw,"VARIABLES=\"x(mm)\" \"y(mm)\" \"z(mm)\" \"nucindex\" \"phi1\" \"Phi\" \"phi2\" \n");
//    fprintf(fw,"ZONE T=\"%d\",I=%i, J=%i, K=%i, F=POINT \n", timeend, xsize, ysize, zsize);
//    for(k = 0; k < zsize; k++){
//        for(j = 0; j < ysize; j++){
//            for(i = 0; i < xsize; i++){
//                fprintf(fw, "%le %le %le %le %le %le %le \n", xc[i], yc[j], zc[k], nucindex[i][j][k], phi1[i][j][k], \
//                  Phi[i][j][k], phi2[i][j][k]);
//            }
//        }
//    }
    fclose(fw);

    Free_1D_Double(xc, xsize);
    Free_1D_Double(yc, ysize);
    Free_1D_Double(zc, zsize);
    Free_3D_Double(temp, xsize, ysize, zsize);
    Free_3D_Double(tempold, xsize, ysize, zsize);
    Free_3D_Double(tempnuc, xsize, ysize, zsize);
    Free_3D_Double(phase, xsize, ysize, zsize);
    Free_3D_Double(nucindex, xsize, ysize, zsize);
    Free_3D_Double(nucx, xsize, ysize, zsize);
    Free_3D_Double(nucy, xsize, ysize, zsize);
    Free_3D_Double(nucz, xsize, ysize, zsize);
    Free_3D_Double(Phi, xsize, ysize, zsize);
    Free_3D_Double(phi2, xsize, ysize, zsize);
    Free_3D_Double(phi1, xsize, ysize, zsize);
    Free_3D_Double(nucr, xsize, ysize, zsize);
    Free_3D_Double(nucrold, xsize, ysize, zsize);
    Free_3D_Double(nucv, xsize, ysize, zsize);

    printf("end of program\n");
}


