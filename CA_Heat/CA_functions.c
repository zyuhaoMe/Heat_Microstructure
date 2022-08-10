//
// Created by zyuhao on 2022-01-18.
// This file stores all the needed functions for microstructure prediction
// using CA algorithm

// store all the functions used in CA

#include "CA_heat.h"

double Gaussian_Distribution(double mean, double sigma) {
    double u1, u2, r;
    u1 = rand() / (1. * RAND_MAX);
    u2 = rand() / (1. * RAND_MAX);
    r = sqrt(-2.0 * log(u1)) * cos(2 * PI * u2);
    r = mean + sigma * r;
    return r;
}

void Seeding(double ***tnuc, double ns, double nb, double Ts, double Tb, double Sigs, double Sigb, int xsize, int ysize, int zsize, double dx) {
    int i, j, k, n;
    int no_s, no_b;

    no_s = xsize * zsize * pow(dx, 2) * ns;
    no_b = xsize * ysize * zsize * pow(dx, 3) * nb;

    // seeding at the surface
    for (n = 0; n < no_s; n++) {
        int random_i = xsize - 1, random_k = zsize - 1;
        double random_x = rand() / (1. * RAND_MAX);
        double random_z = rand() / (1. * RAND_MAX);

        random_x = (xsize - 1) * random_x;
        random_z = (zsize - 1) * random_z;

        for (i = 0; i < xsize; i++) {
            if (i < random_x && i + 1 > random_x) {
                random_i = i;
                break;
            }
        }

        for (k = 0; k < zsize; k++) {
            if (k < random_z && k + 1 > random_z) {
                random_k = k;
                break;
            }
        }

        double random_T = metalts - Gaussian_Distribution(Ts, Sigs);
        // take the smallest undercooling (largest nuc temperature)
        if (tnuc[random_i][0][random_k] < random_T)
            tnuc[random_i][0][random_k] = random_T;
    }

    // seeding at the bulk
    for (n = 0; n < no_b; n++) {
        int random_i = xsize - 1, random_j = ysize - 1, random_k = zsize - 1;
        double random_x = rand() / (1. * RAND_MAX);
        double random_y = rand() / (1. * RAND_MAX);
        double random_z = rand() / (1. * RAND_MAX);

        random_x = (xsize - 1) * random_x;
        random_y = (ysize - 1) * random_y;
        random_z = (zsize - 1) * random_z + 1.;

        for (i = 0; i < xsize; i++) {
            if (i < random_x && i + 1 > random_x) {
                random_i = i;
                break;
            }
        }

        for (j = 0; j < ysize; j++) {
            if (j < random_y && j + 1 > random_y) {
                random_j = j;
                break;
            }
        }

        for (k = 1; k < zsize; k++) {
            if (k < random_z && k + 1 > random_z) {
                random_k = k;
                break;
            }
        }

        double random_T = metalts - Gaussian_Distribution(Tb, Sigb);
        if (tnuc[random_i][random_j][random_k] < random_T)
            tnuc[random_i][random_j][random_k] = random_T;
    }

    //     for (k = 0; k < zsize; k++) {
    //         for (j = 0; j < ysize; j++) {
    //             for (i = 0; i < xsize; i++) {
    //                 if (tnuc[i][j][k] != 300) {
    //                     if (k == 0) {
    //                         printf("surface nucleus: %d %d %lf \n", i, j, metalts - tnuc[i][j][k]);
    //                     } else {
    //                         printf("bulk nucleus: %d %d %d %lf \n", i, j, k, metalts - tnuc[i][j][k]);
    //                     }
    //                 }
    //             }
    //         }
    //     }

    printf("# of surface nucleus = %d \n", no_s);
    printf("# of bulk nucleus = %d \n", no_b);
}

void Crystal_to_Sample(double phi1, double Phi, double phi2, double *ptCrystal, double *ptSample)
{
    double sine, cosine, ptnew1[3], ptnew2[3];

    cosine = cos(phi2*PI/180);
    sine = sin(phi2*PI/180);

    ptnew1[0] = ptCrystal[0]*cosine + ptCrystal[1]*sine;
    ptnew1[1] = -ptCrystal[0]*sine + ptCrystal[1]*cosine;
    ptnew1[2] = ptCrystal[2];

    cosine = cos(Phi*PI/180);
    sine = sin(Phi*PI/180);

    ptnew2[0] = ptnew1[0];
    ptnew2[1] = ptnew1[1]*cosine + ptnew1[2]*sine;
    ptnew2[2] = -ptnew1[1]*sine + ptnew1[2]*cosine;

    cosine = cos(phi1*PI/180);
    sine = sin(phi1*PI/180);

    ptSample[0] = ptnew2[0]*cosine + ptnew2[1]*sine;
    ptSample[1] = -ptnew2[0]*sine + ptnew2[1]*cosine;
    ptSample[2] = ptnew2[2];
}

void CrossProduct(double *v1, double *v2, double *v12){
    v12[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v12[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v12[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

void Subtraction(double *v1, double *v2, double *v12){
    v12[0] = v1[0] - v2[0];
    v12[1] = v1[1] - v2[1];
    v12[2] = v1[2] - v2[2];
}

// double *CrossProduct(double *v1, double *v2)
// {
//    double *v12 = Allocate_1D_Double(3, 0);
//    v12[0] = v1[1]*v2[2] - v1[2]*v2[1];
//    v12[1] = v1[2]*v2[0] - v1[0]*v2[2];
//    v12[2] = v1[0]*v2[1] - v1[1]*v2[0];
//    return(v12);
// }

// double *Subtraction(double *v1, double *v2)
// {
//    double *v12 = Allocate_1D_Double(3, 0);
//    v12[0] = v1[0] - v2[0];
//    v12[1] = v1[1] - v2[1];
//    v12[2] = v1[2] - v2[2];
//    return(v12);
// }

double VectorLength(double *v)
{
    double l;
    l = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return(l);
}

double Min(double a, double b)
{
    double m;
    if(a>=b)
    {
        m=b;
    }
    else
    {
        m=a;
    }
    return(m);
}

double Max(double a, double b)
{
    double m;
    if(a>=b)
    {
        m=a;
    }
    else
    {
        m=b;
    }
    return(m);
}

void Sample_to_Crystal(double phi1, double Phi, double phi2, double *ptSample, double *ptCrystal)
{
    double sine, cosine, ptnew1[3], ptnew2[3];

    cosine = cos(phi1*PI/180);
    sine = sin(phi1*PI/180);

    ptnew1[0] = ptSample[0]*cosine - ptSample[1]*sine;
    ptnew1[1] = ptSample[0]*sine + ptSample[1]*cosine;
    ptnew1[2] = ptSample[2];

    cosine = cos(Phi*PI/180);
    sine = sin(Phi*PI/180);

    ptnew2[0] = ptnew1[0];
    ptnew2[1] = ptnew1[1]*cosine - ptnew1[2]*sine;
    ptnew2[2] = ptnew1[1]*sine + ptnew1[2]*cosine;

    cosine = cos(phi2*PI/180);
    sine = sin(phi2*PI/180);

    ptCrystal[0] = ptnew2[0]*cosine - ptnew2[1]*sine;
    ptCrystal[1] = ptnew2[0]*sine + ptnew2[1]*cosine;
    ptCrystal[2] = ptnew2[2];

}

/*
 * reference for this subroutine
 * [1] http://www.songho.ca/math/line/line.html
 * [2] Gaddin&Rappaz-1997
 * [3] Gaddin&Rappaz-1999
 */
void Calculate_New_Octahedron(double phi1, double Phi, double phi2, double xOct, double yOct, double zOct,
                              double rOct, double xc, double yc, double zc, double dx, double *newx, double *newy, double *newz, double *newr)
{
    double sine, cosine;
    double af,bf,cf,df,dist,min_dist;
    double t,dist1,dist2,dist3;
    double L12, L13, Lmu, LIS1, LIS2, LJS1, LJS3;
    int i,j,k;

    double *pC = Allocate_1D_Double(3, 0);
    double *pS = Allocate_1D_Double(3, 0);
    double *pA = Allocate_1D_Double(3, 0);
    double *pS1 = Allocate_1D_Double(3, 0);
    double *pS2 = Allocate_1D_Double(3, 0);
    double *pS3 = Allocate_1D_Double(3, 0);
    double *S1S2 = Allocate_1D_Double(3, 0);
    double *S1S3 = Allocate_1D_Double(3, 0);
    double *AI = Allocate_1D_Double(3, 0);
    double *AJ = Allocate_1D_Double(3, 0);
    double *AS1 = Allocate_1D_Double(3, 0);
    double *nF = Allocate_1D_Double(3, 0);
    double *pI = Allocate_1D_Double(3, 0);
    double *pJ = Allocate_1D_Double(3, 0);
    double *IS1 = Allocate_1D_Double(3, 0);
    double *IS2 = Allocate_1D_Double(3, 0);
    double *JS1 = Allocate_1D_Double(3, 0);
    double *JS3 = Allocate_1D_Double(3, 0);
    double *S_Cmu = Allocate_1D_Double(3, 0);
    double *C_Cmu = Allocate_1D_Double(3, 0);

    double *tempv1 = Allocate_1D_Double(3, 0);
    double *tempv2 = Allocate_1D_Double(3, 0);

    pC[0] = xc - xOct;
    pC[1] = yc - yOct;
    pC[2] = zc - zOct;

    Crystal_to_Sample(phi1, Phi, phi2, pC,pS);

    if(fabs(pS[0])+fabs(pS[1])+fabs(pS[2])>=rOct)
    {
        printf("this cell is not yet captured \n");
    }

    // step1: find F
    df = -rOct;
    if(pS[0] >= 0)
        af = 1.;
    else
        af = -1;

    if(pS[1] >= 0)
        bf = 1.;
    else
        bf = -1;

    if(pS[2] >= 0)
        cf = 1.;
    else
        cf = -1;

    // step2: find A
    t = -(af*pS[0]+bf*pS[1]+cf*pS[2]+df)/(af*af+bf*bf+cf*cf);
    pA[0] = pS[0] + t*af;
    pA[1] = pS[1] + t*bf;
    pA[2] = pS[2] + t*cf;

    // step3: Find S1, S2 and S3
    dist1 = (pA[0]-af*(-df))*(pA[0]-af*(-df))+pA[1]*pA[1]+pA[2]*pA[2];
    dist2 = pA[0]*pA[0]+(pA[1]-bf*(-df))*(pA[1]-bf*(-df))+pA[2]*pA[2];
    dist3 = pA[0]*pA[0]+pA[1]*pA[1]+(pA[2]-cf*(-df))*(pA[2]-cf*(-df));

    if(dist1<=dist2 && dist1<=dist3)
    {
        pS1[0]=af*(-df);pS1[1]=0.;pS1[2]=0.;
        pS2[0]=0.;pS2[1]=bf*(-df);pS2[2]=0.;
        pS3[0]=0.;pS3[1]=0.;pS3[2]=cf*(-df);
    }
    else if(dist2<=dist1 && dist2<=dist3)
    {
        pS1[0]=0.;pS1[1]=bf*(-df);pS1[2]=0.;
        pS2[0]=af*(-df);pS2[1]=0.;pS2[2]=0.;
        pS3[0]=0.;pS3[1]=0.;pS3[2]=cf*(-df);
    }
    else
    {
        pS1[0]=0.;pS1[1]=0.;pS1[2]=cf*(-df);
        pS2[0]=0.;pS2[1]=bf*(-df);pS2[2]=0.;
        pS3[0]=af*(-df);pS3[1]=0.;pS3[2]=0.;
    }

    // step4: find I and J
    Subtraction(pS2,pS1,S1S2);
    Subtraction(pS3,pS1,S1S3);
    Subtraction(pS1,pA,AS1);
    nF[0] = af; nF[1] = bf; nF[2] = cf;

    CrossProduct(S1S2,nF,AI);
    CrossProduct(S1S3,nF,AJ);

    CrossProduct(AI,S1S2,tempv1);
    CrossProduct(AS1,S1S2,tempv2);
    t = tempv2[0] / tempv1[0];
    pI[0] = pA[0]+t*AI[0]; pI[1] = pA[1]+t*AI[1]; pI[2] = pA[2]+t*AI[2];

    CrossProduct(AJ,S1S3,tempv1);
    CrossProduct(AS1,S1S3,tempv2);
    t = tempv2[0] / tempv1[0];
    pJ[0] = pA[0]+t*AJ[0]; pJ[1] = pA[1]+t*AJ[1]; pJ[2] = pA[2]+t*AJ[2];

    // step5: find Lmu
    Subtraction(pS1,pI,IS1);
    Subtraction(pS2,pI,IS2);
    Subtraction(pS1,pJ,JS1);
    Subtraction(pS3,pJ,JS3);
    LIS1 = VectorLength(IS1);
    LIS2 = VectorLength(IS2);
    LJS1 = VectorLength(JS1);
    LJS3 = VectorLength(JS3);

    L12 = 0.5*(Min(LIS1,sqrt(3.)*dx) + Min(LIS2,sqrt(3.)*dx));
    L13 = 0.5*(Min(LJS1,sqrt(3.)*dx) + Min(LJS3,sqrt(3.)*dx));
    Lmu = sqrt(2./3.)*Max(L12,L13);

    // step6: find Cmu
    dist = VectorLength(pS1);
    tempv1[0] = -pS1[0]/dist; tempv1[1] = -pS1[1]/dist; tempv1[2] = -pS1[2]/dist;
    S_Cmu[0] = pS1[0] + sqrt(3.)*Lmu*tempv1[0];
    S_Cmu[1] = pS1[1] + sqrt(3.)*Lmu*tempv1[1];
    S_Cmu[2] = pS1[2] + sqrt(3.)*Lmu*tempv1[2];
    Sample_to_Crystal(phi1,Phi,phi2,S_Cmu,C_Cmu);

    C_Cmu[0] = C_Cmu[0] + xOct;
    C_Cmu[1] = C_Cmu[1] + yOct;
    C_Cmu[2] = C_Cmu[2] + zOct;

    *newx = C_Cmu[0];
    *newy = C_Cmu[1];
    *newz = C_Cmu[2];
    *newr = sqrt(3.)*Lmu;

    if(sqrt(3.)*Lmu > rOct+1e-6)
    {
        printf("new octahedron is too large!! \n");
    }

    free(pC);
    free(pS);
    free(pA);
    free(pS1);
    free(pS2);
    free(pS3);
    free(S1S2);
    free(S1S3);
    free(AI);
    free(AJ);
    free(AS1);
    free(nF);
    free(pI);
    free(pJ);
    free(IS1);
    free(IS2);
    free(JS1);
    free(JS3);
    free(S_Cmu);
    free(C_Cmu);

}

int Is_in_Octahedron(double phi1, double Phi, double phi2, double xOct, double yOct, double zOct,
                     double rOct, double xc, double yc, double zc)
{
    int it;
    double sine, cosine, ptnew1[3], ptnew2[3];
    double *pS = Allocate_1D_Double(3, 0);
    double *pt = Allocate_1D_Double(3, 0);
    it = 0;

    pt[0] = xc - xOct;
    pt[1] = yc - yOct;
    pt[2] = zc - zOct;

    Crystal_to_Sample(phi1, Phi, phi2, pt,pS);

    if(fabs(pS[0])+fabs(pS[1])+fabs(pS[2])<rOct)
        it = 1;
    free(pt);
    free(pS);
    return(it);
}

void Initialization(double ***phi1, double ***Phi, double ***phi2, double ***nucindex, double *xc, double *yc, double *zc,
                    double ***nucx, double ***nucy, double ***nucz, double grain_size, double dx, int xsize, int ysize, int zsize)
{
    int i, j, k, ii, jj, kk;
    int radius = grain_size/dx, ctgrains=0;
    printf("\n radius = %d",radius);
    for (k=0;k<zsize;k++){
        for (j=0;j<ysize;j++){
            for (i=0;i<xsize;i++){

                int thisi = i;
                int thisj = j;
                int thisk = k;

                int radius1 = radius -1 + 2.95*rand()/(1.*RAND_MAX);
                int radius2 = radius -1 + 2.95*rand()/(1.*RAND_MAX);
                int radius3 = radius -1 + 2.95*rand()/(1.*RAND_MAX);

                int deltai=0,deltaj=0,deltak=0;

                deltai = 0.2*radius*rand()/(1.*RAND_MAX); //thisi%radius;
                deltaj = 0.2*radius*rand()/(1.*RAND_MAX); //thisj%radius;
                deltak= 0.2*radius*rand()/(1.*RAND_MAX);  //thisk%radius;

                if ((thisi+deltai)%radius1==0 && (thisj+deltaj)%radius2==0 && (thisk+deltak)%radius3==0) {
                    //printf("\n %d coor. = %d, %d, %d",id,thisi,thisj,thisk);

                    nucindex[i][j][k] = 100.*rand()/(1.*RAND_MAX);
                    nucx[i][j][k] = xc[i];
                    nucy[i][j][k] = yc[j];
                    nucz[i][j][k] = zc[k];

                    // random orientation
                    Phi[i][j][k]=acos(2.*rand()/(1.*RAND_MAX)-1.)*180./PI;
                    phi2[i][j][k] = (360.*rand()/(1.*RAND_MAX));
                    phi1[i][j][k] = (360.*rand()/(1.*RAND_MAX));

                    // Y fiber
                    // 	  phi1[i][j][k] = 90;
                    // 	  phi2[i][j][k] = 90;
                    // 	  Phi[i][j][k]=acos(2.*rand()/(1.*RAND_MAX)-1.)*180./PI;

                    // Z fiber
                    // 	  phi1[i][j][k] = (360.*rand()/(1.*RAND_MAX));
                    // 	  Phi[i][j][k] = 0;
                    // 	  phi2[i][j][k] = 0;

                    // bi-modal
                    // 	  double which = rand()/(1.*RAND_MAX);
                    // 	  if(which > 0.5)
                    // 	  {
                    // 	    // Y-alligned
                    // 	    Phi[i][j][k]= 0.;
                    // 	    phi2[i][j][k] = 90;
                    // 	    phi1[i][j][k] = 90;
                    // 	    nucindex[i][j][k] = 25;
                    // 	  }
                    // 	  else
                    // 	  {
                    // 	    // Z-alligned
                    // 	    Phi[i][j][k]= 15.;
                    // 	    phi2[i][j][k] = 0.;
                    // 	    phi1[i][j][k] = 0.;
                    // 	    nucindex[i][j][k] = 75;
                    // 	  }

                    //printf("\nnew nuclei %d %d %d %d %f %f %f %d %d %d",ctgrains,i,j,k,nucindex[i][j][k],Phi[i][j][k],phi2[i][j][k],radius1,radius2,radius3);
                    ctgrains++;

                    //lac[i][j][k]=nucindex[i][j][k];

                }

            }
        }
    }

    printf("# of grains = %d \n",ctgrains);
    //getchar();
    int keep=1,cting=0;
    while (cting<2*radius) {
        cting++;
        printf("start %d %d \n",cting,radius);
        keep=0;
        for (k=0;k<zsize;k++){
            for (j=0;j<ysize;j++){
                for (i=0;i<xsize;i++){

                    int thisi = i;
                    int thisj = j;
                    int thisk = k;

                    if (nucindex[i][j][k]>0) {

                        for (ii=i-1;ii<=i+1;ii++){
                            for (jj=j-1;jj<=j+1;jj++){
                                for (kk=k-1;kk<=k+1;kk++){
                                    if (ii>=0 && ii<xsize && jj>=0 && jj<ysize && kk>=0 && kk<zsize && nucindex[ii][jj][kk]<=0.1){
                                        double dist = sqrt((xc[ii]-nucx[i][j][k])*(xc[ii]-nucx[i][j][k]) + (yc[jj]-nucy[i][j][k])*(yc[jj]-nucy[i][j][k]) + (zc[kk]-nucz[i][j][k])*(zc[kk]-nucz[i][j][k]));
                                        //printf("\nfind nb  %d %d %d %f",ii,jj,kk,dist/dx);
                                        if (dist<1.0*dx*cting)
                                        {
                                            nucx[ii][jj][kk] = nucx[i][j][k];
                                            nucy[ii][jj][kk] = nucy[i][j][k];
                                            nucz[ii][jj][kk] = nucz[i][j][k];

                                            nucindex[ii][jj][kk]  = nucindex[i][j][k] ;
                                            Phi[ii][jj][kk]      = Phi[i][j][k];
                                            phi2[ii][jj][kk]       = phi2[i][j][k];
                                            phi1[ii][jj][kk]       = phi1[i][j][k];
                                            keep=1;
                                        }
                                    }
                                }
                            }
                        }
                    }

                }
            }
        }

        if (cting>48 && keep==1) {
            printf("\nerror: still not filling all cells with initial grains!!!!!");
            getchar();
        }

    }

}
