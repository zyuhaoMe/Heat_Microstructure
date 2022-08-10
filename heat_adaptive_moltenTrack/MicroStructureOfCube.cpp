//
// Created by zyuhao on 6/24/22.
//

#include "MicroStructureOfCube.h"

// dynamically allocate 1d array
double* MicroStructureOfCube::Allocate_1D_Double(int size, double value)
{
  double *v;
  int i;

  v = (double *) malloc(size * sizeof(double)) ;

  for (i=0;i<size;i++)
  {
    v[i]=value;
  }
  return (v);
}


// delete 1d dynamic array
void MicroStructureOfCube::Free_1D_Double(double *a, int size)
{
  free((void *)a);
}


// dynamically allocate 3d array
// use openmp to accelerate the value initialization of 3d array
double*** MicroStructureOfCube:: Allocate_3D_Double(int size1, int size2, int size3, double value)
{
  double ***v,*vv;
  int i,j,k;
  vv = (double *)malloc(size1*size2*size3*sizeof(double));
  v = (double ***)malloc(size1*sizeof(double**));

  // MC NOTHR
  //    ANNOTATE_SITE_BEGIN(alloc_3Da);
  for (i=0;  i<size1; i++){
    //	ANNOTATE_ITERATION_TASK(i);
    v[i] = (double **)malloc(size2*sizeof(double *));
    for (j=0; j<size2; j++){
      v[i][j] = vv + i*size2*size3 + j*size3;
    }
  }
  //    ANNOTATE_SITE_END();
  // initialize
  //    ANNOTATE_SITE_BEGIN(alloc_3D);
#pragma omp parallel for private(i,j,k) collapse(2)
  for (i=0; i<size1; i++){
    //	ANNOTATE_ITERATION_TASK(i);
    for (j=0;  j<size2; j++){
      for (k=0; k<size3; k++){
        v[i][j][k] = value;
      }
    }
  }
  //    ANNOTATE_SITE_END();
  return (v);
}


// delete 3d array
void MicroStructureOfCube::Free_3D_Double(double ***a, int size1, int size2, int size3)
{
  int i,j,k;
  double * vv;

  vv = a[0][0];
  for(i=0; i<size1; i++){
    free((void*)a[i]);
  }
  free(a);
  free(vv);
}


// Read a block of data from a whole field
void MicroStructureOfCube::MpiReadBlock(char *gd_name, double ***data, int id, MPI_Status status,
                  int offset_x, int offset_y, int offset_z, int integer,
                  int xsize, int ysize, int zsize, int gd_xsize, int gd_ysize, int gd_zsize)
{
  double *mpi_data; // 3d double array to 1d
  MPI_File fhr; // a file handle in MPI using MPI IO
  int offset_domain; // offset of current block, along x, transfer to 1d
  int i,j,k,l; // loop variable

  // dynamically allocate 1d array (3d array size)
  mpi_data = (double *) malloc((xsize)*(ysize)*(zsize)*sizeof(double));
  // block is divided along x direction
  offset_domain = id*(integer)*gd_ysize*gd_zsize*sizeof(double);

  // open the file in read only mode using mpi
  MPI_File_open(MPI_COMM_WORLD, gd_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhr);
  // set pointer to absolute offset
  MPI_File_seek(fhr, offset_x*gd_ysize*gd_zsize*sizeof(double), MPI_SEEK_SET);
  // set pointer to current position + offset
  MPI_File_seek(fhr, offset_domain, MPI_SEEK_CUR);

  l=0;
  for(i=0;i<xsize;i++)
  {
    MPI_File_seek(fhr, offset_y*gd_zsize*sizeof(double), MPI_SEEK_CUR);
    for(j=0;j<ysize;j++)
    {
      MPI_File_seek(fhr, offset_z*sizeof(double), MPI_SEEK_CUR);
      MPI_File_read(fhr, &mpi_data[l], zsize, MPI_DOUBLE, &status);
      l=l+zsize;
      MPI_File_seek(fhr, (gd_zsize-offset_z-zsize)*sizeof(double), MPI_SEEK_CUR);
    }
    MPI_File_seek(fhr, (gd_ysize-offset_y-ysize)*gd_zsize*sizeof(double), MPI_SEEK_CUR);
  }

  MPI_File_close(&fhr);

  l=0;
  // reconstruct the block from 1d to 3d
  for (i=0;i<xsize;i++){
    for (j=0;j<ysize;j++){
      for (k=0;k<zsize;k++){
        data[i][j][k] = mpi_data[l];
        l = l+1;
      }}}

  free(mpi_data);
}


// write a block of data into a file
// the logic is similar to MPIReadBlock
void MicroStructureOfCube::MpiWriteBlock(char *gd_name, double ***data, int id, MPI_Status status,
                   int offset_x, int offset_y, int offset_z, int integer, int xsize, int ysize, int zsize,
                   int gd_xsize, int gd_ysize, int gd_zsize)
{
  double *mpi_data;
  MPI_File fhw;
  int offset_domain;
  int i,j,k,l;

  mpi_data = (double *) malloc((xsize)*(ysize)*(zsize)*sizeof(double));
  l=0;
  for (i=0;i<xsize;i++){
    for (j=0;j<ysize;j++){
      for (k=0;k<zsize;k++){
        mpi_data[l] = data[i][j][k];
        l = l+1;
      }}}

  offset_domain = id*(integer)*gd_ysize*gd_zsize*sizeof(double);

  MPI_File_open(MPI_COMM_WORLD, gd_name, MPI_MODE_WRONLY, MPI_INFO_NULL, &fhw);
  MPI_File_seek(fhw, offset_x*gd_ysize*gd_zsize*sizeof(double), MPI_SEEK_SET);
  MPI_File_seek(fhw, offset_domain, MPI_SEEK_CUR);

  l=0;
  for(i=0;i<xsize;i++)
  {
    if(i%10 == 0 && id ==0) printf("\n writing %d", i);
    MPI_File_seek(fhw, offset_y*gd_zsize*sizeof(double), MPI_SEEK_CUR);
    for(j=0;j<ysize;j++)
    {
      MPI_File_seek(fhw, offset_z*sizeof(double), MPI_SEEK_CUR);
      MPI_File_write(fhw, &mpi_data[l], zsize, MPI_DOUBLE, &status);
      l=l+zsize;
      MPI_File_seek(fhw, (gd_zsize-offset_z-zsize)*sizeof(double), MPI_SEEK_CUR);
    }
    MPI_File_seek(fhw, (gd_ysize-offset_y-ysize)*gd_zsize*sizeof(double), MPI_SEEK_CUR);
  }

  MPI_File_close(&fhw);
  free(mpi_data);

}



// For one block, send the connecting surface values of a field
// to neighbouring blocks. The first and last block along x direction
// will only send one contacting surface to its neighbouring block.
void double_transfer(int xsize, int ysize, int zsize, int np, int id, double ***data, int TAGin)
{
  MPI_Status status;
  MPI_Request req[2*ysize];

  int j,k;
  int LeftNB, RightNB;


  if (zsize>5000) {
    printf("\narray size is not big energy");
  }

  LeftNB  = id-1;
  RightNB = id+1;

  int Delta = 0;
  if(np != 1){
    for (j=0;j<ysize;j++)
    {
      int TAG = TAGin*1000+j;
      if (LeftNB >= 0){

        MPI_Send(&data[1][j][0], zsize, MPI_DOUBLE, LeftNB, TAG, MPI_COMM_WORLD);
        MPI_Irecv(&data[0][j][0], zsize, MPI_DOUBLE, LeftNB, TAG, MPI_COMM_WORLD, &req[j]);

      }

      if (RightNB < np){

        MPI_Send(&data[xsize-2][j][0], zsize, MPI_DOUBLE, RightNB, TAG, MPI_COMM_WORLD);
        MPI_Irecv(&data[xsize-1][j][0], zsize, MPI_DOUBLE, RightNB, TAG, MPI_COMM_WORLD, &req[ysize+j]);

      }


    }

    if(id == 0)
    {
      for (j=ysize;j<2*ysize;j++)
        MPI_Wait(&req[j],&status);
    }
    else if(id == np - 1)
    {
      for (j=0;j<ysize;j++)
        MPI_Wait(&req[j],&status);
    }
    else
    {
      for (j=0;j<2*ysize;j++)
        MPI_Wait(&req[j],&status);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

}



// For one block, send the connecting surface values of a field
// to neighbouring blocks. The first and last block along x direction
// will only send one contacting surface to its neighbouring block.
void MicroStructureOfCube::double_transfer(int xsize, int ysize, int zsize, int np, int id, double ***data, int TAGin) {
  MPI_Status status;
  MPI_Request req[2 * ysize];

  int j, k;
  int LeftNB, RightNB;


  if (zsize > 5000) {
    printf("\narray size is not big energy");
  }

  LeftNB = id - 1;
  RightNB = id + 1;

  int Delta = 0;
  if (np != 1) {
    for (j = 0; j < ysize; j++) {
      int TAG = TAGin * 1000 + j;
      if (LeftNB >= 0) {

        MPI_Send(&data[1][j][0], zsize, MPI_DOUBLE, LeftNB, TAG,
                 MPI_COMM_WORLD);
        MPI_Irecv(&data[0][j][0], zsize, MPI_DOUBLE, LeftNB, TAG,
                  MPI_COMM_WORLD, &req[j]);

      }

      if (RightNB < np) {

        MPI_Send(&data[xsize - 2][j][0], zsize, MPI_DOUBLE, RightNB, TAG,
                 MPI_COMM_WORLD);
        MPI_Irecv(&data[xsize - 1][j][0], zsize, MPI_DOUBLE, RightNB, TAG,
                  MPI_COMM_WORLD, &req[ysize + j]);

      }


    }

    if (id == 0) {
      for (j = ysize; j < 2 * ysize; j++)
        MPI_Wait(&req[j], &status);
    } else if (id == np - 1) {
      for (j = 0; j < ysize; j++)
        MPI_Wait(&req[j], &status);
    } else {
      for (j = 0; j < 2 * ysize; j++)
        MPI_Wait(&req[j], &status);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}



// Capturing algorithm from Gandin&Rappaz
// *********************************************************

// Start from Sample coordinate system, rotate (phi1, Phi, phi2) by Bunge convention, to Crytal coordinate system
// convert from intrinsic orientation (phi1, Phi, phi2) to extrinsic orientation (phi2, Phi, phi1)
// Right-hand rotation
// Reference of this conversion
// https://en.wikipedia.org/wiki/Davenport_chained_rotations
void MicroStructureOfCube::Sample_to_Crystal(double phi1, double Phi, double phi2, double *ptSample, double *ptCrystal)
{
  double sine, cosine, ptnew1[3], ptnew2[3];

  cosine = cos(phi2*PI/180);
  sine = sin(phi2*PI/180);

  ptnew1[0] = ptSample[0]*cosine - ptSample[1]*sine;
  ptnew1[1] = ptSample[0]*sine + ptSample[1]*cosine;
  ptnew1[2] = ptSample[2];

  cosine = cos(Phi*PI/180);
  sine = sin(Phi*PI/180);

  ptnew2[0] = ptnew1[0];
  ptnew2[1] = ptnew1[1]*cosine - ptnew1[2]*sine;
  ptnew2[2] = ptnew1[1]*sine + ptnew1[2]*cosine;

  cosine = cos(phi1*PI/180);
  sine = sin(phi1*PI/180);

  ptCrystal[0] = ptnew2[0]*cosine - ptnew2[1]*sine;
  ptCrystal[1] = ptnew2[0]*sine + ptnew2[1]*cosine;
  ptCrystal[2] = ptnew2[2];
}



// Start from Crystal coordinate system, apply three Euler angles in reverse to arrive at Sample coordinate system
// convert from intrinsic orientation (-phi2, -Phi, -phi1) to extrinsic orientation (-phi1, -Phi, -phi2)
// https://en.wikipedia.org/wiki/Davenport_chained_rotations
void MicroStructureOfCube::Crystal_to_Sample(double phi1, double Phi, double phi2, double *ptCrystal, double *ptSample)
{
  double sine, cosine, ptnew1[3], ptnew2[3];

  cosine = cos(phi1*PI/180);
  sine = sin(phi1*PI/180);

  ptnew1[0] = ptCrystal[0]*cosine + ptCrystal[1]*sine;
  ptnew1[1] = -ptCrystal[0]*sine + ptCrystal[1]*cosine;
  ptnew1[2] = ptCrystal[2];

  cosine = cos(Phi*PI/180);
  sine = sin(Phi*PI/180);

  ptnew2[0] = ptnew1[0];
  ptnew2[1] = ptnew1[1]*cosine + ptnew1[2]*sine;
  ptnew2[2] = -ptnew1[1]*sine + ptnew1[2]*cosine;

  cosine = cos(phi2*PI/180);
  sine = sin(phi2*PI/180);

  ptSample[0] = ptnew2[0]*cosine + ptnew2[1]*sine;
  ptSample[1] = -ptnew2[0]*sine + ptnew2[1]*cosine;
  ptSample[2] = ptnew2[2];

}

void CrossProduct(double *v1, double *v2, double *v12)
{
  v12[0] = v1[1]*v2[2] - v1[2]*v2[1];
  v12[1] = v1[2]*v2[0] - v1[0]*v2[2];
  v12[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

void Subtraction(double *v1, double *v2, double *v12)
{
  v12[0] = v1[0] - v2[0];
  v12[1] = v1[1] - v2[1];
  v12[2] = v1[2] - v2[2];
}

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

int MicroStructureOfCube::Is_in_Octahedron(double phi1, double Phi, double phi2, double xOct, double yOct, double zOct, \
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

  Crystal_to_Sample(phi1, Phi, phi2, pt, pS);

  if(fabs(pS[0])+fabs(pS[1])+fabs(pS[2])<rOct)
    it = 1;

  free(pS); free(pt);
  return(it);
}

/*
 * reference for this subroutine
 * [1] http://www.songho.ca/math/line/line.html
 * [2] Gaddin&Rappaz-1997
 * [3] Gaddin&Rappaz-1999
 */
void MicroStructureOfCube::Calculate_New_Octahedron(double phi1, double Phi, double phi2, double xOct, double yOct, double zOct, \
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

  Crystal_to_Sample(phi1, Phi, phi2, pC, pS);

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
  Subtraction(pS2,pS1, S1S2);
  Subtraction(pS3,pS1, S1S3);
  Subtraction(pS1,pA, AS1);
  nF[0] = af; nF[1] = bf; nF[2] = cf;

  CrossProduct(S1S2,nF, AI);
  CrossProduct(S1S3,nF, AJ);

  CrossProduct(AI,S1S2, tempv1);
  CrossProduct(AS1,S1S2, tempv2);
  t = tempv2[0] / tempv1[0];
  pI[0] = pA[0]+t*AI[0]; pI[1] = pA[1]+t*AI[1]; pI[2] = pA[2]+t*AI[2];

  CrossProduct(AJ,S1S3, tempv1);
  CrossProduct(AS1,S1S3, tempv2);
  t = tempv2[0] / tempv1[0];
  pJ[0] = pA[0]+t*AJ[0]; pJ[1] = pA[1]+t*AJ[1]; pJ[2] = pA[2]+t*AJ[2];

  // step5: find Lmu
  Subtraction(pS1,pI, IS1);
  Subtraction(pS2,pI, IS2);
  Subtraction(pS1,pJ, JS1);
  Subtraction(pS3,pJ, JS3);
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
  Sample_to_Crystal(phi1,Phi,phi2,S_Cmu, C_Cmu);

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
  free(tempv1);
  free(tempv2);
}


