//
// Created by zyuhao on 6/24/22.

// this is the unfinished portion of CA code

// A copy and rewriting of CA_main.c predicting the microstructure at certain time step
// using the thermal history. All the "dirty work" has been encapsulated as
// private functions. I have also seperated the CA code from the MTEX reading code
// and added more comments based on my understanding. I hope it can help you
// understand the code faster and apply it through treating a bunch of private
// functions as black box.

// Due to the scale gap of both 3d space and time resolution
// (macro scale for thermal hist while micro scale for microstructure),
// trilinear interpolation is performed for 3d space and linear interpolation
// is performed for small time step, in order to provide CA model with required
// thermal history.

// Good luck. If you need any help, please email me through: zyuhao@umich.edu
// I will be glad to share my understanding about the CA code.
// But it is worth mentioning that the true reason for my to quit is exactly that
// I can neither fully understand the parallel version of CA code nor
// just build an interface to simply use those
// functionalities defined in the CA_parallel due to lack of comments
// and shallow understanding of 3d CA model in the ca_parallel source code.

// I have stuck in understanding CA_parallel for 10 months. So I truly recommend you
// just use functions in CA_serial and re-written a module (MicrostructureOfCube class)
// that take in molten pool history (ThermalHistOfCubeMoltenPoolTrack object pointer)
// and re-written parallel computing process by yourself. (It is just a suggestion.)

#ifndef THERMAL_MULTI_MPI_MICROSTRUCTUREOFCUBE_H
#define THERMAL_MULTI_MPI_MICROSTRUCTUREOFCUBE_H

#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <malloc.h>
#include <string>
#include <cstdio>
#include <cmath>
#include "ThermalHistOfCubeMoltenPoolTrack.h"
#include "constants.h"

using namespace std;

class MicroStructureOfCube {
  private:
    // function that dynamically allocate 1d array
    double *Allocate_1D_Double(int size, double value);


    // delete 1d dynamic array
    void Free_1D_Double(double *a, int size);


    // dynamically allocate 3d array
    // use openmp to accelerate the value initialization of 3d array
    double ***Allocate_3D_Double(int size1, int size2, int size3, double value);


    // delete 3d array
    void Free_3D_Double(double ***a, int size1, int size2, int size3);


    // Read a block of data from a whole field based upon the parameter of
    // each block
    // gd_name: the name of the file to store corresponding data
    // data: 3d array of with a size of a block
    // id: current process id
    // offset_x,y,z: offset of current
    // integer:
    // x,y,zsize: number of meshes along each axis of the total field
    // gd_x,y,zszie: number of meshes along each axis of each block
    void MpiReadBlock(char *gd_name, double ***data, int id, MPI_Status status,
                      int offset_x, int offset_y, int offset_z, int integer,
                      int xsize, int ysize, int zsize, int gd_xsize, int gd_ysize, int gd_zsize);


    // write a block of data into a file based upon the parameter of each block
    // the logic is similar to MPIReadBlock
    void MpiWriteBlock(char *gd_name, double ***data, int id, MPI_Status status,
                       int offset_x, int offset_y, int offset_z, int integer, int xsize, int ysize, int zsize,
                       int gd_xsize, int gd_ysize, int gd_zsize);


    // For one block, send the connecting surface values of a field
    // to neighbouring blocks. The first and last block along x direction
    // will only send one contacting surface to its neighbouring block.
    void double_transfer(int xsize, int ysize, int zsize, int np, int id, double ***data, int TAGin);


    // Capturing algorithm from Gandin&Rappaz
    // *********************************************************

    // Start from Sample coordinate system, rotate (phi1, Phi, phi2) by Bunge convention, to Crytal coordinate system
    // convert from intrinsic orientation (phi1, Phi, phi2) to extrinsic orientation (phi2, Phi, phi1)
    // Right-hand rotation
    // Reference of this conversion
    // https://en.wikipedia.org/wiki/Davenport_chained_rotations
    void Sample_to_Crystal(double phi1, double Phi, double phi2, double *ptSample, double *ptCrystal);


    // Start from Crystal coordinate system, apply three Euler angles in reverse to arrive at Sample coordinate system
    // convert from intrinsic orientation (-phi2, -Phi, -phi1) to extrinsic orientation (-phi1, -Phi, -phi2)
    // https://en.wikipedia.org/wiki/Davenport_chained_rotations
    void Crystal_to_Sample(double phi1, double Phi, double phi2, double *ptCrystal, double *ptSample);

    int Is_in_Octahedron(double phi1, double Phi, double phi2, double xOct, double yOct, double zOct,
                         double rOct, double xc, double yc, double zc);

    /*
    * reference for this subroutine
    * [1] http://www.songho.ca/math/line/line.html
    * [2] Gaddin&Rappaz-1997
    * [3] Gaddin&Rappaz-1999
    */
    void Calculate_New_Octahedron(double phi1, double Phi, double phi2, double xOct, double yOct, double zOct,
      double rOct, double xc, double yc, double zc, double dx, double *newx, double *newy, double *newz, double *newr);


  public:
    /* variables and functions */

    // ctor
    MicroStructureOfCube(ThermalHistOfCubeMoltenPoolTrack* inThermalHistCube);


    // pointer of object to calculate the temperature field
    // its member function could also track temperature from the position
    // of laser to the contour of a specific temperature.
    // This object should have been initialized by a scanning strategy (MultiTrackLaser object)
    ThermalHistOfCubeMoltenPoolTrack* thermalHistCube;


    // required micro CA grid size



    // required micro time step


    // function that update the required micro time step

    // function that update the required macro time


    // function that, for a block,
    // Firstly interpolates the temperature filed according to the micro
    // time first. Secondly, interpolates the temperature field to the micro
    // CA mesh center.



};


#endif //THERMAL_MULTI_MPI_MICROSTRUCTUREOFCUBE_H
