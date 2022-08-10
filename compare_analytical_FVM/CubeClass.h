
/*
 * Author: zyuhao 
 * Date: 2022-04-08
 * Content: This is the head file of a CubeClass
 * contains the geometric parameters of a cube (position, length, start pos)
 * and info about meshes (relative x, y, z 1-d array) x/y/zCoorList
 * 
 */

//

#ifndef COMPARE_ANALYTICAL_FVM_CUBECLASS_H
#define COMPARE_ANALYTICAL_FVM_CUBECLASS_H

#include <iostream>
#include <cmath>

using namespace std;

class CubeClass {
  public:
    double xStart, yStart, zStart; // start position of the cube
    double xLen, yLen, zLen; // length of cube along x,y,z direction
    double meshSize; // size of mesh
    int xsize, ysize, zsize; // number of mesh along x,y,z direction, calculated

    double *xc, *yc, *zc; // 1d array storing each discrete position of XYZ

    // default ctor
    CubeClass();

    // default dtor
    ~CubeClass();

    // function clear dynamic arrays
    void clear();

    // function that set the size of mesh, return true if the meshSize if valid
    // return false if the meshSize is invalid
    bool setMeshSize(double myMeshSize);

    // function that set the start position of the cube
    // theoratically, any variable will be valid
    void setStartPositionXYZ(
         double myXStart,
         double myYStart,
         double myZStart
         );

    // function that set the length of the cube along x,y,z direction
    // return true if the length is valid
    // return false if the length is invalid
    bool setLenghtXYZ(
         double myXLen,
         double myYLen,
         double myZLen
         );

    // function that calculate the xc, yc, zc according based on all other
    // variables, can be used to ctor or refresh after changing variables
    void calSizeCoorXYZ();

    // friend function
    friend ostream& operator<<(ostream& os, const CubeClass& myCube);
};

// overloading << for CubeClass
ostream& operator<<(ostream& os, const CubeClass& myCube);


#endif //COMPARE_ANALYTICAL_FVM_CUBECLASS_H
