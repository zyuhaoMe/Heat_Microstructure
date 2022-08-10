/* 
 * Author: zyuhao 
 * Date: 2022-04-08
 * Content: This is the source code file of a CubeClass
 * contains the geometric parameters of a cube (position, length, start pos)
 * and info about meshes (relative x, y, z 1-d array) x/y/zCoorList
 * 
 */

//

#include "CubeClass.h"

// default ctor
CubeClass::CubeClass()
{
  xLen = yLen = 2e-3;
  zLen = 1e-3;
  xStart = yStart = zStart = 0.;
  meshSize = 1e-5;

  // refresh size and coor
  calSizeCoorXYZ();
}

// default dtor
CubeClass::~CubeClass()
{
  clear();
}

// function clear dynamic arrays
void CubeClass::clear()
{
  delete [] xc;
  delete [] yc;
  delete [] zc;
}

// function that set the size of mesh, return true if the meshSize if valid
// return false if the meshSize is invalid
bool CubeClass::setMeshSize(double myMeshSize)
{
  meshSize = myMeshSize;

  // refresh size and coor
  calSizeCoorXYZ();
}

// function that set the start position of the cube
// theoratically, any variable will be valid
void CubeClass::setStartPositionXYZ(
     double myXStart,
     double myYStart,
     double myZStart
     )
{
  xStart = myXStart;
  yStart = myYStart;
  zStart = myZStart;

  // refresh size and coor
  calSizeCoorXYZ();
}


// function that set the length of the cube along x,y,z direction
// return true if the length is valid
// return false if the length is invalid
bool CubeClass::setLenghtXYZ(
     double myXLen,
     double myYLen,
     double myZLen
     )
{
  xLen = myXLen;
  yLen = myYLen;
  zLen = myZLen;

  // refresh size and coor
  calSizeCoorXYZ();
}

// function that calculate the x/y/zSize (number of mesh along each direction)
// and the xc, yc, zc according based on all other
// variables, can be used to ctor or refresh after changing variables
void CubeClass::calSizeCoorXYZ()
{
  xsize = int(ceil(xLen / meshSize));
  ysize = int(ceil(yLen / meshSize));
  zsize = int(ceil(zLen / meshSize));

  // prevent memory leak, delete heap and re-allocate
  clear();

  xc = new double[xsize];
  yc = new double[ysize];
  zc = new double[zsize];
  // allocate each position coordinate
  int i; // loop variable

  for (i = 0; i < xsize; i++)
  {
    xc[i] = xStart + (i +0.5) * meshSize;
  }

  for (i = 0; i < ysize; i++)
  {
    yc[i] = yStart + (i +0.5) * meshSize;
  }

  for (i = 0; i < zsize; i++)
  {
    zc[i] = zStart + (i + 0.5) * meshSize;
  }
}

// overloading << for CubeClass
ostream& operator<<(ostream& os, const CubeClass& myCube)
{
  os << endl << "Cube parameters list: " << endl
       << "  Len (x/y/z):      " << myCube.xLen << " / " << myCube.yLen << " / "
       << myCube.zLen << endl
       << "  Start (x/y/z):    " << myCube.xStart << " / " << myCube.yStart
       << " / " << myCube.zStart << endl
       << "  Mesh size:        " << myCube.meshSize << endl
       << "  Mesh num (x/y/z): " << myCube.xsize << " / " << myCube.ysize
       << " / " << myCube.zsize << endl;

  return os;
}