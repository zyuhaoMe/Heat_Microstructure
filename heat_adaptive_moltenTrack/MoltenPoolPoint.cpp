//
// Created by zyuhao on 6/9/22.
// implementation file of MoltenPoolPoint class
//

#include "MoltenPoolPoint.h"


MoltenPoolPoint::MoltenPoolPoint(
     int inIndexX,
     int inIndexY,
     int inIndexZ
     )
{
  indexX = inIndexX;
  indexY = inIndexY;
  indexZ = inIndexZ;
}

bool MoltenPoolPoint::operator==(const MoltenPoolPoint& rhs) const
{
  if (
       indexX == rhs.indexX &&
       indexY == rhs.indexY &&
       indexZ == rhs.indexZ
      )
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool MoltenPoolPoint::operator < (const MoltenPoolPoint& rhs) const
{
  if (
       (indexX * xsize * ysize + indexY * ysize + indexZ) <
  (rhs.indexX * xsize * ysize + rhs.indexY * ysize + rhs.indexZ)
  )
  {
    return true;
  }
  else
  {
    return false;
  }

}