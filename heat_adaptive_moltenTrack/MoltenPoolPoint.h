//
// Created by zyuhao on 6/9/22.
// Information representing a point in the melt pool
//

#ifndef THERMAL_MULTI_MPI_MOLTENPOOLPOINT_H
#define THERMAL_MULTI_MPI_MOLTENPOOLPOINT_H

class MoltenPoolPoint {

  public:
    int indexX, indexY, indexZ; // int index in 3d coordinates

    static int xsize;
    static int ysize;
    static int zsize;

    // default ctor
    MoltenPoolPoint() = default;

    // ctor
    MoltenPoolPoint(
         int inIndexX,
         int inIndexY,
         int inIndexZ
         );

    // dtor
    ~MoltenPoolPoint() = default;

    // static setter of x y z size before any object has been created
   // static int setXYZsize

    // overwrite == for std::unique
    bool operator == (const MoltenPoolPoint& rhs) const;

    // overwrite < for std::sort
    bool operator < (const MoltenPoolPoint& rhs) const;

};


#endif //THERMAL_MULTI_MPI_MOLTENPOOLPOINT_H
