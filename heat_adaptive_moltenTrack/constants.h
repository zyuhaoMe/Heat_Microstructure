//
// Created by zyuhao on 5/31/22.
// store the global constants

#ifndef THERMAL_MULTI_MPI_CONSTANTS_H
#define THERMAL_MULTI_MPI_CONSTANTS_H
#include <cmath>



/* [0]
 * constants "inherited" from compare_analytical_FVM,
 * used in LaserClass and MultiTrackLaserClass
 */

// pi
const double PI = 3.14159265359;

// index of scanning path: stationary scanning
const int INDEX_SCANNING_PATH_STATIONARY = 0;

// index of scanning path: linear scanning
const int INDEX_SCANNING_PATH_LINEAR = 1;

// index of scanning path: spiral scanning
const int INDEX_SCANNING_PATH_SPIRAL = 2;

// index of scanning path: multi single track scanning
const int INDEX_SCANNING_PATH_MULTI_LINE = 3;

// index of line segment in multi line meaning laser is off
const int INDEX_LASER_IS_OFF = -1;

// time duration of active laser heat source during stationary scanning path
const double TIME_LASER_ACTIVE_STATIONARY = 5e-3;

// time duration of active laser heat source during liear scaning path
const double TIME_LASER_ACTIVE_LINEAR = 1.25e-3;

// default ambient temperature
const double TEMP_AMBIENT_DEFAULT = 300;


/* [0.1] parameter of test cube */

// 3d dimension
const double XLEN_CUBE_TEST = 2e-3;
const double YLEN_CUBE_TEST = 2e-3;
const double ZLEN_CUBE_TEST = 1e-3;
const double MESHSIZE_CUBE_TEST = 1e-5;

// total time for the test cube case
const double TOTAL_TIME_CUBE_TEST = 6e-3;

/*============================================================================*/

/* [1]
 * constants relative to the geometry of 3d build
 * 4 repeated blocks w.r.t layer num + 1 unchanged block
 */


/* [1.1] geometry of one repeated block */

// length of leg 1 (mm)
const double LEN_LEG_1 = 5e-3;

// length of leg 2
const double LEN_LEG_2 = 0.5e-3;

// length of leg 3
const double LEN_LEG_3 = 2.5e-3;

// length of gap
const double LEN_GAP = 2e-3;

// total length of one repeated block
const double LEN_RPT_BLOCK = LEN_LEG_1 + LEN_LEG_2 + LEN_LEG_3 + 3 * LEN_GAP;

// height of leg with unchanged leg length
const double HT_UNCH_LEG_LEN = 5e-3;

// height of leg with changing leg length
const double HT_CHG_LEG_LEN = 2e-3;

// height of block part with no leg
const double HT_NO_LEG_LEN = 5e-3;

// sin(45) and cos(45)



/* [1.2] shared geometry parameter of all 5 blocks */

// width of each block (same for all block)
const double WD_BLOCK = 5e-3;

// height of each block (same for all block)
const double HT_BLOCK = HT_UNCH_LEG_LEN + HT_CHG_LEG_LEN + HT_NO_LEG_LEN;



/* [1.3] geometry of unchanged unique block */

// total length of one unique block
const double LEN_UNIQUE_BLOCK = 19e-3;

// length of chamfer with leading angles of 45 degrees
const double LEN_CHAMFER = WD_BLOCK / 2;

// length of unique block with on chamfer
const double LEN_UNIQUE_BLOCK_NO_CHAMFER = LEN_UNIQUE_BLOCK - LEN_CHAMFER;

// [1.4] absolute x position of blockA, blockB, blockC, blockD, blockE
const double X_POS_ABS_BLOCK_A = 0.;
const double X_POS_ABS_BLOCK_B = LEN_RPT_BLOCK;
const double X_POS_ABS_BLOCK_C = LEN_RPT_BLOCK * 2.;
const double X_POS_ABS_BLOCK_D = LEN_RPT_BLOCK * 3.;
const double X_POS_ABS_BLOCK_E = LEN_RPT_BLOCK * 4.;


/*============================================================================*/


/* [2]
 * constants relative to scanning strategy: track geometry and laser parameters
 */


/* [2.1] laser and geometry parameters shared by both odd and even layers */

// hatch of both odd and even layers
const double HATCH_SPACING = 0.1e-3;

// thickness of each layer
const double LAYER_THK = 0.02e-3;

// number of same parts being manufactured, effect the thermal history
const int NUM_RPT_PARTS = 4;

// speed of laser during contouring for all layers
const double SPEED_CTR = 0.9;

// speed of laser during infilling for all layers
const double SPEED_FILL = 0.8;

// power of laser during contouring for all layers
const double PWR_CTR = 1e2;

// power of laser during infilling for all layers
const double PWR_FILL = 195.;



/* [2.2] contouring scan timing for each part */

/* [2.2.1] layer 1-250 laser-on duration */

// L1, L4, L5, L10, constant contour length

/* [2.2] odd layer scanning parameters */


// even layer scanning parameters

#endif //THERMAL_MULTI_MPI_CONSTANTS_H
