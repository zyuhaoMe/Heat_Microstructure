/* 
 * Author: zyuhao 
 * Date: 2022-04-12
 * Content: This head defines the global constant
 * 
 */

//

#ifndef COMPARE_ANALYTICAL_FVM_CONSTANTS_H
#define COMPARE_ANALYTICAL_FVM_CONSTANTS_H

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

#endif //COMPARE_ANALYTICAL_FVM_CONSTANTS_H
