#ifndef CONFIG_H
#define CONFIG_H

#include <Eigen/Dense>

#define stateSize_t 4
#define commandSize_t 1

// typedef for stateSize types
typedef Eigen::Matrix<double,stateSize_t,1> stateVec_t;                       // stateSize x 1
typedef Eigen::Matrix<double,1,stateSize_t> stateVecTrans_t;                  // 1 x stateSize
typedef Eigen::Matrix<double,stateSize_t,stateSize_t> stateMat_t;               // stateSize x stateSize
typedef Eigen::Matrix<double,stateSize_t,stateSize_t> stateTens_t[stateSize_t];   // stateSize x stateSize x stateSize

// typedef for commandSize types
typedef Eigen::Matrix<double,commandSize_t,1> commandVec_t;                           // commandSize x 1
typedef Eigen::Matrix<double,1,commandSize_t> commandVecTrans_t;                      // 1 x commandSize
typedef Eigen::Matrix<double,commandSize_t,commandSize_t> commandMat_t;                 // commandSize x commandSize
typedef Eigen::Matrix<double,commandSize_t,commandSize_t> commandTens_t[commandSize_t];   // stateSize x commandSize x commandSize



// typedef for mixed stateSize and commandSize types
typedef Eigen::Matrix<double,stateSize_t,commandSize_t> stateR_commandC_t;                          // stateSize x commandSize
typedef Eigen::Matrix<double,stateSize_t,commandSize_t> stateR_commandC_stateD_t[stateSize_t];        // stateSize x commandSize x stateSize
typedef Eigen::Matrix<double,stateSize_t,commandSize_t> stateR_commandC_commandD_t[commandSize_t];    // stateSize x commandSize x commandSize
typedef Eigen::Matrix<double,commandSize_t,stateSize_t> commandR_stateC_t;                          // commandSize x stateSize
typedef Eigen::Matrix<double,commandSize_t,stateSize_t> commandR_stateC_stateD_t[stateSize_t];        // commandSize x stateSize x stateSize
typedef Eigen::Matrix<double,commandSize_t,stateSize_t> commandR_stateC_commandD_t[commandSize_t];    // commandSize x stateSize x commandSize
typedef Eigen::Matrix<double,stateSize_t,stateSize_t> stateR_stateC_commandD_t[commandSize_t];    // stateSize x stateSize x commandSize
typedef Eigen::Matrix<double,commandSize_t,commandSize_t> commandR_commandC_stateD_t[stateSize_t];    // commandSize x commandSize x stateSize


#endif // CONFIG_H
