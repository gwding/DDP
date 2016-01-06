#ifndef COSTFUNCTIONPNEUMATICARMELBOW_H
#define COSTFUNCTIONPNEUMATICARMELBOW_H

#include "config.h"

#include "costfunction.h"

#include <Eigen/Dense>

using namespace Eigen;

class CostFunctionPneumaticarmElbow : public CostFunction
{
public:
    CostFunctionPneumaticarmElbow();
private:
    stateMat_t Q, Qf;
    commandMat_t R;
    stateVec_t lx;
    stateMat_t lxx;
    commandVec_t lu;
    commandMat_t luu;
    commandR_stateC_t lux;
    stateR_commandC_t lxu;
    double dt;
protected:
    // attributes //
public:
private:

protected:
    // methods //
public:
    void computeAllCostDeriv(const stateVec_t& X, const commandVec_t& U);
    void computeFinalCostDeriv(const stateVec_t& X);
private:
protected:
    // accessors //
public:
    stateVec_t getlx();
    stateMat_t getlxx();
    commandVec_t getlu();
    commandMat_t getluu();
    commandR_stateC_t getlux();
    stateR_commandC_t getlxu();
};

#endif
