#ifndef ILQRSOLVER_H
#define ILQRSOLVER_H

#include "dynamicmodel.h"
#include "costfunction.h"
#include <Eigen/Dense>

using namespace Eigen;

class ILQRSolver
{
public:
    struct traj
    {
        VectorXd* xList;
        VectorXd* uList;
        unsigned int iter;
    };

public:
    ILQRSolver(DynamicModel& myDynamicModel, CostFunction& myCostFunction);
private:
protected:
    // attributes //
public:
private:
    DynamicModel* dynamicModel;
    CostFunction* costFunction;
    unsigned int stateNb;
    unsigned int commandNb;
    VectorXd x;
    VectorXd u;
    VectorXd xInit;
    VectorXd xDes;
    unsigned int T;
    unsigned int iter;
    double dt;
    unsigned int iterMax;
    double stopCrit;
    double changeAmount;

    VectorXd* xList;
    VectorXd* uList;
    VectorXd* updatedxList;
    VectorXd* updateduList;
    struct traj lastTraj;
    VectorXd* tmpxPtr;
    VectorXd* tmpuPtr;
    VectorXd zeroCommand;

    VectorXd nextVx;
    MatrixXd nextVxx;
    VectorXd Qx;
    MatrixXd Qxx;
    VectorXd Qu;
    MatrixXd Quu;
    MatrixXd QuuInv;
    MatrixXd Qux;
    VectorXd k;
    MatrixXd K;
    VectorXd* kList;
    MatrixXd* KList;
    double alphaList[5];
    double alpha;

    MatrixXd eye_stateSize;


    double mu;
    MatrixXd muEye;
    unsigned char completeBackwardFlag;

protected:
    // methods //
public:
    void FirstInitSolver(VectorXd myxInit, VectorXd myxDes, unsigned int& myT,
                    double& mydt, unsigned int& myiterMax,double& mystopCrit);
    void initSolver(VectorXd myxInit, VectorXd myxDes);
    void solveTrajectory();
    void initTrajectory();
    void backwardLoop();
    void forwardLoop();
    char isQuudefinitePositive(MatrixXd& Quu);
    struct traj getLastSolvedTrajectory();
private:
protected:

};

#endif // ILQRSOLVER_H
