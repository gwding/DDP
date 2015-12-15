#include "ilqrsolver.h"

/* Debug */
#include <iostream>
using namespace std;
/* */

using namespace Eigen;

ILQRSolver::ILQRSolver(DynamicModel& myDynamicModel, CostFunction& myCostFunction)
{
    dynamicModel = &myDynamicModel;
    costFunction = &myCostFunction;
    stateNb = myDynamicModel.getStateNb();
    commandNb = myDynamicModel.getCommandNb();

    // initializing matrix sizes
    /*x.resize(stateNb);
    u.resize(commandNb);
    xInit.resize(stateNb);
    xDes.resize(stateNb);


    nextVx.resize(stateNb);
    nextVxx.resize(stateNb,stateNb);
    Qx.resize(stateNb);
    Qxx.resize(stateNb,stateNb);
    Qu.resize(commandNb);
    Quu.resize(commandNb,commandNb);
    QuuInv.resize(commandNb,commandNb);
    Qux.resize(commandNb,stateNb);
    k.resize(stateNb);
    K.resize(stateNb,stateNb);*/

    muEye.resize(stateNb,stateNb);

    zeroCommand.resize(commandNb);
    zeroCommand.setZero();

    eye_stateSize.resize(stateNb,stateNb);
    eye_stateSize.setIdentity();
}

void ILQRSolver::FirstInitSolver(VectorXd myxInit, VectorXd myxDes, unsigned int& myT,
                       double& mydt, unsigned int& myiterMax,double& mystopCrit)
{
    xInit = myxInit;
    xDes = myxDes;
    T = myT;
    dt = mydt;
    iterMax = myiterMax;
    stopCrit = mystopCrit;

    xList = new VectorXd[myT+1];
    uList = new VectorXd[myT];
    updatedxList = new VectorXd[myT+1];
    updateduList = new VectorXd[myT];
    kList = new VectorXd[myT];
    KList = new MatrixXd[myT];
    tmpxPtr = new VectorXd[myT+1];
    tmpuPtr = new VectorXd[myT];
    lastTraj.xList = new VectorXd[myT+1];
    lastTraj.uList = new VectorXd[myT];
    /*for(int i=0;i<myT;i++)
    {
        xList[i].resize(stateNb);
        uList[i].resize(commandNb);
        updatedxList[i].resize(stateNb);
        updateduList[i].resize(commandNb);
        kList[i].resize(stateNb);
        KList[i].resize(commandNb,stateNb);
        tmpxPtr[i].resize(stateNb);
        tmpuPtr[i].resize(commandNb);
        lastTraj.xList[i].resize(stateNb);
        lastTraj.uList[i].resize(commandNb);
    }
    xList[myT].resize(stateNb);
    updatedxList[myT].resize(stateNb);
    tmpxPtr[myT].resize(stateNb);
    lastTraj.xList[myT].resize(stateNb);

    k.setZero();
    K.setZero();*/

    alphaList[0] = 1.0;
    alphaList[1] = 0.8;
    alphaList[2] = 0.6;
    alphaList[3] = 0.4;
    alphaList[4] = 0.2;
    alpha = 1.0;
}

void ILQRSolver::initSolver(VectorXd myxInit, VectorXd myxDes)
{
    xInit = myxInit;
    xDes = myxDes;
}

void ILQRSolver::solveTrajectory()
{
    initTrajectory();
    for(iter=0;iter<iterMax;iter++)
    {
        backwardLoop();
        forwardLoop();
        if(changeAmount<stopCrit)
            break;
        tmpxPtr = xList;
        tmpuPtr = uList;
        xList = updatedxList;
        updatedxList = tmpxPtr;
        uList = updateduList;
        updateduList = tmpuPtr;
    }
}

void ILQRSolver::initTrajectory()
{
    xList[0] = xInit;
    zeroCommand.setZero();
    for(unsigned int i=0;i<T;i++)
    {
        uList[i] = zeroCommand;
        xList[i+1] = dynamicModel->computeNextState(dt,xList[i],zeroCommand);
    }
}

void ILQRSolver::backwardLoop()
{
    costFunction->computeFinalCostDeriv(xList[T],xDes);
    nextVx = costFunction->getlx();
    nextVxx = costFunction->getlxx();

    mu = 0.0;
    completeBackwardFlag = 0;

    while(!completeBackwardFlag)
    {
        completeBackwardFlag = 1;
        muEye = mu*stateMat_t::Zero();
        for(int i=T-1;i>=0;i--)
        {
            x = xList[i];
            u = uList[i];

            dynamicModel->computeAllModelDeriv(dt,x,u);
            costFunction->computeAllCostDeriv(x,xDes,u);

            Qx = costFunction->getlx() + dynamicModel->getfx().transpose() * nextVx;
            Qu = costFunction->getlu() + dynamicModel->getfu().transpose() * nextVx;
            Qxx = costFunction->getlxx() + dynamicModel->getfx().transpose() * (nextVxx+muEye) * dynamicModel->getfx();
            Quu = costFunction->getluu() + dynamicModel->getfu().transpose() * (nextVxx+muEye) * dynamicModel->getfu();
            Qux = costFunction->getlux() + dynamicModel->getfu().transpose() * (nextVxx+muEye) * dynamicModel->getfx();

            Qxx += dynamicModel->computeTensorContxx(nextVx);
            Qux += dynamicModel->computeTensorContux(nextVx);
            Quu += dynamicModel->computeTensorContuu(nextVx);


            if(!isQuudefinitePositive(Quu))
            {
                /*
                  To be Implemented : Regularization (is Quu definite positive ?)
                */
                if(mu==0.0) mu += 1e-4;
                else mu *= 10;
                completeBackwardFlag = 0;
                break;
            }

            QuuInv = Quu.inverse();
            k = -QuuInv*Qu;
            K = -QuuInv*Qux;

            nextVx = Qx - K.transpose()*Quu*k;
            nextVxx = Qxx - K.transpose()*Quu*K;

            kList[i] = k;
            KList[i] = K;
        }
    }
}

void ILQRSolver::forwardLoop()
{
    changeAmount = 0.0;
    updatedxList[0] = xInit;
    // Line search to be implemented
    alpha = 1.0;
    for(unsigned int i=0;i<T;i++)
    {
        updateduList[i] = uList[i] + alpha*kList[i] + KList[i]*(updatedxList[i]-xList[i]);
        updatedxList[i+1] = dynamicModel->computeNextState(dt,updatedxList[i],updateduList[i]);
        for(unsigned int j=0;j<commandNb;j++)
        {
            changeAmount += abs(uList[i](j,0) - updateduList[i](j,0));
        }
    }
}

ILQRSolver::traj ILQRSolver::getLastSolvedTrajectory()
{
    lastTraj.xList = updatedxList;
    lastTraj.uList = updateduList;
    lastTraj.iter = iter;
    return lastTraj;
}

char ILQRSolver::isQuudefinitePositive(MatrixXd& Quu)
{
    /*
      Todo : check if Quu is definite positive
    */
    return 1;
}
