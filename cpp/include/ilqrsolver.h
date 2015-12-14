#ifndef ILQRSOLVER_H
#define ILQRSOLVER_H

#include "dynamicmodel.h"
#include "costfunction.h"
#include <Eigen/Dense>

using namespace Eigen;

template<typename precision,int stateSize,int commandSize>
class ILQRSolver
{
    // typedef for stateSize types
    typedef Eigen::Matrix<precision,stateSize,1> stateVec_t;                       // stateSize x 1
    typedef Eigen::Matrix<precision,1,stateSize> stateVecTrans_t;                  // 1 x stateSize
    typedef Eigen::Matrix<precision,stateSize,stateSize> stateMat_t;               // stateSize x stateSize
    typedef Eigen::Matrix<precision,stateSize,stateSize> stateTens_t[stateSize];   // stateSize x stateSize x stateSize

    // typedef for commandSize types
    typedef Eigen::Matrix<precision,commandSize,1> commandVec_t;                           // commandSize x 1
    typedef Eigen::Matrix<precision,1,commandSize> commandVecTrans_t;                      // 1 x commandSize
    typedef Eigen::Matrix<precision,commandSize,commandSize> commandMat_t;                 // commandSize x commandSize
    typedef Eigen::Matrix<precision,commandSize,commandSize> commandTens_t[commandSize];   // stateSize x commandSize x commandSize



    // typedef for mixed stateSize and commandSize types
    typedef Eigen::Matrix<precision,stateSize,commandSize> stateR_commandC_t;                          // stateSize x commandSize
    typedef Eigen::Matrix<precision,stateSize,commandSize> stateR_commandC_stateD_t[stateSize];        // stateSize x commandSize x stateSize
    typedef Eigen::Matrix<precision,stateSize,commandSize> stateR_commandC_commandD_t[commandSize];    // stateSize x commandSize x commandSize
    typedef Eigen::Matrix<precision,commandSize,stateSize> commandR_stateC_t;                          // commandSize x stateSize
    typedef Eigen::Matrix<precision,commandSize,stateSize> commandR_stateC_stateD_t[stateSize];        // commandSize x stateSize x stateSize
    typedef Eigen::Matrix<precision,commandSize,stateSize> commandR_stateC_commandD_t[commandSize];    // commandSize x stateSize x commandSize
    typedef Eigen::Matrix<precision,stateSize,stateSize> stateR_stateC_commandD_t[commandSize];    // stateSize x stateSize x commandSize
    typedef Eigen::Matrix<precision,commandSize,commandSize> commandR_commandC_stateD_t[stateSize];    // commandSize x commandSize x stateSize

public:
    struct traj
    {
        stateVec_t* xList;
        commandVec_t* uList;
        unsigned int iter;
    };

public:
private:
    DynamicModel* dynamicModel;
    CostFunction* costFunction;
    unsigned int stateNb;
    unsigned int commandNb;
    stateVec_t x;
    commandVec_t u;
    stateVec_t xInit;
    stateVec_t xDes;
    unsigned int T;
    unsigned int iter;
    double dt;
    unsigned int iterMax;
    precision stopCrit;
    precision changeAmount;

    stateVec_t* xList;
    commandVec_t* uList;
    stateVec_t* updatedxList;
    commandVec_t* updateduList;
    struct traj lastTraj;

    stateVec_t nextVx;
    stateMat_t nextVxx;
    stateVec_t Qx;
    stateMat_t Qxx;
    commandVec_t Qu;
    commandMat_t Quu;
    commandMat_t QuuInv;
    commandR_stateC_t Qux;
    commandVec_t k;
    commandR_stateC_t K;
    commandVec_t* kList;
    commandR_stateC_t* KList;
    precision alphaList[5];
    precision alpha;

    precision mu;
    stateMat_t muEye;
    unsigned char completeBackwardFlag;

public:
    ILQRSolver(DynamicModel& myDynamicModel, CostFunction& myCostFunction)
    {
        dynamicModel = &myDynamicModel;
        costFunction = &myCostFunction;
        stateNb = myDynamicModel.getStateNb();
        commandNb = myDynamicModel.getCommandNb();
    }

    void FirstInitSolver(stateVec_t& myxInit, stateVec_t& myxDes, unsigned int& myT,
                           precision& mydt, unsigned int& myiterMax,precision& mystopCrit)
    {
        xInit = myxInit;
        xDes = myxDes;
        T = myT;
        dt = mydt;
        iterMax = myiterMax;
        stopCrit = mystopCrit;

        xList = new stateVec_t[myT+1];
        uList = new commandVec_t[myT];
        updatedxList = new stateVec_t[myT+1];
        updateduList = new commandVec_t[myT];
        k.setZero();
        K.setZero();
        kList = new commandVec_t[myT];
        KList = new commandR_stateC_t[myT];

        alphaList[0] = 1.0;
        alphaList[1] = 0.8;
        alphaList[2] = 0.6;
        alphaList[3] = 0.4;
        alphaList[4] = 0.2;
        alpha = 1.0;
    }

    void initSolver(stateVec_t& myxInit, stateVec_t& myxDes)
    {
        xInit = myxInit;
        xDes = myxDes;
    }

    void solveTrajectory()
    {
        stateVec_t* tmpxPtr;
        commandVec_t* tmpuPtr;
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

    void initTrajectory()
    {
        xList[0] = xInit;
        commandVec_t zeroCommand;
        zeroCommand.setZero();
        for(unsigned int i=0;i<T;i++)
        {
            uList[i] = zeroCommand;
            xList[i+1] = dynamicModel->computeNextState(dt,xList[i],zeroCommand);
        }
    }

    void backwardLoop()
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

    void forwardLoop()
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

    ILQRSolver::traj getLastSolvedTrajectory()
    {
        lastTraj.xList = updatedxList;
        lastTraj.uList = updateduList;
        lastTraj.iter = iter;
        return lastTraj;
    }

    char isQuudefinitePositive(commandMat_t& Quu)
    {
        /*
          Todo : check if Quu is definite positive
        */
        return 1;
    }


};

#endif // ILQRSOLVER_H
