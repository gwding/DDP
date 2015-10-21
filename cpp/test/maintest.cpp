#include <iostream>
#include <fstream>

#include <Eigen/Dense>



using namespace std;
using namespace Eigen;

int main()
{
    Eigen::Matrix<int,4,1> A,B,test;
    A << 1,2,3,4;
    B << 2,1,4,3;

    cout << A << endl;
}
