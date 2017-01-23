#ifndef UMATRIX_H
#define UMATRIX_H

#include <Eigen/Dense>

using namespace Eigen;

class UMatrix
{
public:
    UMatrix(VectorXd initialPressure, VectorXd initialVelocity, int N, int L);
    MatrixXd getPressure();
    MatrixXd getVelocity();
    void setPoint(int i, int j, Vector2d point);

private:
    MatrixXd pressure;
    MatrixXd velocity;
};

#endif // UMATRIX_H
