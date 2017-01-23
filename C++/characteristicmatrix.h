#ifndef CHARACTERISTICMATRIX_H
#define CHARACTERISTICMATRIX_H

#include <Eigen/Dense>

using namespace Eigen;

class CharacteristicMatrix
{
public:
    CharacteristicMatrix(MatrixXd bd0, MatrixXd bdL);
    MatrixXd getW();
    MatrixXd getZ();
    void setW(int i, int j, double value);
    void setZ(int i, int j, double value);

private:
    MatrixXd w;
    MatrixXd z;
};

#endif // CHARACTERISTICMATRIX_H
