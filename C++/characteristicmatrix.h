#ifndef CHARACTERISTICMATRIX_H
#define CHARACTERISTICMATRIX_H

#include <Eigen/Dense>

using namespace Eigen;

class CharacteristicMatrix
{
public:
    CharacteristicMatrix(int N);
    double getW(int i, int j);
    double getZ(int i, int j);
    void setW(int i, int j, double value);
    void setZ(int i, int j, double value);


private:
    MatrixXd w;
    MatrixXd z;

};

#endif // CHARACTERISTICMATRIX_H
