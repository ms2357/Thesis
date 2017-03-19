#include "characteristicmatrix.h"
/*defines the matrix that hold the value of the characteristic variables on the boundaries,
 * allows updates via set functions*/

CharacteristicMatrix::CharacteristicMatrix(MatrixXd bd0, MatrixXd bdL)
{
    w.resize(bd0.cols(), 2);
    w.col(0) = bd0.row(0);
    z.resize(bdL.cols(), 2);
    z.col(1) = bdL.row(1);
}

MatrixXd CharacteristicMatrix::getW() { return w; }

MatrixXd CharacteristicMatrix::getZ() { return z; }

void CharacteristicMatrix::setW(int i, int j, double value)
{
    w(i, j) = value;
}

void CharacteristicMatrix::setZ(int i, int j, double value)
{
    z(i, j) = value;
}

