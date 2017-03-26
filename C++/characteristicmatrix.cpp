#include "characteristicmatrix.h"
/*defines the matrix that hold the value of the characteristic variables on the boundaries,
 * allows updates via set functions*/

CharacteristicMatrix::CharacteristicMatrix(int N)
{
    w = MatrixXd::Zero(N + 1, 2);
    z = MatrixXd::Zero(N + 1, 2);
}

double CharacteristicMatrix::getW(int i, int j) { return w(i, j); }

double CharacteristicMatrix::getZ(int i, int j) { return z(i, j); }

void CharacteristicMatrix::setW(int i, int j, double value)
{
    w(i, j) = value;
}

void CharacteristicMatrix::setZ(int i, int j, double value)
{
    z(i, j) = value;
}

