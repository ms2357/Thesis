#include "umatrix.h"
/*set up for the U matrix which holds the pressure and velocity on each edge, fucntion setPoint allows
for easy updating of each node, interior and boundary*/
UMatrix::UMatrix(VectorXd initialPressure, VectorXd initialVelocity, int N, int L)
{
    pressure.resize(N + 1, L + 1);
    pressure.row(0) = initialPressure;
    velocity.resize(N + 1, L + 1);
    velocity.row(0) = initialVelocity;
}

MatrixXd UMatrix::getPressure() { return pressure; }

MatrixXd UMatrix::getVelocity() { return velocity; }

void UMatrix::setPoint(int i, int j, Vector2d point)
{
    pressure(i, j) = point(0);
    velocity(i, j) = point(1);
}
