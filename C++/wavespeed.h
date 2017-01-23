#ifndef WAVESPEED_H
#define WAVESPEED_H

#include <Eigen/Dense>

using namespace Eigen;

class WaveSpeed
{
public:
    WaveSpeed(double density, double bulkModulus);
    Matrix2d getA() const;
    Matrix2d getR() const;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
    Matrix2d A;
    Matrix2d R;
};

#endif // WAVESPEED_H
