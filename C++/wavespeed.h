#ifndef WAVESPEED_H
#define WAVESPEED_H

#include <Eigen/Dense>

using namespace Eigen;

class WaveSpeed
{
public:
    WaveSpeed(double density1, double bulkModulus1, double density2, double bulkModulus2);
    Matrix2d getA() const;
    Matrix2d getR() const;
    Matrix2d getRI() const;
    Matrix2d getZ() const;
    Matrix2d getZcoefficientMatrix() const;
    double getc0w() const;
    double getc0z() const;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
    Matrix2d A;
    Matrix2d R;
    Matrix2d RI;
    Matrix2d Z;
    Matrix2d ZcoefficientMatrix;
    double c0w;
    double c0z;
};

#endif // WAVESPEED_H
