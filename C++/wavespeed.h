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
    Matrix2d getRI() const;
    double getc0w() const;
    double getc0z() const;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
    Matrix2d A;
    Matrix2d R;
    Matrix2d RI;
    double c0w;
    double c0z;
};

#endif // WAVESPEED_H
