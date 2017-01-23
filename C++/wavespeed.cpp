#include "wavespeed.h"

WaveSpeed::WaveSpeed(double density, double bulkModulus)
{
    double impedence = density * bulkModulus;

    A << 0, bulkModulus,
         1/density, 0;

    R << -impedence, impedence,
                  1,         1;
}

Matrix2d WaveSpeed::getA() const { return A; }

Matrix2d WaveSpeed::getR() const { return R; }

