#include "wavespeed.h"

WaveSpeed::WaveSpeed(double density, double bulkModulus)
{
    double impedence = density * bulkModulus;

    A << 0, bulkModulus,
         1/density, 0;

    R << -impedence, impedence,
                  1,         1;

    RI << R.inverse();

    VectorXd initwavespeed = A.eigenvalues().real();

    // Swap these if necessary
    c0w = initwavespeed(1);
    c0z = initwavespeed(0);

}

Matrix2d WaveSpeed::getA() const { return A; }

Matrix2d WaveSpeed::getR() const { return R; }

Matrix2d WaveSpeed::getRI() const { return RI; }

double WaveSpeed::getc0w() const { return c0w; }

double WaveSpeed::getc0z() const { return c0z; }


