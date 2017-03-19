#include "wavespeed.h"
/*setup and store all of the details related to the wave speed, impedence, and coefficient matrices
 * used to solve the system of equations governed by continuity of pressure and conservation of momentum*/
WaveSpeed::WaveSpeed(double density1, double bulkModulus1, double density2, double bulkModulus2)
{
    double impedence1 = density1 * bulkModulus1;
    double impedence2 = density2 * bulkModulus2;
    //this may need to be changed to accomidate for diff A????
    A << 0, bulkModulus1,
         1/density1, 0;
    //this may need to be changed to accomidate for diff R????
    R << -impedence1, impedence1,
                  1,         1;

    RI << R.inverse();

    //coefficient matrix to solve for unk characteristics
    Z << impedence2, impedence1,
                 -1,          1;

    ZcoefficientMatrix << impedence2, impedence1,
                                  -1,          1;

    VectorXd initwavespeed = A.eigenvalues().real();

    // Swap these if necessary
    c0w = initwavespeed(1);
    c0z = initwavespeed(0);

}

Matrix2d WaveSpeed::getA() const { return A; }

Matrix2d WaveSpeed::getR() const { return R; }

Matrix2d WaveSpeed::getRI() const { return RI; }

Matrix2d WaveSpeed::getZ() const { return Z; }

Matrix2d WaveSpeed::getZcoefficientMatrix() const { return ZcoefficientMatrix; }

double WaveSpeed::getc0w() const { return c0w; }

double WaveSpeed::getc0z() const { return c0z; }


