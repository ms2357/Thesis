#include <iostream>
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/LU>

#include "laxfriedrichs.h"
#include "mesh2d.h"
#include "wavespeed.h"
#include "umatrix.h"
#include "characteristicmatrix.h"

int main()
{
    using namespace Eigen;

    /* Oval problem setup */
    /* Setup the mesh */
    double t0 = 0.0, tf = 1.0; int N = 4; // time interval
    double x0 = 0.0, xf = 1.0; int L = 4; // space interval
    Mesh2D mesh(t0, tf, N, x0, xf, L);

    /* Wave speed details */
    double rho0 = 1.0, K0 = 1.0;
    WaveSpeed waveSpeedDetails(rho0, K0);

    /* Setup matrix U */
    VectorXd initialPressure = (M_PI * mesh.getPositionMesh()).array().sin();
    VectorXd initialVelocity = (M_PI * mesh.getPositionMesh()).array().cos();
    UMatrix U(initialPressure, initialVelocity, N, L);

    /* Impose boundary conditions and setup characteristic matrix */
    MatrixXd bd0(2, mesh.getTimeMesh().size());
    bd0.row(0) = (M_PI * mesh.getTimeMesh()).array().sin();
    bd0.row(1) = (M_PI * mesh.getTimeMesh()).array().cos();
    bd0 = waveSpeedDetails.getR() * bd0;
    MatrixXd bdL(2, mesh.getTimeMesh().size());
    bdL.row(0) = (M_PI * (L - mesh.getTimeMesh().array())).sin();
    bdL.row(1) = (M_PI * (L - mesh.getTimeMesh().array())).cos();
    bdL = waveSpeedDetails.getR() * bdL;
    CharacteristicMatrix C(bd0, bdL);

    laxfriedrichs(mesh, waveSpeedDetails, U, C);
    std::cout << "Pressure matrix:" << std::endl;
    std::cout << U.getPressure() << std::endl;
    std::cout << "Velocity matrix:" << std::endl;
    std::cout << U.getVelocity() << std::endl;

    return 0;

}
