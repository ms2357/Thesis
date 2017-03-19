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

    /* Oval problem setup, initial parameters */
    /* Setup the mesh */
    double t0 = 0.0, tf = 1.0;             // time interval
    double x0 = 0.0, xf = 1.0; int L = 4; // space interval
    Mesh2D mesh(x0, xf, L,t0, tf);
    double N = mesh.getTimeLength()-1;     //number of time intervals


    /* Wave speed details */
    double rho1 = 1.0, K1 = 1.0, rho2 = 1.0, K2 = 1.0;
    WaveSpeed waveSpeedDetails(rho1, K1, rho2, K2);

    /* Setup matrix U, holds pressure and velocity */
    VectorXd initialPressure = (2 * M_PI * mesh.getPositionMesh()).array().sin();
    VectorXd initialVelocity = MatrixXd::Ones(1 , mesh.getPositionMesh().size());
    UMatrix UEdge1(initialPressure, initialVelocity, N, L);
    UMatrix UEdge2(initialPressure, initialVelocity, N, L);

    /* Setup characteristic matrix */
    MatrixXd CharacteristicMatrixEdge1W = MatrixXd::Zero(N,2);
    MatrixXd CharacteristicMatrixEdge1Z = MatrixXd::Zero(N,2);
    MatrixXd CharacteristicMatrixEdge2W = MatrixXd::Zero(N,2);
    MatrixXd CharacteristicMatrixEdge2Z = MatrixXd::Zero(N,2);

    /*Run LxF method for solving interterior nodes */
    laxfriedrichs(mesh, waveSpeedDetails, UEdge1, UEdge2, CharacteristicMatrixEdge1W, CharacteristicMatrixEdge1Z, CharacteristicMatrixEdge2W, CharacteristicMatrixEdge2Z);
    //std::cout << "Pressure matrix:" << std::endl;
    //std::cout << U.getPressure() << std::endl;
    //std::cout << "Velocity matrix:" << std::endl;
    //std::cout << U.getVelocity() << std::endl;

    return 0;

}
