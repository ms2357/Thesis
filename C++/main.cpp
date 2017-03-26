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
    VectorXd initialVelocity(mesh.getPositionMesh().size());
    initialVelocity.setOnes();
    //VectorXd initialVelocity = VectorXd::Ones(mesh.getPositionMesh().size());
    UMatrix UEdge1(initialPressure, initialVelocity, N, L);
    UMatrix UEdge2(initialPressure, initialVelocity, N, L);

    /* Setup characteristic matrix */
    CharacteristicMatrix CharacteristicsEdge1(N);
    CharacteristicMatrix CharacteristicsEdge2(N);


    /*Run LxF method for solving interterior nodes */
    laxfriedrichs(mesh, waveSpeedDetails, UEdge1, UEdge2, CharacteristicsEdge1, CharacteristicsEdge2);
    std::cout << "Pressure matrix:" << std::endl;
    std::cout << UEdge1.getPressure() << std::endl;
    std::cout << "Velocity matrix:" << std::endl;
    std::cout << UEdge1.getVelocity() << std::endl;

    return 0;

}
