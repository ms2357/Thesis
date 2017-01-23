#ifndef MESH2D_H
#define MESH2D_H

#include <Eigen/Dense>

using namespace Eigen;

class Mesh2D
{
public:
    Mesh2D(double timeStart, double timeEnd, int timeLength,
           double positionStart, double positionEnd, int positionLength);
    VectorXd getTimeMesh() const;
    VectorXd getPositionMesh() const;
    double getTimeStep() const;
    double getPositionStep() const;

private:
    double timeStep;
    VectorXd timeMesh;
    double positionStep;
    VectorXd positionMesh;
};

#endif // MESH2D_H
