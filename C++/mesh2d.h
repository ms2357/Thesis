#ifndef MESH2D_H
#define MESH2D_H

#include <Eigen/Dense>

using namespace Eigen;

class Mesh2D
{
public:
    Mesh2D(double positionStart, double positionEnd, int positionLength,
           double timeStart, double timeEnd);

    VectorXd getTimeMesh() const;
    VectorXd getPositionMesh() const;
    double getTimeStep() const;
    double getPositionStep() const;
    double getTimeLength() const;

private:
    double timeStep;
    VectorXd timeMesh;
    double positionStep;
    VectorXd positionMesh;
    double timeLength;
};

#endif // MESH2D_H
