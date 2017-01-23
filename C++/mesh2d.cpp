#include "mesh2d.h"

Mesh2D::Mesh2D(double timeStart, double timeEnd, int timeLength,
               double positionStart, double positionEnd, int positionLength)
{
    timeStep = 1.0 / timeLength;
    timeMesh = VectorXd::LinSpaced(((timeEnd - timeStart) / timeStep) + 1,
                                   timeStart, timeStart + timeStep * timeLength);

    positionStep = 1.0 / positionLength;
    positionMesh = VectorXd::LinSpaced(((positionEnd - positionStart) / positionStep) + 1,
                                       positionStart, positionStart + positionStep * positionLength);
}

VectorXd Mesh2D::getTimeMesh() const { return timeMesh; }

VectorXd Mesh2D::getPositionMesh() const { return positionMesh; }

double Mesh2D::getTimeStep() const { return timeStep; }

double Mesh2D::getPositionStep() const { return positionStep; }

