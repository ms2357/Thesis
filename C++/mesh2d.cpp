#include "mesh2d.h"

Mesh2D::Mesh2D(double positionStart, double positionEnd, int positionLength,
               double timeStart, double timeEnd)
{
    positionMesh = VectorXd::LinSpaced(positionLength + 1 , positionStart , positionEnd);
    positionStep = positionMesh(1)-positionMesh(0);

    //multiple of timeLength is the denom of timeStep
    timeLength = 2 * positionLength + 1;
    timeStep = positionStep / 2.0;
    timeMesh = VectorXd::LinSpaced(((timeEnd - timeStart) / timeStep) + 1,
                                   timeStart, timeStart + timeStep * timeLength);


}

VectorXd Mesh2D::getTimeMesh() const { return timeMesh; }

VectorXd Mesh2D::getPositionMesh() const { return positionMesh; }

double Mesh2D::getTimeStep() const { return timeStep; }

double Mesh2D::getPositionStep() const { return positionStep; }

double Mesh2D::getTimeLength() const { return timeLength; }
