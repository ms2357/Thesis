#ifndef LAXFRIEDRICHS
#define LAXFRIEDRICHS

#include <Eigen/Dense>
#include <Eigen/LU>

#include "umatrix.h"
#include "characteristicmatrix.h"
#include "mesh2d.h"
#include "wavespeed.h"

using namespace Eigen;

void laxfriedrichs(const Mesh2D& mesh, const WaveSpeed& waveSpeedDetails, UMatrix& UEdge1, UMatrix& UEdge2,
                   CharacteristicMatrix& CharacteristicsEdge1, CharacteristicMatrix& CharacteristicsEdge2);

MatrixXd calculateCFL(const Mesh2D& mesh, const Matrix2d& A);

void updateInteriorPoint(int i, int j, UMatrix& U, const MatrixXd& mu);

void calcKnownCharacteristics(int i, int L, UMatrix& U, CharacteristicMatrix& Characteristics,
                              const Matrix2d& R, double c0w, double c0z, double TimeStep, double PositionStep);

void JunctionConditions(int i, UMatrix & UEdge1, CharacteristicMatrix& CharacteristicsEdge1, CharacteristicMatrix& CharacteristicsEdge2,
                        const WaveSpeed& waveSpeedDetails);

void updateBoundaryPoints(int i, int L, UMatrix& UEdge1, MatrixXd CharacteristicMatrixEdge1,
                          MatrixXd CharacteristicMatrixEdge2, const MatrixXd& R, const MatrixXd& RI,
                          const MatrixXd& mu);


#endif // LAXFRIEDRICHS

