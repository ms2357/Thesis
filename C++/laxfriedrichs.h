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
                   MatrixXd CharacteristicMatrixEdge1W, MatrixXd CharacteristicMatrixEdge1Z,
                   MatrixXd CharacteristicMatrixEdge2W, MatrixXd CharacteristicMatrixEdge2z);

MatrixXd calculateCFL(const Mesh2D& mesh, const Matrix2d& A);

void updateInteriorPoint(int i, int j, UMatrix& U, const MatrixXd& mu);

void calculateKnownCharacteristics(int i, int L, UMatrix& U, MatrixXd& CharacteristicMatrixW,
                              MatrixXd& CharacteristicMatrixZ, const MatrixXd& R,
                              double c0w, double c0z, double TimeStep, double PositionStep);

void calculateUnkownCharacteristics(int i, UMatrix & U, MatrixXd& CharacteristicMatrixW1, MatrixXd& CharacteristicMatrixZ1,
                                    MatrixXd& CharacteristicMatrixW2, MatrixXd& CharacteristicMatrixZ2,
                                    const MatrixXd& Z, const MatrixXd& ZcoefficientMatrix, const MatrixXd& RI);

void updateBoundaryPoints(int i, int L, UMatrix& UEdge1, UMatrix& UEdge2, MatrixXd CharacteristicMatrixEdge1,
                          MatrixXd CharacteristicMatrixEdge2, const MatrixXd& R, const MatrixXd& RI,
                          const MatrixXd& mu);


#endif // LAXFRIEDRICHS

