#ifndef LAXFRIEDRICHS
#define LAXFRIEDRICHS

#include <Eigen/Dense>
#include <Eigen/LU>

#include "umatrix.h"
#include "characteristicmatrix.h"
#include "mesh2d.h"
#include "wavespeed.h"

using namespace Eigen;

void laxfriedrichs(const Mesh2D& mesh, const WaveSpeed& waveSpeedDetails, UMatrix& U, CharacteristicMatrix& C);
MatrixXd calculateCFL(const Mesh2D& mesh, const Matrix2d& A);
void updateInteriorPoint(int i, int j, UMatrix& U, const MatrixXd& mu);
void updateBoundaryPoints(int i, int L, UMatrix& U, CharacteristicMatrix& C, const MatrixXd& R, const MatrixXd& mu);


#endif // LAXFRIEDRICHS

