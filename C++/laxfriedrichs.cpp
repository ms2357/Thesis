#include "laxfriedrichs.h"

void laxfriedrichs(const Mesh2D& mesh, const WaveSpeed& waveSpeedDetails, UMatrix& U, CharacteristicMatrix& C)
{
    MatrixXd mu = calculateCFL(mesh, waveSpeedDetails.getA());
    int N = mesh.getTimeMesh().rows() - 1;
    int L = mesh.getPositionMesh().rows()- 1;

    for (int i = 0; i < N; i++) {
        for (int j = 1; j < L; j++) {
            updateInteriorPoint(i, j, U, mu);
        }
        updateBoundaryPoints(i, L, U, C, waveSpeedDetails.getR(), mu);
    }
}

MatrixXd calculateCFL(const Mesh2D& mesh, const Matrix2d& A)
{
    return A * (mesh.getTimeStep() / mesh.getPositionStep());
}


void updateInteriorPoint(int i, int j, UMatrix& U, const MatrixXd& mu)
{
    Vector2d previous(U.getPressure()(i, j - 1),
                      U.getVelocity()(i, j - 1));
    Vector2d next(U.getPressure()(i, j + 1),
                  U.getVelocity()(i, j + 1));
    MatrixXd interiorPoint = (previous + next) / 2 - (mu / 2) * (next - previous);
    U.setPoint(i + 1, j, interiorPoint);
}

void updateBoundaryPoints(int i, int L, UMatrix& U, CharacteristicMatrix& C, const MatrixXd& R, const MatrixXd& mu)
{
    int j = 0;
    Vector2d previous(U.getPressure()(i, j), U.getVelocity()(i, j));
    previous = R * previous;
    Vector2d next(U.getPressure()(i, j + 1), U.getVelocity()(i, j + 1));
    next = R * next;
    Vector2d current = mu * (previous - next) + previous;
    C.setZ(i, j, current(1));
    Vector2d newC(C.getW()(i, j), C.getZ()(i, j));
    MatrixXd newU = R.lu().solve(newC);
    U.setPoint(i + 1, j, newU);

    j = L;
    previous << U.getPressure()(i, j), U.getVelocity()(i, j);
    previous = R * previous;
    next << U.getPressure()(i, j - 1), U.getVelocity()(i, j - 1);
    next = R * next;
    current = mu * (previous - next) + previous;
    C.setW(i, 1, current(0));
    newC(0) = C.getW()(i, 1);
    newC(1) = C.getZ()(i, 1);
    newU = R.lu().solve(newC);
    U.setPoint(i + 1, j, newU);
}

