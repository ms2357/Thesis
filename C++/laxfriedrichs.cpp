#include "laxfriedrichs.h"
/*Implements the Lax Friedrichs method for solving the interior nodes, while the boundaries for each edge are
handled using the junction conditions via the functions for calculating the known/unk characteristic variables
then solving the system of equations governed by continuity of pressre and conservation of momentum*/
/*LxF Method*/
void laxfriedrichs(const Mesh2D& mesh, const WaveSpeed& waveSpeedDetails, UMatrix& UEdge1, UMatrix& UEdge2,
                   MatrixXd CharacteristicMatrixEdge1W, MatrixXd CharacteristicMatrixEdge1Z,
                   MatrixXd CharacteristicMatrixEdge2W, MatrixXd CharacteristicMatrixEdge2Z)
{
    MatrixXd mu = calculateCFL(mesh, waveSpeedDetails.getA());
    int N = mesh.getTimeMesh().rows() - 1;
    int L = mesh.getPositionMesh().rows()- 1;
    for (int i = 0; i < N; i++) {
        for (int j = 1; j < L; j++) {
            updateInteriorPoint(i, j, UEdge1, mu);
            updateInteriorPoint(i, j, UEdge2, mu);

        }
        calculateKnownCharacteristics(i, L, UEdge1, CharacteristicMatrixEdge1W, CharacteristicMatrixEdge1Z,
                                 waveSpeedDetails.getR(), waveSpeedDetails.getc0w(), waveSpeedDetails.getc0z(),
                                 mesh.getTimeStep(), mesh.getPositionStep());
        calculateKnownCharacteristics(i, L, UEdge2, CharacteristicMatrixEdge2W, CharacteristicMatrixEdge2Z,
                                 waveSpeedDetails.getR(), waveSpeedDetails.getc0w(), waveSpeedDetails.getc0z(),
                                 mesh.getTimeStep(), mesh.getPositionStep());
        calculateUnkownCharacteristics(i, UEdge1, CharacteristicMatrixEdge1W, CharacteristicMatrixEdge1Z,
                                         CharacteristicMatrixEdge2W, CharacteristicMatrixEdge2Z,
                                        waveSpeedDetails.getZ(), waveSpeedDetails.getZcoefficientMatrix(),
                                        waveSpeedDetails.getRI());
        calculateUnkownCharacteristics(i, UEdge2, CharacteristicMatrixEdge1W, CharacteristicMatrixEdge1Z,
                                         CharacteristicMatrixEdge2W, CharacteristicMatrixEdge2Z,
                                        waveSpeedDetails.getZ(), waveSpeedDetails.getZcoefficientMatrix(),
                                        waveSpeedDetails.getRI());

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

void calculateKnownCharacteristics(int i, int L, UMatrix &U, MatrixXd& CharacteristicMatrixW, MatrixXd &CharacteristicMatrixZ,
                              const MatrixXd &R, double c0w, double c0z, double TimeStep, double PositionStep)
{
    //calc w char using pres/vel at vertex a
    int a = 0;
    Vector2d previousw(U.getPressure()(i, a),
                      U.getVelocity()(i, a));
    Vector2d previousW = R * previousw;
    Vector2d nextw(U.getPressure()(i, a + 1),
                  U.getVelocity()(i, a + 1));
    Vector2d nextW = R * nextw;
    Vector2d currentW = c0w * (TimeStep / PositionStep) * (previousW - nextW) + previousW;
    //update w char value at vertex a
    CharacteristicMatrixW(i + 1, a) = currentW(0);

    //calc  z char using pres/vel at vertex b
    int b = L;
    Vector2d previousz(U.getPressure()(i, b - 1),
                      U.getVelocity()(i, b - 1));
    Vector2d previousZ = R * previousz;
    Vector2d nextz(U.getPressure()(i, b),
                  U.getVelocity()(i, b));
    Vector2d nextZ = R * nextz;
    Vector2d currentZ = c0z * (TimeStep / PositionStep) * (previousZ - nextZ) + nextZ;
    //update z char value at vertex b
    CharacteristicMatrixZ(i + 1, 1) = currentZ(1);

}

void calculateUnkownCharacteristics(int i, UMatrix & U, MatrixXd& CharacteristicMatrixW1, MatrixXd& CharacteristicMatrixZ1,
                                    MatrixXd& CharacteristicMatrixW2, MatrixXd& CharacteristicMatrixZ2,
                                    const MatrixXd& Z, const MatrixXd& ZcoefficientMatrix, const MatrixXd& RI)
{
    //calc unk chars at vertex using juction conditions. Solve system Ax=b using interpolated chars at vertex a
    Vector2d knowncharacteristics(CharacteristicMatrixW1(i + 1, 0), CharacteristicMatrixZ2(i +1, 0));
    Vector2d KnownCharacteristics = Z * knowncharacteristics;
    Vector2d UnknownCharacteristics = ZcoefficientMatrix * KnownCharacteristics;

    //set new chars at vertex
    CharacteristicMatrixW2(i + 1, 0) = UnknownCharacteristics(0);
    CharacteristicMatrixZ1(i + 1, 0) = UnknownCharacteristics(1);
    //calc press/vel using chars at vertex a
    Vector2d newcharacteristics(CharacteristicMatrixW1(i + 1, 0), CharacteristicMatrixZ1(i + 1, 0));
    Vector2d NewCharacteristics = RI * newcharacteristics;
    U.setPoint(i + 1,0, NewCharacteristics);



}
