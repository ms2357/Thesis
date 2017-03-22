#include "laxfriedrichs.h"

/*Implements the Lax Friedrichs method for solving the interior nodes, while the boundaries for each edge are
handled using the junction conditions via the functions for calculating the known/unk characteristic variables
then solving the system of equations governed by continuity of pressure and conservation of momentum*/
/*LxF Method*/
void laxfriedrichs(const Mesh2D& mesh, const WaveSpeed& waveSpeedDetails, UMatrix& UEdge1, UMatrix& UEdge2,
                   CharacteristicMatrix& CharacteristicsEdge1, CharacteristicMatrix& CharacteristicsEdge2)
{
    MatrixXd mu = calculateCFL(mesh, waveSpeedDetails.getA());
    int N = mesh.getTimeMesh().rows() - 1;
    int L = mesh.getPositionMesh().rows()- 1;
    for (int i = 0; i < N; i++) {
        for (int j = 1; j < L; j++) {
            updateInteriorPoint(i, j, UEdge1, mu);
            updateInteriorPoint(i, j, UEdge2, mu);

        }
        calcKnownCharacteristics(i, L, UEdge1, CharacteristicsEdge1, waveSpeedDetails.getR(), waveSpeedDetails.getc0w(), waveSpeedDetails.getc0z(),
                                 mesh.getTimeStep(), mesh.getPositionStep());
        calcKnownCharacteristics(i, L, UEdge2, CharacteristicsEdge2, waveSpeedDetails.getR(), waveSpeedDetails.getc0w(), waveSpeedDetails.getc0z(),
                                 mesh.getTimeStep(), mesh.getPositionStep());
        JunctionConditions(i, UEdge1, CharacteristicsEdge1, CharacteristicsEdge2, waveSpeedDetails.getZ(),
                           waveSpeedDetails.getZcoefficientMatrix(), waveSpeedDetails.getRI());
        JunctionConditions(i, UEdge2, CharacteristicsEdge1, CharacteristicsEdge2, waveSpeedDetails.getZ(),
                           waveSpeedDetails.getZcoefficientMatrix(), waveSpeedDetails.getRI());
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

void calculateKnownCharacteristics(int i, int L, UMatrix &U, CharacteristicMatrix& CharacteristicsEdge1,
                                   const MatrixXd &R, double c0w, double c0z, double TimeStep, double PositionStep)
{
    //calc w char using pres/vel at right vertex, a
    int a = 0;
    Vector2d previousw(U.getPressure()(i, a),
                      U.getVelocity()(i, a));
    Vector2d previousW = R * previousw;
    Vector2d nextw(U.getPressure()(i, a + 1),
                  U.getVelocity()(i, a + 1));
    Vector2d nextW = R * nextw;
    Vector2d currentW = c0w * (TimeStep / PositionStep) * (previousW - nextW) + previousW;
    //update w char value at vertex a
    CharacteristicsEdge1.setW(i + 1, a, currentW(0));

    //calc  z char using pres/vel at left vertex, b
    int b = L;
    Vector2d previousz(U.getPressure()(i, b - 1),
                      U.getVelocity()(i, b - 1));
    Vector2d previousZ = R * previousz;
    Vector2d nextz(U.getPressure()(i, b),
                  U.getVelocity()(i, b));
    Vector2d nextZ = R * nextz;
    Vector2d currentZ = c0z * (TimeStep / PositionStep) * (previousZ - nextZ) + nextZ;
    //update z char value at vertex b
    CharacteristicsEdge1.setZ(i + 1, 1, currentZ(1));

}

void JunctionConditions(int i, UMatrix& U, CharacteristicMatrix& CharacteristicsEdge1, CharacteristicMatrix& CharacteristicsEdge2,
                        const MatrixXd& Z, const MatrixXd& ZcoefficientMatrix, const MatrixXd& RI)
{
    /*calc unk chars at vertex using juction conditions. Solve system Ax=b using interpolated chars at right hand
     * vertex, a*/

    Vector2d knowncharacteristicsA(CharacteristicsEdge1.getW(i + 1, 0), CharacteristicsEdge2.getZ(i +1, 0));
    Vector2d KnownCharacteristicsA = Z * knowncharacteristicsA;
    Vector2d UnknownCharacteristicsA = ZcoefficientMatrix * KnownCharacteristicsA;

    //set new chars at vertex
    CharacteristicsEdge2.setW(i + 1, 0, UnknownCharacteristicsA(0));
    CharacteristicsEdge1.setZ(i + 1, 0, UnknownCharacteristicsA(1));

    //calc pressure/velocity using chars at vertex a
    Vector2d newcharacteristicsA(CharacteristicsEdge1.getW(i + 1, 0), CharacteristicsEdge1.getZ(i + 1, 0));
    Vector2d NewCharacteristicsA = RI * newcharacteristicsA;
    U.setPoint(i + 1, 0, NewCharacteristicsA);


    /*calc unk chars at vertex using juction conditions. Solve system Ax=b using interpolated chars at lefy hand
     * vertex, b*/

    Vector2d knowncharacteristicsB(CharacteristicsEdge2.getW(i + 1, 1), CharacteristicsEdge1.getZ(i +1, 1));
    Vector2d KnownCharacteristicsB = Z * knowncharacteristicsB;
    Vector2d UnknownCharacteristicsB = ZcoefficientMatrix * KnownCharacteristicsB;

    //set new chars at vertex
    CharacteristicsEdge1.setW(i + 1, 1, UnknownCharacteristicsA(0));
    CharacteristicsEdge2.setZ(i + 1, 1, UnknownCharacteristicsA(1));

    //calc pressure/velocity using chars at vertex a
    Vector2d newcharacteristicsB(CharacteristicsEdge1.getW(i + 1, 1), CharacteristicsEdge1.getZ(i + 1, 1));
    Vector2d NewCharacteristicsB = RI * newcharacteristicsB;
    U.setPoint(i + 1, 0, NewCharacteristicsB);



}
