#include "stddef.h"
#include "stdarg.h"
#include "math.h"
#include "CVODEmex25.h"
#include "MEX.h"
#include "mexsplineaddon.h"
#include "mexmathaddon.h"
#include "kineticformulas.h"

double time;

void model(double time_local, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector)
{
    double X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16;
    double r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,k1,k2,k3,k4;
    double k5,k6,k7,k8;

    time = time_local;

    X1 = stateVector[0];
    X2 = stateVector[1];
    X3 = stateVector[2];
    X4 = stateVector[3];
    X5 = stateVector[4];
    X6 = stateVector[5];
    X7 = stateVector[6];
    X8 = stateVector[7];
    X9 = stateVector[8];
    X10 = stateVector[9];
    X11 = stateVector[10];
    X12 = stateVector[11];
    X13 = stateVector[12];
    X14 = stateVector[13];
    X15 = stateVector[14];
    X16 = stateVector[15];
    r1 = paramdataPtr->parametervector[0]; /* 0.0125 */
    r2 = paramdataPtr->parametervector[1]; /* 0.025 */
    r3 = paramdataPtr->parametervector[2]; /* 0.0375 */
    r4 = paramdataPtr->parametervector[3]; /* 0.05 */
    r5 = paramdataPtr->parametervector[4]; /* 0.0625 */
    r6 = paramdataPtr->parametervector[5]; /* 0.075 */
    r7 = paramdataPtr->parametervector[6]; /* 0.0875 */
    r8 = paramdataPtr->parametervector[7]; /* 0.1 */
    r9 = paramdataPtr->parametervector[8]; /* 0.1125 */
    r10 = paramdataPtr->parametervector[9]; /* 0.125 */
    r11 = paramdataPtr->parametervector[10]; /* 0.1375 */
    r12 = paramdataPtr->parametervector[11]; /* 0.15 */
    r13 = paramdataPtr->parametervector[12]; /* 0.1625 */
    r14 = paramdataPtr->parametervector[13]; /* 0.175 */
    r15 = paramdataPtr->parametervector[14]; /* 0.1875 */
    r16 = paramdataPtr->parametervector[15]; /* 0.2 */
    k1 = paramdataPtr->parametervector[16]; /* 0.02 */
    k2 = paramdataPtr->parametervector[17]; /* 0.02 */
    k3 = paramdataPtr->parametervector[18]; /* 0.02 */
    k4 = paramdataPtr->parametervector[19]; /* 0.02 */
    k5 = paramdataPtr->parametervector[20]; /* 0.02 */
    k6 = paramdataPtr->parametervector[21]; /* 0.02 */
    k7 = paramdataPtr->parametervector[22]; /* 0.02 */
    k8 = paramdataPtr->parametervector[23]; /* 0.02 */
    if (DOflag == DOFLAG_DDT) {
    	DDTvector[0] = X1*(r1-k1*X2);
    	DDTvector[1] = X2*(r2-k2*X3);
    	DDTvector[2] = X3*(r3-k3*X4);
    	DDTvector[3] = X4*(r4-k4*X5);
    	DDTvector[4] = X5*(r5-k5*X6);
    	DDTvector[5] = X6*(r6-k6*X7);
    	DDTvector[6] = X7*(r7-k7*X8);
    	DDTvector[7] = X8*(r8-k8*X9);
    	DDTvector[8] = X9*(r9-k1*X10);
    	DDTvector[9] = X10*(r10-k2*X11);
    	DDTvector[10] = X11*(r11-k3*X12);
    	DDTvector[11] = X12*(r12-k4*X13);
    	DDTvector[12] = X13*(r13-k5*X14);
    	DDTvector[13] = X14*(r14-k6*X15);
    	DDTvector[14] = X15*(r15-k7*X16);
    	DDTvector[15] = X16*(r16-k8*X1);
    } else if (DOflag == DOFLAG_VARREAC) {
    } else if (DOflag == DOFLAG_EVENTS) {
    } else if (DOflag == DOFLAG_EVENTASSIGN) {
    }
}


/* Function for initial condition calculation */
void calc_ic_model(double *icVector, ParamData *paramdataPtr)
{
    double X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16;
    double r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,k1,k2,k3,k4;
    double k5,k6,k7,k8;
    r1 = paramdataPtr->parametervector[0]; /* 0.0125 */
    r2 = paramdataPtr->parametervector[1]; /* 0.025 */
    r3 = paramdataPtr->parametervector[2]; /* 0.0375 */
    r4 = paramdataPtr->parametervector[3]; /* 0.05 */
    r5 = paramdataPtr->parametervector[4]; /* 0.0625 */
    r6 = paramdataPtr->parametervector[5]; /* 0.075 */
    r7 = paramdataPtr->parametervector[6]; /* 0.0875 */
    r8 = paramdataPtr->parametervector[7]; /* 0.1 */
    r9 = paramdataPtr->parametervector[8]; /* 0.1125 */
    r10 = paramdataPtr->parametervector[9]; /* 0.125 */
    r11 = paramdataPtr->parametervector[10]; /* 0.1375 */
    r12 = paramdataPtr->parametervector[11]; /* 0.15 */
    r13 = paramdataPtr->parametervector[12]; /* 0.1625 */
    r14 = paramdataPtr->parametervector[13]; /* 0.175 */
    r15 = paramdataPtr->parametervector[14]; /* 0.1875 */
    r16 = paramdataPtr->parametervector[15]; /* 0.2 */
    k1 = paramdataPtr->parametervector[16]; /* 0.02 */
    k2 = paramdataPtr->parametervector[17]; /* 0.02 */
    k3 = paramdataPtr->parametervector[18]; /* 0.02 */
    k4 = paramdataPtr->parametervector[19]; /* 0.02 */
    k5 = paramdataPtr->parametervector[20]; /* 0.02 */
    k6 = paramdataPtr->parametervector[21]; /* 0.02 */
    k7 = paramdataPtr->parametervector[22]; /* 0.02 */
    k8 = paramdataPtr->parametervector[23]; /* 0.02 */
    X1 = 1.0;
    X2 = 1.0;
    X3 = 1.0;
    X4 = 1.0;
    X5 = 1.0;
    X6 = 1.0;
    X7 = 1.0;
    X8 = 1.0;
    X9 = 1.0;
    X10 = 1.0;
    X11 = 1.0;
    X12 = 1.0;
    X13 = 1.0;
    X14 = 1.0;
    X15 = 1.0;
    X16 = 1.0;
    icVector[0] = X1;
    icVector[1] = X2;
    icVector[2] = X3;
    icVector[3] = X4;
    icVector[4] = X5;
    icVector[5] = X6;
    icVector[6] = X7;
    icVector[7] = X8;
    icVector[8] = X9;
    icVector[9] = X10;
    icVector[10] = X11;
    icVector[11] = X12;
    icVector[12] = X13;
    icVector[13] = X14;
    icVector[14] = X15;
    icVector[15] = X16;
}

