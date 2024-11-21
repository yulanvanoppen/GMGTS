#include "mex.h"

const int NRSTATES = 16;
const int NRPARAMETERS = 24;
const int NRVARIABLES = 0;
const int NRREACTIONS = 0;
const int NREVENTS = 0;

const int hasOnlyNumericICs = 1;
double defaultICs_num[16] = {
	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
char *defaultICs_nonnum[1];

double defaultParam[24] = {
	0.0125,0.025,0.0375,0.05,0.0625,0.075,0.0875,0.1,0.1125,0.125,0.1375,0.15,0.1625,0.175,0.1875,0.2,0.02,0.02,0.02,0.02,
	0.02,0.02,0.02,0.02};
char *stateNames[16] = {
	"X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16"};
char *parameterNames[24] = {
	"r1","r2","r3","r4","r5","r6","r7","r8","r9","r10","r11","r12","r13","r14","r15","r16","k1","k2","k3","k4",
	"k5","k6","k7","k8"};
char *variableNames[1];
char *variableFormulas[1];
char *reactionNames[1];
char *eventNames[1];

void model(double time, double *stateVector, double *DDTvector, ParamData *paramdataPtr, int DOflag, double *variableVector, double *reactionVector, double *gout, int *eventVector);
void calc_ic_model(double *icVector, ParamData *paramdataPtr);

void CVODEmex25(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    CVODEmex25(nlhs, plhs, nrhs, prhs);
}
