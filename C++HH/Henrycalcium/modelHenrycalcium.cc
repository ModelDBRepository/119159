#include "modelHenrycalcium.h"
#include "math.h"

void modelHenrycalcium(double t,double *x,double *dx,double *parameter,double
*extra)
{
double So;
 So = 1. / 2 * (1 + tanh (120 * (extra[1] - .1)));
dx[0]=(1./parameter[2])*(So-x[0])/(parameter[3]-So);
}

