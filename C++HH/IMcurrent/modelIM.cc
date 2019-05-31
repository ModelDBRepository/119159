#include "modelIM.h"
#include "math.h"

void modelIM(double t,double *x,double *dx,double *parameter,double
*extra)
{
double v,winf,tauinf;
v=extra[1];
winf=1.0/(1.0+exp(-(v+35)/10.0));
tauinf=100/(3.3*exp((v+35)/20.0)+exp(-(v+35)/20.0));
dx[0]=(winf-x[0])/tauinf;
/*
double alpham,betam;
alpham=.016 / exp (-(v + 52.7) / 23);
betam=.016 / exp ((v + 52.7) / 18.8); // Old model for IM current
dx[0]=alpham*(1-x[0])-betam*x[0];
*/
}

