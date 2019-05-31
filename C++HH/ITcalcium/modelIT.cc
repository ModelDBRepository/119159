#include "modelIT.h"
#include "math.h"

void modelIT(double t,double *x,double *dx,double *parameter,double *extra)
{
double v;
double mco,tmco,hco,thco;
v=extra[1];
mco= 1./(1+exp(-(52+v)/6.2));
tmco=.204+.333/(exp(-(131+v)/16.7)+exp((15+v)/18.2));
hco=1./(1+exp((72+v)/4));
if(v<-81)
thco=.333*exp((v+466)/66.6);
else
thco=9.32+.333*exp(-(v+21)/10.5);

dx[0]=(mco-x[0])/tmco;
dx[1]=(hco-x[1])/thco;
}

