#include "modelIh.h"
#include "math.h"

void modelIh(double t,double *x,double *dx,double *parameter,double
*extra)
{
double v,hinf,tinf;
v=extra[1];
hinf=1./(1+exp(-(v+35.)/10.0));
tinf=3900.0/(exp(-7.68-.086*v)+exp(5.04+.0701*v));
dx[0]=(hinf-x[0])/tinf;


}

/*m1=1./(1+exp(-(v+60)/8.5));
m2=1./(1+exp(-(v+36)/20));
tma=(1./(exp((v+35.82)/19.69)+exp(-(v+79.69)/12.7))+.37);
h1=1./(1+exp((v+78)/6));
tadef=1./(exp((v+46.05)/5.)+exp(-(v+238.4)/37.45));
if(v<-63)
ta1=tadef;
else
ta1=19;
if(v<-73)
ta2=tadef;
else
ta2=60;

dx[0]=(m1-x[0])/tma;
dx[1]=(m2-x[1])/tma;
dx[2]=(h1-x[2])/ta1;
dx[3]=(h1-x[3])/ta2;*/
