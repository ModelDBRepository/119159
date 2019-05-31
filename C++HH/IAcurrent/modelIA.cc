#include "modelIA.h"
#include "math.h"

void modelIA(double t,double *x,double *dx,double *parameter,double
*extra)
{
double v,alphaa,alphab,betaa,betab;
v=extra[1];
alphaa=-.05 * (v + 20) / (exp (-(v + 20) / 15) - 1);
betaa=.1 * (v + 10) / (exp ((v + 10) / 8) - 1);
alphab=.00015 / exp ((v + 18) / 15);
betab=.06 / (exp (-(v + 73) / 12) + 1);

dx[0]=alphaa*(1-x[0])-betaa*x[0];
dx[1]=alphab*(1-x[1])-betab*x[1];


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
