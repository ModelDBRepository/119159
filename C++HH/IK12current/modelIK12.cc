#include "modelIK12.h"
#include "math.h"

void modelIK12(double t,double *x,double *dx,double *parameter,double
*extra)
{
double m1,tm1,h1,th1,th2,v;
v=extra[1];
m1=1./pow((1+exp(-(v+43)/7)),4.);
tm1=(1./(exp((v-80.98)/25.64)+exp(-(v+132)/17.953))+9.9);
h1=1./(1+exp((v+58)/10.6));
th1=(1./(exp((v-1329)/200)+exp(-(v+129.7)/7.143))+120);
if(v<-70)
th2=8930;
else
th2=th1;
dx[0]=(m1-x[0])/tm1;
dx[1]=(h1-x[1])/th1;
dx[2]=(h1-x[2])/th2;
}
