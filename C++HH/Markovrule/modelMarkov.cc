#include "modelMarkov.h"
#include "math.h"
#include <iostream.h>
void modelMarkov(double t,double *x,double *dx,double *parameter,double
*extra)
{
double CA,fp,fd,L,M,epsilon_p,epsilon_d,a,b,n,bp,bd;
double gain,loss,alp;

L=parameter[0];
M=parameter[1];
a=parameter[2];
b=parameter[3];
epsilon_p=parameter[4];
epsilon_d=parameter[5];
n=parameter[6];
bp=parameter[7];
bd=parameter[8];
CA=extra[1];
if (CA<0) CA=0;

fp=a*pow(CA,L)/(pow(CA,L)+pow(epsilon_p,L));
fd=b*pow(CA,M)/(pow(CA,M)+pow(epsilon_d,M));
alp=parameter[11];

gain=x[0]*pow(x[1],4.);
loss=x[1]*pow(x[0],4.);


dx[0]=fp*(1-x[0])-(1./bp)*x[0];
dx[1]=fd*(1-x[1])-(1./bd)*x[1];
dx[2]=-1*gain*x[2]+loss*x[3];
dx[3]=-(loss+gain)*x[3]+gain*x[2] +alp*gain*x[4];
dx[4]=gain*x[3]-alp*gain*x[4];
}

