#include "modelPD.h"
#include "math.h"
#include <iostream.h>
void modelPD(double t,double *x,double *dx,double *parameter,double
*extra)
{
double Vpre,Vpost,fp,fd,ap,bp,ad,bd,n,F,G,a,b;
double *tpre;
double *tpost;
tpre=new double[1];
tpost=new double[1];

double Scale;

n=parameter[0];
ap=parameter[1];
bp=parameter[2];
ad=parameter[3];
bd=parameter[4];
a=parameter[5];
b=parameter[6];


Vpre=extra[1];
Vpost=extra[2];

fp=.5*(1+tanh(120*(Vpre-.1)));
fd=.5*(1+tanh(120*(Vpost-.1)));


dx[0]=ap*fp*(1-x[0])-bp*x[0];
dx[1]=ad*fd*(1-x[1])-bd*x[1];

F=x[0]*pow(x[1],n);
G=pow(x[0],n)*x[1];

dx[2]=-0.1*F*x[2]+1.5*G*(1-x[2]-x[3]);
dx[3]=-a*F*x[3]+b*F*(1-x[2]-x[3]);

dx[4]=F-G;

}

//Modified the PD code to develop a continuous learning rule that is bounded by the contrainsts of the model based on my PRE paper and Henrys PNAS paper. 05/20/07