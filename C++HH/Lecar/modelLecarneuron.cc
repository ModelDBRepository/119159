
#include "modelLecarneuron.h"
#include "math.h"
#include "iostream.h"
void
modelLecarneuron (double t, double *x, double *dx, double *parameter,
		 double *extra)
{
  double gca,gk,vca,vk,gl,vl,iapp,v,v1,v2,v3,v4,phi,CM;
  double minf,winf,tinf;
 
 gca=parameter[0];
 vca=parameter[1];
 gk=parameter[2];
 vk=parameter[3];
 gl=parameter[4];
 vl=parameter[5];
 iapp=parameter[6];
v1=parameter[7];
v2=parameter[8];
v3=parameter[9];
v4=parameter[10];
phi=parameter[11];
CM=parameter[12];

 v=x[0];

minf=.5*(1+tanh((v-v1)/v2));
winf=.5*(1+tanh((v-v3)/v4));
tinf=1./cosh((v-v3)/(2*v4));

dx[0] =(-gca*minf*(v-vca)-gk*x[1]*(v-vk)-gl*(v-vl)+iapp+extra[1])/CM;
dx[1]=phi*(winf-x[1])/tinf;

  

}
