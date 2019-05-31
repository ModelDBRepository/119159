#include "modelRose.h"
#include "math.h"
#include "iostream.h"
void
modelRose (double t, double *x, double *dx, double *parameter,
		 double *extra)
{
  double v,A,B,C,R,S,I;
  
    A=parameter[0];
	B=parameter[1];
	C=parameter[2];
	R=parameter[3];
	S=parameter[4];
	I=parameter[5];

 v=x[0];

dx[0] =x[1]+A*x[0]*x[0]-pow(x[0],3.0)-x[2]+I+extra[1];
dx[1]=1-B*x[0]*x[0]-x[1];
dx[2]=R*(-x[2]+S*(x[0]-C));

  

}
