
#include "modelQuadraticNeuron.h"
#include "math.h"
#include "iostream.h"
void
modelQuadraticNeuron (double t, double *x, double *dx, double *parameter,
		 double *extra)
{
  double v,  b, vreset,vpeak,k,c,vrest;
  double type;
  b = parameter[0];
  vreset = parameter[1];
  vpeak = parameter[2];
 k=parameter[3];
 vrest=parameter[4];
 c=parameter[5];
x[1]=0; 
dx[0] =k*pow((x[0]-parameter[4]),2.0)+parameter[5]*b+extra[1];

  if (x[0]>vpeak)
  {
  x[0]=vreset;
  x[1]=1;
  }

}
