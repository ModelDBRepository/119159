#include "modelLorenz.h"
#include "math.h"
#include "iostream.h"
void
modelLorenz (double t, double *x, double *dx, double *parameter,
		 double *extra)
{
  double v,rho,beta,sigma;
  
    rho=parameter[0];
	beta=parameter[1];
	sigma=parameter[2];
	
dx[0]=sigma*(x[1]-x[0]);
dx[1]=x[0]*(rho-x[2])-x[1];
dx[2]=x[0]*x[1]-beta*x[2];
}
