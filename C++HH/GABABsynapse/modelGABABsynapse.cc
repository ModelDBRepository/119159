#include "modelGABABsynapse.h"
#include "math.h"


void modelGABABsynapse(double t,double *x,double *dx,double *parameter,double
*extra)
{
double T;
T=1./(1+exp(-(extra[1]-parameter[4])/parameter[5]));
dx[0]=parameter[2]*T*(1-x[0])-parameter[3]*x[0];
dx[1]=parameter[6]*x[0]-parameter[7]*x[1];
}


//the 'extra' parameter signifies the the ode of the object needs input from other
//object at everytime tims step. 
