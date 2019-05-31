#include "modelAMPAsynapse.h"
#include "math.h"


void modelAMPAsynapse(double t,double *x,double *dx,double *parameter,double
*extra)
{
double T,tau,s1,so;
if (extra[0]==0) //extra[0] characterises the step or voltage transfer, of
                  //synapse, in the sense of distinguish from henry and ampa
		  //synapse
{
T=1./(1+exp(-(extra[1]-parameter[4])/parameter[5]));
dx[0]=parameter[2]*T*(1-x[0])-parameter[3]*x[0];
}
else
{
tau=(parameter[2]-parameter[3])/(parameter[2]*parameter[3]);
s1=parameter[2]/(parameter[2]-parameter[3]);
so=1. / 2 * (1 + tanh (120 * (extra[1] - .1)));
dx[0]=(1/tau)*(so-x[0])/(s1-so);
}
}


//the 'extra[1]' parameter signifies the the ode of the object needs input from other
//object at everytime tims step. 
