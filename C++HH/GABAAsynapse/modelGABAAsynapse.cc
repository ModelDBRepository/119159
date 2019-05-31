#include "modelGABAAsynapse.h"
#include "math.h"
#include <iostream.h>

void modelGABAAsynapse(double t,double *x,double *dx,double *parameter,double
*extra)
{
double T,tau,so,s1;
if (extra[0]==0)
{
T=1./(1+exp(-(extra[1]-parameter[4])/parameter[5]));
dx[0]=parameter[2]*T*(1-x[0])-parameter[3]*x[0];
}
else if (extra[0]==1.)
{

//tau=parameter[2]/(parameter[3]*(parameter[3]+parameter[2]));
//s1=(parameter[2]+parameter[3])/parameter[2];

tau=parameter[2];
s1=parameter[3];
//cout<<tau<<" "<<s1<<endl;
so=0.5 * (1 + tanh (120 * (extra[1] - .1)));
dx[0]=(1/tau)*(so-x[0])/(s1-so);
//cout<<tau<<" "<<" "<<s1<<" "<<so<<" "<<x[0]<<endl;

}
}


//the 'extra[1]' parameter signifies the the ode of the object needs input from other
//object at everytime tims step. 
