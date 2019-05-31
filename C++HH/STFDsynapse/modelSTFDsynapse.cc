#include "modelSTFDsynapse.h"
#include "math.h"
#include "iostream.h"

void modelSTFDsynapse(double t,double *x,double *dx,double *parameter,double
*extra)
{
double taurec,tauin,taufac,USE,so,T;
USE=parameter[1];
tauin=parameter[2];
taurec=parameter[3];
taufac=parameter[4];

if (extra[0]==1)
{
T=x[3]*(1-USE)+USE; //for facilitation
//cout<<"A"<<endl;
}
else
{
	T=USE; //model depression
}

so=1. / 2 * (1 + tanh (120 * (extra[1] - .1)));

dx[0]=x[2]/taurec-T*x[0]*so;
dx[1]=-x[1]/tauin+T*x[0]*so;
dx[2]=-x[2]/taurec+x[1]/tauin;
dx[3]=-x[3]/taufac+USE*(1-x[3])*so;



//cout<<x[0]<<" "<<T<<" "<<T*x[0]<<endl;
//cout<<-x[2]/taurec+x[1]/tauin<<" "<<-x[1]/tauin+T*x[0]*so<<" "<<x[2]/taurec-T*x[0]*so<<" "<<T<<endl;
}

//the 'extra[1]' parameter signifies the the ode of the object needs input from other
//object at everytime tims step. 
