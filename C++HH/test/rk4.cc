/*****************************************************************/
/*Sachin S Talathi
Date 03/27/04 */
/*4th order fixed time step Runge Kutta Solver for ODE */
/******************************************************************/

//#include "define.h"
#include "rk4.h"
#include "math.h"
#define Max 50
void step(double timestep,double *x,double *dx,double *parameter,double
*extra,double *xnew,int NDIM,double *tpoint,void (*model)(double, double *,double *,double *,double *))
{
double half_tau,t_half,t_full,xtemp[NDIM];
double F1[NDIM],F2[NDIM],F3[NDIM],F4[NDIM];
half_tau=.5*timestep;
model(*tpoint,x,F1,parameter,extra);
t_half=*tpoint+half_tau;
for(int i=0;i<NDIM;i++)
xtemp[i]=x[i]+F1[i]*half_tau;
model(t_half,xtemp,F2,parameter,extra);
for(int i=0;i<NDIM;i++)
xtemp[i]=x[i]+F2[i]*half_tau;
model(t_half,xtemp,F3,parameter,extra);
t_full=*tpoint+timestep;
for(int i=0;i<NDIM;i++)
xtemp[i]=x[i]+F3[i]*timestep;
model(t_full,xtemp,F4,parameter,extra);
for(int i=0;i<NDIM;i++)
xnew[i]=x[i]+timestep/6.0*(F1[i]+F4[i]+2*(F2[i]+F3[i]));

}
