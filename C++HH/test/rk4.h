/*****************************************************************/
/*Sachin S Talathi
Date 03/27/04 */
/*4th order fixed time step Runge Kutta Solver for ODE */
/******************************************************************/
#ifndef _rk4_H_
#define _rk4_H_
#define epsilon .01
void step(double timestep,double *x,double *dx,double *parameter,double
*extra,double *xnew,int NDIM,double *tpoint,void (*)(double, double *,double *,double *,double *));
#endif
