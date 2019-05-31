#include "modelTwoDsynapse.h"
#include "math.h"
#include <iostream.h>
#include<fstream.h>
/*
This is a 2D model for synapse when the rise time is much slower compared 
to the action potential width, as for example is the case in the situation
of julies model.
*/
void
modelTwoDsynapse (double t, double *x, double *dx, double *parameter, double
		*extra)
{
  double tau, s1=0, so;
  double tau1,s11;
  double Theta=0; 
  double s2=0;
  
  tau =parameter[3]-parameter[2]; 
  s1 = parameter[3]/tau;

tau1=8.0-parameter[2];
s11=8.0/tau1;

Theta=extra[1];
  so = 1. / 2 * (1 + tanh (200 * (extra[1] - .1)));
  dx[0] = (1 / tau) * (so - x[0]) / (s1 - so);
  dx[1]= (1 / tau1) * (so - x[1]) / (s11 - so);

/*
Theta=extra[1];
s1=1. / 2 * (1 + tanh (120 * (Theta - .1)));
s2=1. / 2 * (1 + tanh (120 * (-Theta + .1)));

dx[0]=(1-x[0])*s1/parameter[2]-x[0]*s2/parameter[3];
*/
}
