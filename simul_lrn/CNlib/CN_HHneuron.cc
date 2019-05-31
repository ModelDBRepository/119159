/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_HHNEURON_CC
#define CN_HHNEURON_CC

#include "CN_neuron.cc"
#include "CN_HHneuron.h"

HHneuron::HHneuron(int inlabel, double *inp= HH_p):
  neuron(inlabel, HH_IVARNO, HHNEURON, inp, HH_PNO)
{
}

HHneuron::HHneuron(int inlabel, tnvector<int> inpos, double *inp= HH_p):
  neuron(inlabel, HH_IVARNO, HHNEURON, inpos, inp, HH_PNO)
{
}

inline double HHneuron::E(double *x)
{
  return x[idx];
}

void HHneuron::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den_it) {
    Isyn+= den_it->c_value()->Isyn(x);
  }
  
  // differential eqn for E, the membrane potential
  dx[idx]= -(ipower(x[idx+1],3)*x[idx+2]*p[0]*(x[idx]-p[1])+
	       ipower(x[idx+3],4)*p[2]*(x[idx]-p[3])+p[4]*(x[idx]-p[5])
	       -Isyn)/p[6];

  // diferential eqn for m, the probability for Na channel activation
  _a= (3.5+0.1*x[idx]) / (1.0-exp(-3.5-0.1*x[idx]));
  _b= 4.0*exp(-(x[idx]+60.0)/18.0);
  dx[idx+1]= _a*(1.0-x[idx+1])-_b*x[idx+1];

  // differential eqn for h, the probability for Na channel inactivation
  _a= 0.07*exp(-x[idx]/20.0-3.0);   
  _b= 1.0 / (exp(-3.0-0.1*x[idx])+1.0);
  dx[idx+2]= _a*(1.0-x[idx+2])-_b*x[idx+2];

  // differential eqn for n, the probability for K channel activation
  _a= (-0.5-0.01*x[idx]) / (exp(-5.0-0.1*x[idx])-1.0); 
  _b= 0.125*exp(-(x[idx]+60.0)/80.0);
  dx[idx+3]= _a*(1.0-x[idx+3])-_b*x[idx+3];
}

#endif
