/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_VALNEURON_CC
#define CN_VALNEURON_CC

#include "CN_neuron.cc"

Valneuron::Valneuron(int inlabel, double *the_p= Val_p):
  neuron(inlabel, Val_IVARNO, VALNEURON, the_p, Val_PNO)
{
}

Valneuron::Valneuron(int inlabel, tnvector<int> inpos, double *the_p= Val_p):
  neuron(inlabel, Val_IVARNO, VALNEURON, inpos, the_p, Val_PNO)
{
}

inline double Valneuron::E(double *x)
{
  assert(enabled);
  return x[idx];
}

void Valneuron::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall(den_it) {
    Isyn+= den_it->c_value()->Isyn(x);
  }
  
  // differential eqn for E, the membrane potential
  dx[idx]= -(ipower(x[idx+1],3)*x[idx+2]*p[0]*(x[idx]-p[1]) +
		      ipower(x[idx+3],4)*p[2]*(x[idx]-p[3])+
		      p[4]*(x[idx]-p[5])+p[6]*(x[idx]-p[7])-Isyn)/p[9];

  // diferential eqn for m, the probability for one Na channel activation
  // particle
  _a= 0.32*(13.0-x[idx]-p[8]) / (exp((13.0-x[idx]-p[8])/4.0)-1.0);
  _b= 0.28*(x[idx]+p[8]-40.0)/(exp((x[idx]+p[8]-40.0)/5.0)-1.0);
  dx[idx+1]= (_a*(1.0-x[idx+1])-_b*x[idx+1])/p[10];

  // differential eqn for h, the probability for the Na channel blocking
  // particle to be absent
  _a= 0.128*exp((17.0-x[idx]-p[8])/18.0);   
  _b= 4.0 / (exp((40-x[idx]-p[8])/5.0)+1.0);
  dx[idx+2]= (_a*(1.0-x[idx+2])-_b*x[idx+2])/p[10];

  // differential eqn for n, the probability for one K channel activation
  // particle
  _a= .032*(15.0-x[idx]-p[8]) / (exp((15.0-x[idx]-p[8])/5.0)-1.0); 
  _b= 0.5*exp((10.0-x[idx]-p[8])/40.0);
  dx[idx+3]= (_a*(1.0-x[idx+3])-_b*x[idx+3])/p[10];
}

#endif
