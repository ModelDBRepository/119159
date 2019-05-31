/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_TIMENEURON_CC
#define CN_TIMENEURON_CC

#include "CN_neuron.cc"
#include "CN_TimeNeuron.h"

TimeNeuron::TimeNeuron():
  neuron(-1, TIME_IVARNO, TIMENEURON, NULL, TIME_PNO)
{
}

inline double TimeNeuron::E(double *x)
{
  return x[idx];
}

void TimeNeuron::derivative(double *x, double *dx)
{
  dx[idx]= 1.0;
}

#endif
