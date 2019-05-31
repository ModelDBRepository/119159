/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-02-14
  
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
   This an abstract class which is a wrapper for all types of input
   neurons (multifire_inputneuron, noisy input neuron etc)
--------------------------------------------------------------------------*/


#ifndef CN_INPUTNEURON_H
#define CN_INPUTNEURON_H

#include "CN_neuron.h"

class inputneuron: public neuron
{
 public:
  int fno;
  double *tb, *te;
  double theE, t_last;
  
  inputneuron(int, int, int, tnvector<int>, double *, int);
  ~inputneuron();
  virtual void init(double *, double *);
  virtual void set_input(int, double *, double *, double);
  virtual double E(double *)= 0;
  virtual void derivative(double *, double *) { }
};


#endif





