/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_RALLSYNAPSE_H
#define CN_RALLSYNAPSE_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define RPNO 4

#define RIVARNO 2

double RASYN_INIVARS[RIVARNO]= {
  0.0,                    // 0 - f: internal variable 
  0.0                     // 1 - g: internal variable
};

class Rallsynapse: public synapse
{
  
 public:
  Rallsynapse(neuron *, neuron *, double, double, double, double,
	      int, int, int);
  Rallsynapse(neuron *, neuron *, double, double, double, double);
  Rallsynapse(neuron *, neuron *, double *);
  virtual ~Rallsynapse();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void derivative(double *, double *);
};

#endif
