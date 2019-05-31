/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-01-25
  
--------------------------------------------------------------------------*/


#ifndef SYNAPSEASTRID_H
#define SYNAPSEASTRID_H

#include "CN_synapse.h"
#include "CN_neuron.h"

#define SYNASPNO 5

#define SYNASIVARNO 1


double SYNAS_INIVARS[SYNASIVARNO]= {
  0.0                    // 0 - S: internal variable 
};


class synapseAstrid: public synapse
{
 private:
  double sinf, tau;
 public:
  synapseAstrid(neuron *, neuron *, double, double, double, double, double,
		int, int, int);
  synapseAstrid(neuron *, neuron *, double, double, double, double, double);
  synapseAstrid(neuron *, neuron *, double *);
  virtual ~synapseAstrid();
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn(double *);
  virtual void derivative(double *, double *);
};

#endif
