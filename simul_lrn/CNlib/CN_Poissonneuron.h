/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_POISSONNEURON_H
#define CN_POISSONNEURON_H

#include "CN_neuron.h"
#include "randomGen.h"

#define POI_IVARNO 0
#define POI_PNO 5

randomGen R;

double stdPOI_p[POI_PNO]= {
  2,                       // 0 - spike time of mf_poisson inputneuron
  10.0,                    // 1 - refractory period + spike time
  -60.0,                   // 2 - input neuron resting potential
  50.0,                    // 3 - input neuron potential when firing
  0.1                      // 4 - firing rate Lambda [1/ms]=[10^3 Hz]
};

double *POI_p= stdPOI_p;

class Poissonneuron: public neuron
{
 private:
  double fire_t;
  int firing;
  int refract;
  double *myx, *myxn;
  
 public:
  Poissonneuron(int, tnvector<int>, double *);
  ~Poissonneuron();
  virtual double E(double *);
  virtual void validate_E(double *, double);
  virtual void step();
  virtual void derivative(double *, double *) { }
  virtual void init(double *, double *);
};

#endif





