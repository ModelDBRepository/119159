/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-18
  
--------------------------------------------------------------------------*/

#ifndef CN_MULTIFIRE_INPUTNEURON_H
#define CN_MULTIFIRE_INPUTNEURON_H

#include "CN_inputneuron.h"

#define MF_I_IVARNO 0
#define MF_I_PNO 5

double stdINPUT_p[MF_I_PNO]= {
  2,                       // spike time of multifire inputneuron
  10.0,                    // refractory period + spike time
  -60.0,                   // input neuron resting potential
  50.0,                    // input neuron potential when firing
  10.0                     // period of the periodic input signal
};

double *INPUT_p= stdINPUT_p;

class multifire_inputneuron: public inputneuron
{
 private:
  double red_tx;
  double difft;
  
public:
  multifire_inputneuron(int, tnvector<int>, double *);
  ~multifire_inputneuron();
  virtual double E(double *);
  virtual void integrate(double *, double *) { }
};


#endif





