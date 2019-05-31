/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
   absynapse with inhibitory synapse learning rule of the Entorhnal cortex
   ... no saturation mechanism...
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSEECPLAST2_H
#define CN_ABSYNAPSEECPLAST2_H

#include "queue.h"
#include "CN_absynapse.h"

#define ABECPLAST2PNO 7

#define ABECPLAST2IVARNO 1

double ABECPLAST2_INIVARS[ABECPLAST2IVARNO]= {
  0.0                    // 1 - S: internal variable
};

class absynapseECplast2: public absynapse
{  
 public:
  int synapse_change;
  absynapseECplast2(neuron *, neuron *, double, double, double, double,
			     double, double, double); 
  virtual ~absynapseECplast2();
  virtual void update_gsyn(double *);
  virtual double STDP_func(double);
};

#endif
