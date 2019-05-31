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
   ... with sigmoidal filter saturation mechanism...
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSEECPLAST1_H
#define CN_ABSYNAPSEECPLAST1_H

#include "queue.h"
#include "CN_absynapse.h"

#define ABECPLAST1PNO 10

#define ABECPLAST1IVARNO 1

double ABECPLAST1_INIVARS[ABECPLAST1IVARNO]= {
  0.0                    // 1 - S: internal variable
};

class absynapseECplast1: public absynapse
{  
 public:
  int synapse_change;
  absynapseECplast1(neuron *, neuron *, double, double, double, double,
			     double, double, double, double, double, double); 
  virtual ~absynapseECplast1();
  virtual double rgsyn();
  virtual void set_rgsyn(double);
  virtual double gsyn();
  virtual void set_gsyn(double);
  virtual double Isyn();
  virtual void update_gsyn(double *);
  virtual double STDP_func(double);
};

#endif
