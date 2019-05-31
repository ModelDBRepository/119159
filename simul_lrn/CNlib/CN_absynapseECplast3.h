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
   ... no saturation mechanism ... but delayed action of plasticity ...
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSEECPLAST3_H
#define CN_ABSYNAPSEECPLAST3_H

#include "queue.h"
#include "CN_absynapse.h"

#define ABECPLAST3PNO 8

#define ABECPLAST3IVARNO 1

double ABECPLAST3_INIVARS[ABECPLAST3IVARNO]= {
  0.0                    // 1 - S: internal variable
};

class absynapseECplast3: public absynapse
{  
 private:
  queue<double> chngTq;
  queue<double> dgq;
  double nextT;
  double nextdg;
  
 public:
  int synapse_change;
  absynapseECplast3(neuron *, neuron *, double, double, double, double,
			     double, double, double, double); 
  virtual ~absynapseECplast3();
  virtual void update_gsyn(double *);
  virtual double STDP_func(double);
};

#endif
