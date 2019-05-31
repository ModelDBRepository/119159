//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institute for Nonlinear Dynamics
//            University of California San Diego
//            La Jolla, CA 92093-0402
//
// email to:  tnowotny@ucsd.edu
//
// initial version: 2005-08-17
//
//--------------------------------------------------------------------------


#ifndef CN_TIMENEURON_H
#define CN_TIMENEURON_H

#include "CN_neuron.h"

#define TIME_IVARNO 1
#define TIME_PNO 0

double TIME_INIVARS[TIME_IVARNO]= {
  0.0                       // 0 - time
};

char *TIME_INIVARSTEXT[TIME_IVARNO]= {
  "0 - time"
};

// the time neuron class itself

class TimeNeuron: public neuron
{
 public:
  TimeNeuron();
  ~TimeNeuron() { }
  inline virtual double E(double *);
  void derivative(double *, double *);
};

#endif
