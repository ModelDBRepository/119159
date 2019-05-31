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


#ifndef CN_HHNEURON_H
#define CN_HHNEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define HH_IVARNO 4
#define HH_PNO 7

double stdHH_p[HH_PNO]= {
  120.0,         // 0 - gNa: Na conductance in 1/(mOhms * cm^2)
  55.0,          // 1 - ENa: Na equi potential in mV
  36.0,          // 2 - gK: K conductance in 1/(mOhms * cm^2)
  -72.0,         // 3 - EK: K equi potential in mV
  0.3,           // 4 - gl: leak conductance in 1/(mOhms * cm^2)
  -50.0,         // 5 - El: leak equi potential in mV
  1.0          // 6 - Cmem: membr. capacity density in muF/cm^2
};

double *HH_p= stdHH_p;

char *HH_p_text[HH_PNO]= {
  "0 - gNa: Na conductance in 1/(mOhms * cm^2)",
  "1 - ENa: Na equi potential in mV",
  "2 - gK: K conductance in 1/(mOhms * cm^2)",
  "3 - EK: K equi potential in mV",
  "4 - gl: leak conductance in 1/(mOhms * cm^2)",
  "5 - El: leak equi potential in mV",
  "6 - Cmem: membr. capacity density in muF/cm^2"
};

double HH_INIVARS[HH_IVARNO]= {
  -60.0,                       // 0 - membrane potential E
  0.0529324,                   // 1 - prob. for Na channel activation m
  0.3176767,                   // 2 - prob. for not Na channel blocking h
  0.5961207                    // 3 - prob. for K channel activation n
};

char *HH_INIVARSTEXT[HH_IVARNO]= {
  "0 - membrane potential E",
  "1 - prob. for Na channel activation m",
  "2 - prob. for not Na channel blocking h",
  "3 - prob. for K channel activation n"
};

// the HH neuron class itself

class HHneuron: public neuron
{
 private:
  double Isyn;
  double _a, _b;
 public:
  HHneuron(int, double *);
  HHneuron(int, tnvector<int>, double *);
  ~HHneuron() { }
  inline virtual double E(double *);
  void derivative(double *, double *);
};

#endif
