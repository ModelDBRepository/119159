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


#ifndef CN_VALNEURON_H
#define CN_VALNEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the HH neuron, they are identical for all neurons used
// (and therefore made global to save memory)

#define Val_IVARNO 4
#define Val_PNO 10

double stdVal_p[Val_PNO]= {
  7.15,          // 0 - gNa: Na conductance in 1/(mOhms * cm^2)
  50.0,          // 1 - ENa: Na equi potential in mV
  1.43,          // 2 - gK: K conductance in 1/(mOhms * cm^2)
  -95.0,         // 3 - EK: K equi potential in mV
  0.021,         // 4 - gl: leak conductance in 1/(mOhms * cm^2)
  -55.0,         // 5 - El: leak equi potential in mV
  0.00572,       // 6 - gKl: potassium leakage conductivity
  -95.0,         // 7 - EKl: potassium leakage equi pot in mV
  65.0,          // 8 - V0: ~ total equi potential (?)
  0.143        // 9 - Cmem: membr. capacity density in muF/cm^2
};

double *Val_p= stdVal_p;

char *Val_p_text[Val_PNO]= {
  "0 - gNa: Na conductance in 1/(mOhms * cm^2)",
  "1 - ENa: Na equi potential in mV",
  "2 - gK: K conductance in 1/(mOhms * cm^2)",
  "3 - EK: K equi potential in mV",
  "4 - gl: leak conductance in 1/(mOhms * cm^2)",
  "5 - El: leak equi potential in mV",
  "6 - gKl: potassium leakage conductivity",
  "7 - EKl: potassium leakage equi pot in mV",
  "8 - V0: ~ total equi potential (?)",
  "9 - Cmem: membr. capacity density in muF/cm^2"
};

double Val_INIVARS[Val_IVARNO]= {
  -60.0,                       // 0 - membrane potential E
  0.0529324,                   // 1 - prob. for Na channel activation m
  0.3176767,                   // 2 - prob. for not Na channel blocking h
  0.5961207                    // 3 - prob. for K channel activation n
};

char *Val_INIVARSTEXT[Val_IVARNO]= {
  "0 - membrane potential E",
  "1 - prob. for Na channel activation m",
  "2 - prob. for not Na channel blocking h",
  "3 - prob. for K channel activation n"
};


// Valentins HH neuron class itself

class Valneuron: public neuron
{
 private:
  double Isyn;
  double _a, _b;
 public:
  Valneuron(int, double *);
  Valneuron(int, tnvector<int>, double *);
  ~Valneuron() { }
  inline virtual double E(double *);
  virtual void derivative(double *, double *);
};

#endif



