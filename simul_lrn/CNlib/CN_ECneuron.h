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


#ifndef DC_ECNEURON_H
#define DC_ECNEURON_H

#include "CN_neuron.h"
#include <cmath>

// parameters of the Entorhinal cortex stellate cell

#define ECN_IVARNO 6
#define ECN_PNO 13

double stdECN_p[ECN_PNO]= {
  7.15,          // 0 - gNa: Na conductance in 1/(mOhms * cm^2)
  50.0,          // 1 - ENa: Na equi potential in mV
  1.43,          // 2 - gK: K conductance in 1/(mOhms * cm^2)
  -95.0,         // 3 - EK: K equi potential in mV
  0.021,         // 4 - gl: leak conductance in 1/(mOhms * cm^2)
  -55.0,         // 5 - El: leak equi potential in mV
  0.035, //0.00572,  // 6 - gKl: potassium leakage conductivity
  -95.0,         // 7 - EKl: potassium leakage equi pot in mV
  65.0,          // 8 - V0: ~ total equi potential (?)
  0.286, //0.143,    // 9 - Cmem: membr. capacity density in muF/cm^2
  0.0185, //1.85,    // 10 - gh1
  0.01, // 1.0,      // 11 - gh2
  20.0           // 12 - Vh
};

double *ECN_p= stdECN_p;

char *ECN_p_text[ECN_PNO]= {
  "0 - gNa: Na conductance in 1/(mOhms * cm^2)",
  "1 - ENa: Na equi potential in mV",
  "2 - gK: K conductance in 1/(mOhms * cm^2)",
  "3 - EK: K equi potential in mV",
  "4 - gl: leak conductance in 1/(mOhms * cm^2)",
  "5 - El: leak equi potential in mV",
  "6 - gKl: potassium leakage conductivity",
  "7 - EKl: potassium leakage equi pot in mV",
  "8 - V0: ~ total equi potential (?)",
  "9 - Cmem: membr. capacity density in muF/cm^2",
  "10 - gh1",
  "11 - gh2",
  "12 - Vh"
};

double ECN_INIVARS[ECN_IVARNO]= {
  -64.1251,                       // 0 - membrane potential E
  0.0176331,                   // 1 - prob. for Na channel activation m
  0.994931,                   // 2 - prob. for not Na channel blocking h
  0.0433969,                   // 3 - prob. for K channel activation n
  0.443961,                         // 4 - Ih1 activation
  0.625308                         // 5 - Ih2 activation  
};

char *ECN_INIVARSTEXT[ECN_IVARNO]= {
  "0 - membrane potential E",
  "1 - prob. for Na channel activation m",
  "2 - prob. for not Na channel blocking h",
  "3 - prob. for K channel activation n",
  "4 - Ih1 activation",
  "5 - Ih2 activation"
};


// ECneuron class itself

class ECneuron: public neuron
{
 private:
  double Isyn;
  double _a, _b;
 public:
  ECneuron(int, double *);
  ECneuron(int, tnvector<int>, double *);
  ~ECneuron() { }
  inline virtual double E(double *);
  virtual void derivative(double *, double *);
};

#endif



