//--------------------------------------------------------------------------
// Author: Thomas Nowotny
//
// Institute: Institute for Nonlinear Dynamics
//            University of California San Diego
//            La Jolla, CA 92093-0402
//
// email to:  tnowotny@ucsd.edu
//
// initial version: 2005-08-19
//
//--------------------------------------------------------------------------

// this is an adaptation of Astrid Prinz's model of an LP neuron to two
// compartments and a few changed std parameters, done by Marcello Reyes

#ifndef CN_LPANEURON_H
#define CN_LPANEURON_H

#include "CN_neuron.h"
#include <cmath>

#define LPA_PNO 17
#define LPA_IVARNO 13

//--------------------------------------------------------------------------
// parameters of the LPA neuron (typical set, can be modified at any time
// later on. Modifications take effect immediately in tha next evaluation)

double stdLPA_p[LPA_PNO]= {
  1.0, // 0 - Capacit
  200.0, // 1 - tau_Ca
  1.4960, // 2 - f
  0.05, // 3 - C_Ca_0
  50.0, // 4 - E_Na
  -80.0, // 5 - E_K
  -20.0, // 6 - E_H
  -60.0, // 7 - E_leak
  0.01, // 8 - g_H
  0.01, // 9 - g_leak
  120.0, // 10 - g_Na
  1.0, // 11 - g_CaT
  11.0, // 12 - g_CaS
  50.0, // 13 - gA
  1.22, // 14 - gKCa
  36.0, // 15 - gKd
  0.628 // 16 - Area
};

double *LPA_p= stdLPA_p;

char *LPA_p_text[LPA_PNO]= {
  "0 - Capacit",
  "1 - tau_Ca",
  "2 - f",
  "3 - C_Ca_0",
  "4 - E_Na",
  "5 - E_K",
  "6 - E_H",
  "7 - E_leak",
  "8 -  g_H",
  "9 - g_leak",
  "10 - g_Na",
  "11 - g_CaT",
  "12 - g_CaS",
  "13 - gA",
  "14 - gKCa",
  "15 - gKd",
  "16 - Area"
};
  

//----------------------------------------------------------------------
// the following initial values reflect the steady state with no input

double stdLPA_INIVARS[LPA_IVARNO]= {
  -52.6896657132625989600,   // 0 - V: membrane potential 
  0.0056935725624609398,     // 1 - mNa: Na channel activation 
  0.9,                       // 2 - hNa: Na channel unblocking 
  0.6836680498672979000,     // 3 - mCa1: Ca1 channel activation
  0.0260241567723770370,     // 4 - hCa1: Ca1 channel unblocking
  0.9842652361470087818,     // 5 - mCa2: Ca2 channel activation
  0.0708515410780189148,     // 6 - hCa2: Ca2 channel unblocking
  0.3398902413785575005,     // 7 - mA: IA channel activation
  0.0474145467346798327,     // 8 - hA: IA channel unblocking
  0.3654065423783397493,     // 9 - mKCa: KCa channel activation
  0.0526824330148126588,     // 10 - mKd: Kd channel activation
  0.0302975756966933421,     // 11 - mh: Ih channel activation
  2.1790550738448457580      // 12 - CCa: Ca concentration
};
 
double *LPA_INIVARS= stdLPA_INIVARS;

char *LPA_INIVARSTEXT[LPA_IVARNO]= {
  "0 - V: membrane potential",
  "1 - mNa: Na channel activation", 
  "2 - hNa: Na channel unblocking", 
  "3 - mCa1: Ca1 channel activation",
  "4 - HCa1: Ca1 channel unblocking",
  "5 - mCa2: Ca2 channel activation",
  "6 - hCa2: Ca2 channel unblocking",
  "7 - mA: IA channel activation",
  "8 - hA: IA channel unblocking",
  "9 - mKCa: KCa channel activation",
  "10 - mKd: Kd channel activation",
  "11 - mh: Ih channel activation",
  "12 - CCa: Ca concentration"
};

// the LPA neuron class itself

class LPAneuron: public neuron
{
 private:
  double Isyn, I_Na, I_CaT, I_CaS, I_A, I_KCa, I_Kd, I_H, I_leak;
  double _E_Ca, _minf, _taum, _hinf, _tauh;

 public:
  LPAneuron(int, double *);
  LPAneuron(int, tnvector<int>, double *);
  virtual ~LPAneuron() { }
  inline virtual double E(double *);
  virtual void currents(ostream &, double *);
  virtual void derivative(double *, double *);
};

// for fitting purposes: type == 0 means additive changes
//                       type == 1 means multiplicative changes

double LPA_p_type[LPA_PNO]= {
  1, // 0 - Capacit
  1, // 1 - tau_Ca
  1, // 2 - f
  1, // 3 - C_Ca_0
  0, // 4 - E_Na
  0, // 5 - E_K
  0, // 6 - E_H
  0, // 7 - E_leak
  1, // 8 - g_H
  1, // 9 - g_leak
  1, // 10 - g_Na
  1, // 11 - g_CaT
  1, // 12 - g_CaS
  1, // 13 - gA
  1, // 14 - gKCa
  1, // 15 - gKd
  1  // 16 - Area
};

// the sensitivity is taken from LPneuronMarcello ... might need to be
// adjusted

double LPA_p_sens[LPA_PNO]= {
  3.64202e-10,
  5.76844e-09,
  8.09137e-10,
  9.52267e-07,
  2.52966e-08,
  6.75621e-09,
  3.75613e-06,
  3.23909e-08,
  9.21422e-08,
  2.01029e-09,
  3.99026e-10,
  1.49669e-10,
  4.76687e-10,
  2.43645e-10,
  4.44189e-10,
  4.13484e-10,
  4.0e-10
};


#endif

