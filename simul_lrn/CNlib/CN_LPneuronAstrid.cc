/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-19
  
--------------------------------------------------------------------------*/

#ifndef CN_LPANEURON_CC
#define CN_LPANEURON_CC

#include "CN_neuron.cc"

LPAneuron::LPAneuron(int inlabel, double *the_p= LPA_p):
  neuron(inlabel, LPA_IVARNO, LPANEURON, the_p, LPA_PNO)
{
}

LPAneuron::LPAneuron(int inlabel, tnvector<int> inpos, double *the_p= LPA_p):
  neuron(inlabel, LPA_IVARNO, LPANEURON, inpos, the_p, LPA_PNO)
{
}

inline double LPAneuron::E(double *x)
{
  return x[idx];
}

#define efunc(X,Y,Z) (1.0/(1.0+expl(((X)+(Y))/(Z))))
#define E x[idx]

#define Capacit p[0]
#define tau_Ca p[1]
#define f p[2]
#define C_Ca_0 p[3]
#define E_Na p[4]
#define E_K p[5]
#define E_H p[6]
#define E_leak p[7]
#define g_H p[8]
#define g_leak p[9]
#define g_Na p[10]
#define g_CaT p[11]
#define g_CaS p[12]
#define g_A p[13]
#define g_KCa p[14]
#define g_Kd p[15]
#define Area p[16]

void LPAneuron::currents(ostream &os, double *x)
{
  _E_Ca=(12.5*log(3000.0/x[idx+12]));
  
  I_Na=  g_Na*ipower(x[idx+1],3)*x[idx+2]*(E-E_Na)*Area;
  I_CaT= g_CaT*ipower(x[idx+3],3)*x[idx+4]*(E-_E_Ca)*Area;
  I_CaS= g_CaS*ipower(x[idx+5],3)*x[idx+6]*(E-_E_Ca)*Area;
  I_A=   g_A*ipower(x[idx+7],3)*x[idx+8]*(E-E_K)*Area;
  I_KCa= g_KCa*ipower(x[idx+9],4)*(E-E_K)*Area;
  I_Kd=  g_Kd*ipower(x[idx+10],4)*(E-E_K)*Area;
  I_H=   g_H*x[idx+11]*(E-E_H)*Area;
  I_leak= g_leak*(E-E_leak)*Area;
  
  os << I_Na << " ";
  os << I_CaT << " ";
  os << I_CaS << " ";
  os << I_A << " ";
  os << I_KCa << " ";
  os << I_Kd << " ";
  os << I_H << " ";
  os << I_leak << " ";
}

void LPAneuron::derivative(double *x, double *dx)
{
  Isyn= 0.0;
  forall (den_it) {
    Isyn+= den_it->c_value()->Isyn(x);
  }

  _E_Ca=(12.5*log(3000.0/x[idx+12]));
  
  I_Na=  g_Na*ipower(x[idx+1],3)*x[idx+2]*(E-E_Na);
  I_CaT= g_CaT*ipower(x[idx+3],3)*x[idx+4]*(E-_E_Ca);
  I_CaS= g_CaS*ipower(x[idx+5],3)*x[idx+6]*(E-_E_Ca);
  I_A=   g_A*ipower(x[idx+7],3)*x[idx+8]*(E-E_K);
  I_KCa= g_KCa*ipower(x[idx+9],4)*(E-E_K);
  I_Kd=  g_Kd*ipower(x[idx+10],4)*(E-E_K);
  I_H=   g_H*x[idx+11]*(E-E_H);
  I_leak= g_leak*(E-E_leak);
  
  dx[idx]= Isyn/Capacit -(I_Na + I_Kd + I_CaT + I_CaS + I_A + I_KCa +
			       I_H + I_leak)*Area/Capacit;
  
  // differential eqn for mNa
  _minf=efunc(E,25.5,-5.29);
  _taum=2.64-2.52*efunc(E,120.0,-25.0);
  dx[idx+1]=(_minf-x[idx+1])/_taum;
  // differential eqn for hNa
  _hinf= efunc(E,48.9,5.18);
  _tauh= (1.34*efunc(E,62.9,-10.0))*(1.5+efunc(E,34.9,3.6));
  if (_tauh < 1e-2) _tauh= 1e-2; // this is a horrible numerical instability
  // if the neuron is hyperpolarized, _tauh goes to 0 amplifying the already
  // clear disadvantageous rounding errors in (order of 1 - order of 1)
  // in mNa *a lot*. This leads to an explosion and eventually nan ...
  //  cerr << _hinf << " " << _tauh << " ";
  dx[idx+2]=(_hinf-x[idx+2])/_tauh;

  // differential eqn for mCa1
  _minf= efunc(E,27.1,-7.2);
  _taum= 43.4-42.6*efunc(E,68.1,-20.5);
  dx[idx+3]=(_minf-x[idx+3])/_taum;
  // differential eqn for hCa1
  _hinf= efunc(E,32.1,5.5);
  _tauh= 210.0-179.6*efunc(E,55.0,-16.9);
  dx[idx+4]=(_hinf-x[idx+4])/_tauh;

  // differential eqn for mCa2
  _minf= efunc(E,33.0,-8.1);
  _taum= 2.8+(14.0/(exp((E+27.0)/10.0)+exp((E+70.0)/-13.0)));
  dx[idx+5]=(_minf-x[idx+5])/_taum;
  // differential eqn for hCa2
  _hinf= efunc(E,60,6.2);
  _tauh= 120+(300.0/(exp((E+55.0)/9.0)+exp((E+65.0)/-16)));
  dx[idx+6]=(_hinf-x[idx+6])/_tauh;

  // differential eqn for mA
  _minf= efunc(E,27.2,-8.7);
  _taum= 23.2-20.8*efunc(E,32.9,-15.2);
  dx[idx+7]=(_minf-x[idx+7])/_taum;
  // differential eqn for hA
  _hinf= efunc(E,56.9,4.9);
  _tauh= 77.2-58.4*efunc(E,38.9,-26.5);
  dx[idx+8]=(_hinf-x[idx+8])/_tauh;

  // differential eqn for mKCa
  _minf= (x[idx+12]/(x[idx+12]+3.0))*efunc(E,28.3,-12.6);
  _taum= 180.6-150.2*efunc(E,46.0,-22.7);
  dx[idx+9]=(_minf-x[idx+9])/_taum;
  
  // differential eqn for mKd
  _minf= efunc(E,12.3,-11.8);
  _taum= 14.4-12.8*efunc(E,28.3,-19.2);
  dx[idx+10]=(_minf-x[idx+10])/_taum;
  
  // differential eqn for mh
  _minf= efunc(E,75.0,5.5);
  _taum= 2.0/(exp((E+169.7)/-11.6)+exp((E-26.7)/14.3));
  if (_taum < 1e-2) _taum= 1e-2; // same instability as for _tauh_Na ...
  dx[idx+11]=(_minf-x[idx+11])/_taum;

  // differential eqn for o, prob. of Ca mediated K channel activation 
  dx[idx+12]=(-f*(I_CaT+I_CaS)-x[idx+12]+C_Ca_0)/tau_Ca;

  // the Isyn is calculated elsewhere
  dx[idx+13]= 0.0; 
}

#undef efunc
#undef E
#undef Capacit
#undef tau_Ca
#undef f
#undef C_Ca_0
#undef E_Na
#undef E_K
#undef E_H
#undef E_leak
#undef g_H
#undef g_leak
#undef g_Na
#undef g_CaT
#undef g_CaS
#undef g_A
#undef g_KCa
#undef g_Kd
#undef Area

#endif
