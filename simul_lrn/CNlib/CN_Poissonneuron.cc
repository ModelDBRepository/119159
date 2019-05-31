/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/


#ifndef CN_POISSONNEURON_CC
#define CN_POISSONNEURON_CC

#include "CN_neuron.cc"

Poissonneuron::Poissonneuron(int inlabel, tnvector<int> inpos,
			     double *the_p= POI_p):
  neuron(inlabel, POI_IVARNO, POISSONNEURON,
	      inpos, the_p, POI_PNO)
{
  firing= 0;
  refract= 0;
  fire_t= 0.0;
  myx= new double[1];
  myxn= new double[1];
  myx[0]= p[2];
  myxn[0]= p[2];
  setIdx(-1);
}

Poissonneuron::~Poissonneuron()
{
  delete[] myx;
  delete[] myxn;
}

double Poissonneuron::E(double *x)
{
  return myx[0];
}

void Poissonneuron::validate_E(double *x, double dt)
{
  if (firing) {
    if (x[0] - fire_t > p[0]) { // remember: x[0] is the time
      firing= 0;
      refract= 1;
    }
  }
  else {
    if (refract) {
      if (x[0] - fire_t > p[1]) refract=0;
    }
    else {
      if (R.n() <= p[4]*dt) {
	firing= 1;
	fire_t= x[0];
      }
    }
  }
  if (firing) myxn[0]= p[3];
  else myxn[0]= p[2];
}

void Poissonneuron::step()
{
  myx[0]= myxn[0];
}

void Poissonneuron::init(double *x, double *iniVars)
{
  myx[0]= iniVars[0];
  myxn[0]= iniVars[0];
}

#endif





