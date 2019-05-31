/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_ABSYNAPSE_CC
#define CN_ABSYNAPSE_CC

#include "CN_synapse.cc"

// This is the constructor to be used by derived classes passing the new
// internal var number, parameter number and type tag

absynapse::absynapse(neuron *insource, neuron *intarget,
		     double ingsyn, double inEsyn, double inEpre,
		     double inasyn, double inbsyn, double inVslope,
		     int inIVARNO, int inPNO, int inTYPE):
  synapse(insource, intarget, inIVARNO, inPNO, inTYPE)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= inasyn;           // alpha timescale in 1/msec
  p[4]= inbsyn;           // beta timescale in 1/msec
  p[5]= inVslope;         // steepness of activation curve as func of Vpre  
} 

// This is the constructor to be used directly ...

absynapse::absynapse(neuron *insource, neuron *intarget,
		     double ingsyn, double inEsyn, double inEpre,
		     double inasyn, double inbsyn, double inVslope):
  synapse(insource, intarget, ABSYNIVARNO, ABSYNPNO, ABSYN)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= inasyn;           // alpha timescale in 1/msec
  p[4]= inbsyn;           // beta timescale in 1/msec
  p[5]= inVslope;         // steepness of activation curve as func of Vpre  
} 

absynapse::absynapse(neuron *insource, neuron *intarget, double *inp):
  synapse(insource, intarget, ABSYNIVARNO, ABSYNPNO, ABSYN)
{
  set_p(inp);
} 

absynapse::~absynapse()
{
}

double absynapse::gsyn()
{
  return p[0];
}

void absynapse::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double absynapse::Isyn(double *x)
{
  return -p[0]*x[idx]*(target->E(x)-p[1]);
}


void absynapse::derivative(double *x, double *dx)
{
  dx[idx]= p[3]*(1.0-x[idx])*(1.0+tanh((source->E(x)-p[2])/p[5]))/2.0
    -p[4]*x[idx];
}

// end of class implementation

#endif


