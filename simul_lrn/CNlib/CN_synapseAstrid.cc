/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-19
  
--------------------------------------------------------------------------*/

#ifndef CN_SYNAPSEASTRID_CC
#define CN_SYNAPSEASTRID_CC

#include "CN_synapse.cc"

// This is the constructor to be used by derived classes passing the new
// internal var number, parameter number and type tag

synapseAstrid::synapseAstrid(neuron *insource, neuron *intarget,
		     double ingsyn, double inEsyn, double inEpre,
		     double indelta, double inkminus,
		     int inIVARNO, int inPNO, int inTYPE):
  synapse(insource, intarget, inIVARNO, inPNO, inTYPE)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= indelta;          // delta Voltage in mV
  p[4]= inkminus;         // kminus rate constant
} 

// This is the constructor to be used directly ...

synapseAstrid::synapseAstrid(neuron *insource, neuron *intarget,
		     double ingsyn, double inEsyn, double inEpre,
		     double indelta, double inkminus):
  synapse(insource, intarget, SYNASIVARNO, SYNASPNO, SYNAS)
{
  p[0]= ingsyn;           // gsyn strength of synapse
  p[1]= inEsyn;           // Esyn reversal potential in mV
  p[2]= inEpre;           // Epre presyn threshold potential in mV
  p[3]= indelta;          // delta Voltage in mV
  p[4]= inkminus;         // kminus rate constant
} 

synapseAstrid::synapseAstrid(neuron *insource, neuron *intarget, double *inp):
  synapse(insource, intarget, SYNASIVARNO, SYNASPNO, SYNAS)
{
  set_p(inp);
} 

synapseAstrid::~synapseAstrid()
{
}

double synapseAstrid::gsyn()
{
  return p[0];
}

void synapseAstrid::set_gsyn(double ingsyn)
{
  p[0]= ingsyn;
}

double synapseAstrid::Isyn(double *x)
{
  return -p[0]*x[idx]*(target->E(x)-p[1]);
}


// differential eqn for f, the driving variable for the pre-synaptic
// conductance g and g itself

void synapseAstrid::derivative(double *x, double *dx)
{
  sinf= 1.0/(1.0+exp((p[2]-source->E(x))/p[3]));
  tau= (1.0-sinf)/p[4];
  if (tau < 1e-2) tau= 1e-2;
  dx[idx]= (sinf-x[idx])/tau;
}

// end of class implementation

#endif


