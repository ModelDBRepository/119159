/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2005-08-17
  
--------------------------------------------------------------------------*/

#ifndef CN_NEURONMODEL
#define CN_NEURONMODEL

#include "CN_NeuronModel.h"
#include "CN_TimeNeuron.h"
#include "CN_TimeNeuron.cc"

NeuronModel::NeuronModel(slist<neuron *> *neurs, slist<synapse *> *syns,
			 int &N, ostream &msgos)
{
  int iVarCnt= 0;
  neuron *n;
  synapse *s;
  
  niter= neurs->iterator();
  siter= syns->iterator();
  n= new TimeNeuron();
  neurs->prepend(n);

  // count all variables
  // & deliver the variable positions to everyone
  forall(niter) {
    n= niter->c_value();
    n->setIdx(iVarCnt);
    iVarCnt+= n->iVarNo;
  }
  forall(siter) {
    s= siter->c_value();
    s->setIdx(iVarCnt);
    iVarCnt+= s->iVarNo;
  }
  N= iVarCnt;
  msgos << "# we have " << N << " variables ..." << endl;
}

void NeuronModel::derivative(double *x, double *dx)
{
  forall(niter) {
    niter->c_value()->derivative(x, dx);
  }
  forall(siter) {
    siter->c_value()->derivative(x, dx);
  }
}

#endif
