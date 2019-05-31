/*--------------------------------------------------------------------------
   Author: Thomas Nowotny
  
   Institute: Institute for Nonlinear Dynamics
              University of California San Diego
              La Jolla, CA 92093-0402
  
   email to:  tnowotny@ucsd.edu
  
   initial version: 2002-02-14
  
--------------------------------------------------------------------------*/


#ifndef CN_INPUTNEURON_CC
#define CN_INPUTNEURON_CC

#include "CN_inputneuron.h"
#include "CN_neuron.cc"

inputneuron::inputneuron(int inlabel, int iniVarno, int intype,
			 tnvector<int> inpos, double *inp, int inpno):
  neuron(inlabel, iniVarno, intype, inpos, inp, inpno)
{
  tb= new double[1];
  te= new double[1];
  theE= p[2];
  t_last= -1.0;
  setIdx(-1);
}

inputneuron::~inputneuron()
{
  delete[] tb;
  delete[] te;
}

void inputneuron::init(double *x, double *iniVars)
{
  neuron::init(x,iniVars);
  theE= p[2];
  t_last= x[0]-1.0;
}

void inputneuron::set_input(int infno, double *intb,
			    double *inte, double instime= -1.0)
{
  if (instime > 0.0) p[4]= instime;
  fno= infno;             // no. of bursts in pattern
  if (fno > 0) {
    delete[] tb;
    delete[] te;
    tb= new double[fno]; 
    te= new double[fno];
    
    for (int i= 0; i < fno; i++)
    {
      intb[i]= intb[i] - ((int) ((intb[i])/p[4]))*p[4];
      tb[i]= intb[i];
      inte[i]= inte[i] - ((int) ((inte[i])/p[4]))*p[4];
      te[i]= inte[i];
    }
  }
}

#endif
